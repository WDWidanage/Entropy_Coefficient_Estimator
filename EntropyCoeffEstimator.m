classdef EntropyCoeffEstimator < handle
    % Class to design multilevel temperature signals and estimate the
    % entropy coefficient of a li-ion battery reversible heat source term.
    %
    % Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K.
    %               Oliver Queisser - TVT, Karlsruhe Institute of
    %               Technology, Germany
    %               Sabine Paarmann - TVT, Karlsruhe Institute of
    %               Technology, Germany
    % All Rights Reserved
    % Software may be used freely for non-comercial purposes only
    % 04/10/2021 (Quite)

    properties
        estimationSettings struct = struct();
        cellThermalProperties struct = struct();
    end

    properties (SetAccess = private)
        refSig struct = struct();
        measData struct = struct();
        processedData struct =  struct();
        results struct =  struct();
    end

    methods
        % Constructor method
        function obj = EntropyCoeffEstimator(nvp)

            arguments
                % Signal desgin related properties
                nvp.Tmin_degC(1,1) double {mustBeNumeric} = 10;         % Specify minimum temperature [degC]
                nvp.Tmax_degC(1,1) double {mustBeNumeric} = 50;         % Specify maximum temperature [degC]
                nvp.Tp_H(1,1) double {mustBePositive} = 4;              % Specify signal period [H]
                nvp.timeStep_mins(1,1) double {mustBePositive} = 5;     % Each temperature step will at least be held for timeStep minutes before switching [mins]
                nvp.Ts_s(1,1) double {mustBeNumeric} = 1;               % Sampling rate of data acquisiton [s]
                nvp.Hsupp double = 2;                                   % Suppressed harmonics
                nvp.numLevels = 5;                                      % Sepcify the number of temperature levels
                nvp.startTemp = 20;                                     % Specify the starting temperature [degC]
                nvp.fMax_Hz(1,1) double {mustBeNumeric} = 1.5E-3;       % Maximum desired frequency [Hz]
                nvp.Periods(1,1) double {mustBeInteger} = 1;            % Default number of periods
                nvp.fileName string = "fileName.csv";                   % CSV file name to save P eriods of the reference signal as two columns, time [s] and reference temperature [degC]

                % Cell thermal propeties
                nvp.cellConductivity_wpmk (1,1) double {mustBeNumeric} = 0.86;           % Cell thermal conductance [W/m/K]
                nvp.cellThickness_m (1,1) double {mustBeNumeric} = 7.5E-3;               % Cell thickness [m]
                nvp.cellSpecificHeatCapacity_Jp (1,1) double {mustBeNumeric} = 789.654;  % Cell specific heat capacity [J/kg/K]
                nvp.cellDenisty_kgpm3 (1,1) double {mustBeNumeric} = 2893.916;           % Cell density [kg/m^3]
                nvp.cellVolume_m3 (1,1) double {mustBeNumeric}= 3.302E-5;                % Cell volume [m^3]
                nvp.copperConductivity_wpmk (1,1) double {mustBeNumeric}= 398.71;        % Copper plate conductivity [W/m/K]
                nvp.copperThickness_m (1,1) double {mustBeNumeric}= 3E-3;                % Copper plate thickness [m]
                nvp.GeoFactor (1,1) double {mustBeInteger} = 2;                          % Geometric factor. 2 for plate

                % Measured data
                nvp.ocvData timetable = timetable.empty()               % Place holder for measured OCV data
                nvp.tempData timetable = timetable.empty()              % Place holder for measured temperature data

            end

            % Populate reference signal fields
            obj.refSig.Tmin_degC = nvp.Tmin_degC;
            obj.refSig.Tmax_degC = nvp.Tmax_degC;
            obj.refSig.Tp_H = nvp.Tp_H;
            obj.refSig.timeStep_mins = nvp.timeStep_mins;
            obj.refSig.Ts_s = nvp.Ts_s;
            obj.refSig.Hsupp = nvp.Hsupp;
            obj.refSig.fMax_Hz = nvp.fMax_Hz;
            obj.refSig.startTemp = nvp.startTemp;
            obj.refSig.numLevels = nvp.numLevels;
            obj.refSig.numPeriods = nvp.Periods;
            obj.refSig.textFileName = nvp.fileName;

            % Populate cell thermal propeties
            obj.cellThermalProperties.cellConductivity_wpmk = nvp.cellConductivity_wpmk;
            obj.cellThermalProperties.cellThickness_m = nvp.cellThickness_m;
            obj.cellThermalProperties.cellSpecificHeatCapacity_Jp = nvp.cellSpecificHeatCapacity_Jp;
            obj.cellThermalProperties.cellDenisty_kgpm3 = nvp.cellDenisty_kgpm3;
            obj.cellThermalProperties.cellVolume_m3 = nvp.cellVolume_m3;
            obj.cellThermalProperties.copperConductivity_wpmk = nvp.copperConductivity_wpmk;
            obj.cellThermalProperties.copperThickness_m = nvp.copperThickness_m;
            obj.cellThermalProperties.GeoFactor = nvp.GeoFactor;

            % Populate measured data as null
            obj.measData.ocvData = nvp.ocvData;
            obj.measData.tempData = nvp.tempData;

            obj.results.cntr = 0;

            % Create mltilevel temperature profile
            obj = obj.MultiLevTempProfile();

            % Create P periods and save a text file for Peltier
            % Create two periods for CSV file
            P = obj.refSig.numPeriods;
            refSigTmp = repmat(obj.refSig.refTempSig,P,1);
            timeVecTmp = [1:length(refSigTmp)]'*obj.refSig.Ts_s - obj.refSig.Ts_s;
            fileName = obj.refSig.textFileName;
            idx0 = [1;find([0;diff(refSigTmp)]~=0)];
            refSig = [refSigTmp(idx0);refSigTmp(end)];
            timeVec = [timeVecTmp(idx0);timeVecTmp(end)];
            csvwrite(fileName,[timeVec,refSig]);
            fileNameNoExt = extractBefore(fileName,".csv");
            save(fileNameNoExt +".mat","obj")

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Import a legacy reference signal
        % (This function is not required if the reference signal is
        % generated with via the EntropyCoeffEstimator class)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ImportLegacyRefSig(obj,nvp)
            arguments
                obj EntropyCoeffEstimator
                nvp.filePth string = [];
            end
            if isempty(nvp.filePth)
                [ocvFile,OCVFilePth] = uigetfile({'*.mat';},...
                    'Select legacy reference signal MAT file');
                fullFilePth = fullfile(OCVFilePth,ocvFile);
            else
                fullFilePth = nvp.filePth;
            end
            legacyRefSig = load(fullFilePth);
            UpdateRefSig(obj,legacyRefSig);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add measured data functionality
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ImportRefSig(obj,nvp)
            arguments
                obj EntropyCoeffEstimator
                nvp.filePth string = [];
            end
            if isempty(nvp.filePth)
                [ocvFile,OCVFilePth] = uigetfile({'*.mat';},...
                    'Select reference signal MAT file');
                fullFilePth = fullfile(OCVFilePth,ocvFile);
            else
                fullFilePth = nvp.filePth;
            end
            refSignalObj = load(fullFilePth);
            refSigName = fieldnames(refSignalObj);
            obj.refSig = refSignalObj.(refSigName{1}).refSig;
        end

        function ImportOCVData(obj,nvp)

            arguments
                obj EntropyCoeffEstimator
                nvp.filePth string = [];
                nvp.HdrLine = 1;
                nvp.HdrNames = ["Time_s","OCV_V"];
                nvp.DecimalSep = ",";
            end

            if isempty(nvp.filePth)
                [ocvFile,OCVFilePth] = uigetfile({'*.txt';'*.csv';'*.mat';},...
                    'Select OCV file');
                fullFilePth = fullfile(OCVFilePth,ocvFile);
            else
                fullFilePth = nvp.filePth;
            end
            imOpts = detectImportOptions(fullFilePth,'VariableNamesLine',nvp.HdrLine);


            % Specify column names and types
            imOpts.VariableNames =  nvp.HdrNames;
            imOpts.VariableTypes = repmat("double",size(imOpts.VariableNames));

            % Specify file level properties
            imOpts.ExtraColumnsRule = "ignore";
            imOpts.EmptyLineRule = "read";

            % Specify variable properties
            imOpts = setvaropts(imOpts, imOpts.VariableNames, "DecimalSeparator", nvp.DecimalSep);

            ocvData = readtable(fullFilePth, imOpts);
            ocvData.Time_s = seconds(ocvData.(1) - ocvData.(1)(1)); % Rest time to 0 seconds
            obj.measData.ocvData = table2timetable(ocvData);

        end

        function ImportTemperatureData(obj,nvp)
            arguments
                obj EntropyCoeffEstimator
                nvp.filePth string = [];
                nvp.HdrLine = 1;
                nvp.HdrNames = strings();
                nvp.DecimalSep = ",";
            end

            if isempty(nvp.filePth)
                [ocvFile,OCVFilePth] = uigetfile({'*.txt';'*.csv';'*.mat';},...
                    'Select temperature file');
                fullFilePth = fullfile(OCVFilePth,ocvFile);
            else
                fullFilePth = nvp.filePth;
            end
            imOpts = detectImportOptions(fullFilePth,'VariableNamesLine',nvp.HdrLine);

            % Specify column names and types
            if ~isempty(nvp.HdrNames)
                imOpts.VariableNames = nvp.HdrNames;
            end
            imOpts.VariableTypes = repmat("double",size(imOpts.VariableNames));

            % Specify file level properties
            imOpts.ExtraColumnsRule = "ignore";
            imOpts.EmptyLineRule = "read";

            % Specify variable properties
            imOpts = setvaropts(imOpts, imOpts.VariableNames, "DecimalSeparator", nvp.DecimalSep);

            tempData = readtable(fullFilePth, imOpts);
            tempData.Time_s = seconds(tempData.(1) - tempData.(1)(1));% Rest time to 0 seconds
            obj.measData.tempData = table2timetable(tempData);

            obj.CaloricCellTemperature();

        end

        function ImportExpData(obj,nvp)
            % Added to handle single text file with OCV and temperature readings 
            arguments
                obj EntropyCoeffEstimator
                nvp.filePth string = [];
                nvp.HdrLine = 2;
                nvp.HdrNames = strings();
                nvp.DecimalSep = ".";
            end

            if isempty(nvp.filePth)
                [ocvFile,OCVFilePth] = uigetfile({'*.txt';'*.csv';'*.mat';},...
                    'Select temperature file');
                fullFilePth = fullfile(OCVFilePth,ocvFile);
            else
                fullFilePth = nvp.filePth;
            end
            imOpts = detectImportOptions(fullFilePth,'VariableNamesLine',nvp.HdrLine);

            % Specify column names and types
            if ~isempty(nvp.HdrNames)
                imOpts.VariableNames = nvp.HdrNames;
            end
            imOpts.VariableTypes = repmat("double",size(imOpts.VariableNames));

            % Specify file level properties
            imOpts.ExtraColumnsRule = "ignore";
            imOpts.EmptyLineRule = "read";

            % Specify variable properties
            imOpts = setvaropts(imOpts, imOpts.VariableNames, "DecimalSeparator", nvp.DecimalSep);

            tempData = readtable(fullFilePth, imOpts);
            Time_s = seconds(tempData.(1) - tempData.(1)(1)); % Reset time to 0 seconds
            tempData.Time_s = Time_s;
            obj.measData.tempData = table2timetable(tempData);

            % Call interpolation function
            obj.dataInterPolate();

            obj.CaloricCellTemperature();
            
%             obj.measData.ocvData = timetable(Time_s,tempData.U,'VariableNames',{'OCV_V'});
            obj.measData.ocvData = obj.measData.tempData(:,end);
            obj.measData.ocvData.Properties.VariableNames = {'OCV_V'};

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimate dUdT. Estimate non-parameteric FRF and parameterise to
        % determine dudT (the steady-state gain value)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = EstimateEntropyCoeff(obj,nvp)  % Estimate frequency response (kernel function) and parameterise it with a transfer function in Fourier space
            arguments
                obj
                % dUdT estimation settings
                nvp.transientOnOff (1,1) {mustBeTextScalar} = " ";       % Set to "on" if OCV hasn't reached steady state
                nvp.modelOrder_num = [];                                % Transfer function numerator model order
                nvp.modelOrder_denom = [];                              % Transfer function denomenator model order
                nvp.usePeriods = [];                                    % Use the following periods for kernel function fit
            end

            % Populate estimation settings
            obj.estimationSettings.transientOnOff = nvp.transientOnOff;
            obj.estimationSettings.modelOrder_num = nvp.modelOrder_num;
            obj.estimationSettings.modelOrder_denom = nvp.modelOrder_denom;
            obj.estimationSettings.usePeriods = nvp.usePeriods;

            obj.ReshapeAveragePeriods();
            obj.estimateFRF();
            obj.fitFRF()
        end

        %%%%%%%%%%%%%%%%%%%%
        % Ploting functions
        %%%%%%%%%%%%%%%%%%%%
        function PlotRefSig(obj)
            % Frequency domain
            refHarm = 1:floor(obj.refSig.Nref/2);
            refFreqVec = refHarm/(obj.refSig.Nref*obj.refSig.Ts_s);
            UTmp = fft(obj.refSig.refTempSig);

            % Excited components
            UrefHarm = UTmp(refHarm+1);
            UrefExc = UTmp(obj.refSig.excHarmonics+1);

            % Suppressed componenets
            suppHarmMul = obj.refSig.Hsupp; % Suppressed harmonic multiples
            suppHarmsTmp  = [];
            for hh = 1:length(suppHarmMul)
                suppHarmsTmp = [suppHarmsTmp ,suppHarmMul(hh):suppHarmMul(hh): obj.refSig.excHarmonics(end)];
            end
            suppHarms = sort(suppHarmsTmp); % Sort the suppressed harmonics in ascending order
            suppFreq_Hz = suppHarms/(obj.refSig.Nref*obj.refSig.Ts_s);
            UrefSupp = UTmp(suppHarms+1);

            figure()
            subplot(2,1,1)
            stairs(obj.refSig.refTimeVec_s/3600,obj.refSig.refTempSig,'. -');
            xlabel('Time (H)'); ylabel({'Reference temperature';'Signal [deg^\circC]'})

            subplot(2,1,2)
            semilogx(obj.refSig.excFreq_Hz*1000,abs(UrefExc),'o r',refFreqVec*1000,abs(UrefHarm),'* b'); hold on;
            xlabel('Freq (mHz)'); ylabel({' Temperature FFT'; 'magnitude [-]'})
            legend('Excited frequencies','All frequencies')
            if ~isempty(suppHarms)
                semilogx(suppFreq_Hz*1000,abs(UrefSupp),'o g'); hold on;
                legend('Excited frequencies','All frequencies', 'Suppressed frequencies')
            end
        end

        function PlotMeasSig(obj)
            % Plot all measured temperature signals and OCV agaist time
            figure
            stackedplot(hours(obj.measData.tempData.Time_s),obj.measData.tempData.Variables,"DisplayLabels",obj.measData.tempData.Properties.VariableNames);
            xlabel('Time [s]')
            figure
            stackedplot(hours(obj.measData.ocvData.Time_s),obj.measData.ocvData.Variables,"DisplayLabels",obj.measData.ocvData.Properties.VariableNames);
            xlabel('Time [s]')
        end

        function PlotProcessedSig(obj)

            obj.ReshapeAveragePeriods();

            % Frequency domain
            refHarm = 1:floor(obj.refSig.Nref/2);
            refFreqVec = refHarm/(obj.refSig.Nref*obj.refSig.Ts_s);
            URefTmp = fft(obj.refSig.refTempSig);
            UMeasTmp = fft(obj.processedData.CaloricPeriod_degC.Caloric_Period);
            OMeasTmp = fft(obj.processedData.ocvPeriod_V.OCV_Period);

            % Reference signal excited components
            UrefHarm = URefTmp(refHarm+1);
            UmeasHarm = UMeasTmp(refHarm+1);
            OMeasHarm = OMeasTmp(refHarm+1);
            UrefExc = URefTmp(obj.refSig.excHarmonics+1);
            UmeasExc = UMeasTmp(obj.refSig.excHarmonics+1);
            OmeasExc = OMeasTmp(obj.refSig.excHarmonics+1);

            % Suppressed componenets
            suppHarmMul = obj.refSig.Hsupp; % Suppressed harmonic multiples
            suppHarmsTmp  = [];
            for hh = 1:length(suppHarmMul)
                suppHarmsTmp = [suppHarmsTmp ,suppHarmMul(hh):suppHarmMul(hh): obj.refSig.excHarmonics(end)];
            end
            suppHarms = sort(suppHarmsTmp); % Sort the suppressed harmonics in ascending order
            suppFreq_Hz = suppHarms/(obj.refSig.Nref*obj.refSig.Ts_s);
            UrefSupp = URefTmp(suppHarms+1);
            UmeasSupp = UMeasTmp(suppHarms+1);
            OmeasSupp = OMeasTmp(suppHarms+1);


            % Reference and measured temperature time and frequency plots
            figure()
            subplot(3,1,1)
            stairs(hours(obj.processedData.CaloricPeriod_degC.Time),[obj.refSig.refTempSig,obj.processedData.CaloricPeriod_degC.Caloric_Period],'. -');
            xlabel("Time [H]"); ylabel("Temperature [degC]")

            subplot(3,1,2)
            semilogx(obj.refSig.excFreq_Hz*1000,abs(UrefExc),'o r',refFreqVec*1000,abs(UrefHarm),'* b'); hold on;
            xlabel('Freq (mHz)'); ylabel({' Temperature FFT'; 'magnitude [-]'}); title(" Reference temperature spectrum")
            legend('Excited frequencies','All frequencies')
            if ~isempty(suppHarms)
                semilogx(suppFreq_Hz*1000,abs(UrefSupp),'o g'); hold on;
                legend('Excited frequencies','All frequencies', 'Suppressed frequencies')
            end
            subplot(3,1,3)
            semilogx(obj.refSig.excFreq_Hz*1000,abs(UmeasExc),'o r',refFreqVec*1000,abs(UmeasHarm),'* b'); hold on;
            xlabel('Freq (mHz)'); ylabel({' Temperature FFT'; 'magnitude [-]'}); title("Measured temperature spectrum")
            legend('Excited frequencies','All frequencies')
            if ~isempty(suppHarms)
                semilogx(suppFreq_Hz*1000,abs(UmeasSupp),'o g'); hold on;
                legend('Excited frequencies','All frequencies', 'Suppressed frequencies')
            end

            % OCV time and frequency plots
            figure
            subplot(2,1,1)
            stairs(hours(obj.processedData.ocvPeriod_V.Time),obj.processedData.ocvPeriod_V.OCV_Period,'. -');
            xlabel("Time [H]"); ylabel("OCV [V]")
            subplot(2,1,2)
            semilogx(obj.refSig.excFreq_Hz*1000,abs(OmeasExc),'o r',refFreqVec*1000,abs(OMeasHarm),'* b'); hold on;
            xlabel('Freq (mHz)'); ylabel({' Temperature FFT'; 'magnitude [-]'}); title("Measured OCV spectrum")
            legend('Excited frequencies','All frequencies')
            if ~isempty(suppHarms)
                semilogx(suppFreq_Hz*1000,abs(OmeasSupp),'o g'); hold on;
                legend('Excited frequencies','All frequencies', 'Suppressed frequencies')
            end
        end

        function PlotKernel(obj)

            freqExc = obj.refSig.excFreq_Hz*1000;

            figure
            db = @(x)20*log10(abs(x));
            phG = @(G) unwrap(angle(G))*180/pi;     % Function for transfer function phase
            ax(1) = subplot(3,1,1);
            semilogx(freqExc,db(obj.results.kernel),'. -')
            xlabel('Excited frequencies [mHz]'); ylabel('Magnitude [dB]'); legend('Kernel')
            ax(2) = subplot(3,1,2);
            semilogx(freqExc,phG(obj.results.kernel),'. -')
            xlabel('Excited frequencies [mHz]'); ylabel('Phase [deg]')
            ax(3) = subplot(3,1,3);
            semilogx(freqExc,db(obj.results.varKernel)/2,'. -')
            xlabel('Excited frequencies [mHz]'); ylabel('Kernel std [dB]'); legend('Kernel std')
            linkaxes(ax,'x')
        end

        function PlotKernelFit(obj)

            freqExc = obj.refSig.excFreq_Hz*1000;
            modelResp = obj.results.optimumFRF;

            figure
            db = @(x)20*log10(abs(x));
            phG = @(G) unwrap(angle(G))*180/pi;     % Function for transfer function phase

            kernel_phase = phG(obj.results.kernel);
            TF_phase = phG(modelResp);

            if any(kernel_phase < 0)
                kernel_phase = kernel_phase + 360;
            end
            if any(TF_phase < 0)
                TF_phase = TF_phase + 360;
            end
            
            titleStr = sprintf("Model Order: num %d denom %d.  Model fit: %.1f%%",obj.estimationSettings.modelOrder_num,obj.estimationSettings.modelOrder_denom,obj.results.fitMetrics.FitPercent);
            ax(1) = subplot(2,1,1);
            semilogx(freqExc,db(obj.results.kernel),'. -',freqExc,db(modelResp), '-')
            xlabel('Excited frequencies [mHz]'); ylabel('Magnitude [dB]'); legend('Kernel','TF fit'); title(titleStr)
            ax(2) = subplot(2,1,2);
            semilogx(freqExc,kernel_phase,'. -',freqExc,TF_phase,'-')
            xlabel('Excited frequencies [mHz]'); ylabel('Phase [deg]')
            linkaxes(ax,'x')

        end
    end

    methods (Access = private)

        %%%%%%%%%%%%%%%%%%%
        % Signal generation
        %%%%%%%%%%%%%%%%%%%
        function obj = MultiLevTempProfile(obj)

            % Design multilevel signal
            [Nref,harmVec,timeVec] = GenMultiLevN(obj);
            xSigTmp = FreqTimeOptimiser(Nref,harmVec,obj); % Dimensionless temperature
            quantLvls = linspace(-1,1,obj.refSig.numLevels);
            tempLvls = linspace(obj.refSig.Tmin_degC,obj.refSig.Tmax_degC,obj.refSig.numLevels);

            for ii = 1:length(xSigTmp)
                [~, idx] = min(abs(xSigTmp(ii) - quantLvls));
                tempSigTmp(ii,1) = tempLvls(idx);
            end
            %             tempSigTmp = multilev_new(Nref,harmVec,ones(size(harmVec)),obj.refSig.Hsupp,obj.refSig.numLevels,1,1000,[1:.1:4]);

            % Map to a time axis with Ts sampling
            Np = Nref*obj.refSig.timeStep_mins*60/obj.refSig.Ts_s;
            refTimeVec = [0:Np-1]'*obj.refSig.Ts_s;
            obj.refSig.refTimeVec_s = refTimeVec;
            obj.refSig.refTempSig = interp1(timeVec,tempSigTmp,refTimeVec,'previous','extrap');
            obj.refSig.Nref = Np;
            obj.refSig.excHarmonics = harmVec(:);
            obj.refSig.excFreq_Hz = harmVec(:)/(Np*obj.refSig.Ts_s);
        end

        function [NpClk,harmVec,timeVec,freqVec] = GenMultiLevN(obj)
            % Add doc string blurb

            Tp_s = obj.refSig.Tp_H*3600;                                % Convert the time period to seconds;
            clkTs_s = obj.refSig.timeStep_mins*60;                      % Convert the minimum step time to seconds;
            Ts_s = obj.refSig.Ts_s;
            suppHarmMul = obj.refSig.Hsupp;
            fMax = obj.refSig.fMax_Hz;


            NTmp = floor(Tp_s/clkTs_s);        % Number of clock pulses
            hMinClk = 1;

            if suppHarmMul == 0
                NpClk = NTmp;
            else
                prodHarmMul = prod(suppHarmMul);
                rem = mod(NTmp,prodHarmMul);
                NpClk = NTmp - rem;                     % With hamornc supression NpClk should be a multiple of the suppressed harmonics
            end

            f0 = 1/(NpClk*clkTs_s);                     % Fundemental frequency
            hMaxClk = floor(fMax/f0);


            if suppHarmMul == 0
                harmVec = [hMinClk:hMaxClk]';
            else
                harmVecTmp = [hMinClk:hMaxClk]';
                for ii = 1:length(suppHarmMul)
                    hM = suppHarmMul(ii);
                    harmVecTmp = setdiff(harmVecTmp,[hM:hM:hMaxClk]);
                    harmVec = harmVecTmp';
                end
            end

            timeVec = [0:NpClk-1]'*clkTs_s;
            freqVec = harmVec/(NpClk*clkTs_s);

            if NpClk < 2*hMaxClk
                error('Maximum specified frequency must be less than %2.2E Hz',NpClk*f0/2); % Np/(2*Tp*3600)
            end
            if mod(clkTs_s,Ts_s) ~= 0
                tStep = clkTs_s - mod(clkTs_s,Ts_s);
                error('Sampling time (Ts) is not a multiple of timeStep.\nNearest multiple for timeStep_mins: %3ds',tStep);
            end
            if length(harmVec) < 6
                Tp_lowerLimit = 15/freqVec(end); % (fmax - f0)/f0 >= 14. Need at least seven excited harmonics to estimate kernel function with transients and times two becuase even harmonics are suppressed
                Tp_lowerLimit_H = ceil(Tp_lowerLimit/clkTs_s)*clkTs_s/3600; % Round up to closest minimum signal time step.
                error('Number of excited harmonics is less than 6\nEither increase Tp_H to %.2f hours',Tp_lowerLimit_H);
            end

        end

        function xQV = FreqTimeOptimiser(Nref,harmVec,obj)
            % Swtich between time and frequency to generate the desired
            % signal

            XVec0 = zeros(Nref,1);                  % Initialise a zero vector to hold the frequency components
            amp = ones(length(harmVec),1);          % Unit amplitude is desired
            phi = -harmVec'.*(harmVec' - 1)*pi/length(harmVec);  % Initialise with Schroeder phases

            XVec0(harmVec+1) = amp.*exp(1i*phi);
            xVec0 = 2*real(ifft(XVec0));
            xVec0 = xVec0/max(abs(xVec0));          % Scale to unity

            xVec = xVec0; % Initialise

            for ii = 1:1
                xQV = QuantizeSignal(xVec,obj);
                XVec = fft(xQV);                    % FFT of quntised signal to subsitute flat spectrum
                phi = angle(XVec);                  % Get phases of discretised signal
                XVec(harmVec+1) = amp.*exp(1i*phi(harmVec+1));
                xVec = 2*real(ifft(XVec));
                xVec = xVec/max(abs(xVec));          % Scale to unity
            end
        end

        function xQV = QuantizeSignal(xVec,obj)
            % Quantise a signal to the desired levels and enforce
            % contraints
            quantLvls = linspace(-1,1,obj.refSig.numLevels);
            tempLvls = linspace(obj.refSig.Tmin_degC,obj.refSig.Tmax_degC,obj.refSig.numLevels);
            [~, idxTstart] = min(abs(obj.refSig.startTemp - tempLvls));

            if tempLvls(idxTstart) ~= obj.refSig.startTemp
                warning("Based on the specified Tmin: %.1fdegC, Tmax: %.1fdegC and number of signal levels: %d, actual start temperature will be: %.1fdegC, requested start temperature is: %.1fdegC",obj.refSig.Tmin_degC,obj.refSig.Tmax_degC,obj.refSig.numLevels,tempLvls(idxTstart), obj.refSig.startTemp )
            end

            % Generate quantised signal
            for ii = 1:length(xVec)
                [~, idx] = min(abs(xVec(ii) - quantLvls));
                xQVTmp(ii,1) = quantLvls(idx);
            end

            % Translate to requested starting (or closest) temperature
            idxShift = find(xQVTmp == quantLvls(idxTstart),1);
            xQV = [xQVTmp(idxShift:length(xVec));xQVTmp(1:idxShift-1)];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update reference signal propeties based on any preivous ref
        % signals (legacy data)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = UpdateRefSig(obj,legacyRefSig)
            arguments
                obj EntropyCoeffEstimator
                legacyRefSig struct
            end
            obj.refSig.fMax_Hz = legacyRefSig.u.fMax;
            obj.refSig.Hsupp = legacyRefSig.u.Hsupp;
            obj.refSig.Nref = legacyRefSig.u.Nref;
            obj.refSig.numLevels = legacyRefSig.u.Nl;
            obj.refSig.refTempSig = legacyRefSig.u.refTempSig;
            obj.refSig.refTimeVec_s = legacyRefSig.u.refTimeVec;
            obj.refSig.startTemp = legacyRefSig.u.refTempSig(1);
            obj.refSig.timeStep_mins = legacyRefSig.u.timeStep_s/60;
            obj.refSig.Tmax_degC = legacyRefSig.u.Tmax;
            obj.refSig.Tmin_degC = legacyRefSig.u.Tmin;
            obj.refSig.Tp_H = legacyRefSig.u.Tp;
            obj.refSig.Ts_s = legacyRefSig.u.Ts;
            obj.refSig.excHarmonics = legacyRefSig.u.harmVec(:);
            obj.refSig.excFreq_Hz = legacyRefSig.u.harmVec(:)/(legacyRefSig.u.Nref * legacyRefSig.u.Ts);
        end

        function obj = ReshapeAveragePeriods(obj)
            numSigSamples = length(obj.measData.ocvData.Time_s);
            numPeriods = floor(numSigSamples/obj.refSig.Nref);
            numOffSamples = rem(numSigSamples,obj.refSig.Nref);
            numOfTempChannels = size(obj.measData.tempData,2);

            if numOffSamples > 0
                if numPeriods > 1 singPlrStr = "periods"; else singPlrStr = "period"; end
                wrnMsg = ["Measured data are not an integer multiple of the reference signal period\nMeasured data has %d %s and %d extra samples"];
                warning(wrnMsg,numPeriods,singPlrStr,numOffSamples)
            end

            totalNumSamples =  numPeriods*obj.refSig.Nref;
            idxRng = numSigSamples - totalNumSamples +1 : numSigSamples;

            % Get integer number of periods
            timeVecTmp = obj.measData.tempData.Time_s(idxRng);

            % If more that one period is available then discard first
            % period and average remaining periods
            try obj.estimationSettings.usePeriods;
                if isempty(obj.estimationSettings.usePeriods)
                    if numPeriods > 1
                        pAveIdx = 2:numPeriods;
                    else
                        pAveIdx = 1;
                    end
                else
                    if obj.estimationSettings.usePeriods > 1
                        pAveIdx = 2:obj.estimationSettings.usePeriods;
                    else
                        pAveIdx = 1;
                    end

                end
            catch
                if numPeriods > 1
                    pAveIdx = 2:numPeriods;
                else
                    pAveIdx = 1;
                end
            end

            % Get integer number of periods and reshape to Nref x P
            ocvVecTmp = obj.measData.ocvData.OCV_V(idxRng);
            ocvMatrix = reshape(ocvVecTmp,obj.refSig.Nref,numPeriods);

            % Average OCV over periods
            obj.processedData.ocvPeriod_V = timetable(timeVecTmp(1:obj.refSig.Nref),mean(ocvMatrix(:,pAveIdx),2),'VariableNames',{'OCV_Period'});


            % Get integer number of periods and reshape to Nref x P
            caloricTmp = obj.processedData.Caloric.Caloric(idxRng);
            caloricMatrix = reshape(caloricTmp,obj.refSig.Nref,numPeriods);

            obj.processedData.numPeriods = numPeriods;

            % Average Caloric temperature over periods
            obj.processedData.CaloricPeriod_degC = timetable(timeVecTmp(1:obj.refSig.Nref),mean(caloricMatrix(:,pAveIdx),2),'VariableNames',{'Caloric_Period'});

            % Get integer number of periods and reshape to Nref x P and
            % average measured temperatures
            tempTimeTable = timetable();
            for tt = 1 : numOfTempChannels
                tempVecTmp = obj.measData.tempData.(tt)(idxRng);
                tempChannelName = [obj.measData.tempData.Properties.VariableNames(tt) + "_Period"];
                tempMatrix = reshape(tempVecTmp,obj.refSig.Nref,numPeriods);
                tempVec = mean(tempMatrix(:,pAveIdx),2);
                tempTimeTable = [tempTimeTable, timetable(timeVecTmp(1:obj.refSig.Nref),tempVec,'VariableNames',tempChannelName)];
            end

            obj.processedData.tempData = tempTimeTable;
        end

        function obj = CaloricCellTemperature(obj)

            alphaCuCell = obj.cellThermalProperties.copperConductivity_wpmk/obj.cellThermalProperties.copperThickness_m; %Approximation heat conductivity Cu-plate to Cell, W/(m^2*K)
            KeyFigures.Nu2Infinity = (pi^2)/2;
            KeyFigures.Bi = (alphaCuCell * obj.cellThermalProperties.cellThickness_m)/ obj.cellThermalProperties.cellConductivity_wpmk; %Biot-Number
            KeyFigures.NuiInfinity = (4+obj.cellThermalProperties.GeoFactor + KeyFigures.Bi)/(1+(KeyFigures.Bi/KeyFigures.Nu2Infinity));

            % % IMPROVE: Write a generic function to select which channels to average
            try
                MeanTemp.TECRef = (obj.measData.tempData.TECRefOben + obj.measData.tempData.TECRefUnten)/2;
                MeanTemp.Caloric = (obj.measData.tempData.LIBMitteOben + obj.measData.tempData.LIBPlusOben + obj.measData.tempData.LIBMinusOben...
                                    + obj.measData.tempData.LIBMitteUnten + obj.measData.tempData.LIBMinusUnten)/5; %für erste Zeile kann wegen t=0 keine kalorische Mitteltemperatur berechnet werden
                setTemp = obj.measData.tempData.SollTemp;
            catch
                MeanTemp.TECRef = (obj.measData.tempData.TEC1 + obj.measData.tempData.TEC2)/2;
                MeanTemp.Caloric = (obj.measData.tempData.SurfaceTopCenter + obj.measData.tempData.SurfaceTopCathode + obj.measData.tempData.SurfaceTopAnode...
                                    + obj.measData.tempData.SurfaceBottomCenter + obj.measData.tempData.SurfaceBottomAnode)/5; %für erste Zeile kann wegen t=0 keine kalorische Mitteltemperatur berechnet werden
                numPeriods = floor(length(MeanTemp.TECRef)/obj.refSig.Nref) + 1;
                tmpSetTemp = repmat(obj.refSig.refTempSig,numPeriods,1);
                setTemp = tmpSetTemp(1:length(MeanTemp.TECRef));
            end
            Zeilen = size(MeanTemp.TECRef,1);
            t = 2;
            Caloric(1,1) = 0;
            for ii = 2:Zeilen
                if isequal(setTemp(ii), setTemp(ii-1))
                    %Kalorische Mitteltemperatur
                    KeyFigures.Fo(ii) = (obj.cellThermalProperties.cellConductivity_wpmk*t)/(obj.cellThermalProperties.cellThickness_m^2 * obj.cellThermalProperties.cellDenisty_kgpm3 * obj.cellThermalProperties.cellSpecificHeatCapacity_Jp);
                    KeyFigures.Nui0(ii) = (sqrt(pi)+ 10*KeyFigures.Bi*sqrt(KeyFigures.Fo(ii)))/(1+5*KeyFigures.Bi*sqrt(pi*KeyFigures.Fo(ii)))*1/sqrt(KeyFigures.Fo(ii));
                    KeyFigures.Nui(ii) = sqrt((KeyFigures.NuiInfinity^2)-(0.4^2)+(KeyFigures.Nui0(ii) + 0.4)^2);
                    KeyFigures.NTUi(ii) = (obj.cellThermalProperties.GeoFactor * KeyFigures.Fo(ii))/((1/KeyFigures.Bi) + (1/KeyFigures.Nui(ii)));
                    Caloric(ii,1) = (MeanTemp.TECRef(ii) + (MeanTemp.Caloric(ii-(t/2))-MeanTemp.TECRef(ii))*exp(-KeyFigures.NTUi(ii)));
                    t = t+2;
                else
                    t = 2;
                    %Kalorische Mitteltemperatur
                    KeyFigures.Fo(ii) = (obj.cellThermalProperties.cellConductivity_wpmk*t)/(obj.cellThermalProperties.cellThickness_m^2 * obj.cellThermalProperties.cellDenisty_kgpm3 * obj.cellThermalProperties.cellSpecificHeatCapacity_Jp);
                    KeyFigures.Nui0(ii) = ((sqrt(pi)+ 10*KeyFigures.Bi*sqrt(KeyFigures.Fo(ii)))/(1+5*KeyFigures.Bi*sqrt(pi*KeyFigures.Fo(ii))))*(1./sqrt(KeyFigures.Fo(ii)));
                    KeyFigures.Nui(ii) = sqrt(KeyFigures.NuiInfinity^2-0.4^2.+(KeyFigures.Nui0(ii) + 0.4)^2);
                    KeyFigures.NTUi(ii) = (obj.cellThermalProperties.GeoFactor * KeyFigures.Fo(ii))/((1/KeyFigures.Bi) + (1/KeyFigures.Nui(ii)));
                    Caloric(ii,1) =(MeanTemp.TECRef(ii) + (MeanTemp.Caloric(ii-1) - MeanTemp.TECRef(ii))*exp(-KeyFigures.NTUi(ii)));
                end

            end
            obj.processedData.Caloric = timetable(obj.measData.tempData.Time_s,Caloric,'VariableNames',"Caloric");
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Functions for FRF and transfer function estimation. Gain of
        % transfer function is dUdT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = estimateFRF(obj)
            % Make sure that the data is already loaded
            if isempty(fieldnames(obj.processedData))
                error(sprintf("No temperature or OCV data is imported to perform analysis. Import data via:\n ImportOCVData() and \n ImportTemperatureData()"))
            end

            UMeasTmp = fft(obj.processedData.CaloricPeriod_degC.Caloric_Period);
            OMeasTmp = fft(obj.processedData.ocvPeriod_V.OCV_Period);

            idxExc = obj.refSig.excHarmonics + 1;
            UmeasExc = UMeasTmp(idxExc);
            OmeasExc = OMeasTmp(idxExc);


            if strcmp(obj.estimationSettings.transientOnOff," ") && obj.processedData.numPeriods > 1
                if isempty(obj.estimationSettings.usePeriods)
                    tranOnOff = 0;
                    obj.estimationSettings.transientOnOff = "off";
                elseif obj.estimationSettings.usePeriods(1) > 1
                    tranOnOff = 0;
                    obj.estimationSettings.transientOnOff = "off";
                else
                    tranOnOff = 1;
                    obj.estimationSettings.transientOnOff = "on";
                end
            else
                tranOnOff = 1;
                obj.estimationSettings.transientOnOff = "on";
            end


            [G,~,~,Cg] = LPMethod(obj,UmeasExc,OmeasExc,idxExc,'transient',tranOnOff);

            obj.results.kernel = G;        % Estimated temperature to OCV frequency response (kernel)
            obj.results.varKernel = Cg;    % Estimated temperature to OCV frequency response variance
        end

        function fitFRF(obj)

            kernel = obj.results.kernel;
            Cg = obj.results.varKernel;
            angFreq_rads = 2*pi*obj.refSig.excFreq_Hz;
            Ts = obj.refSig.Ts_s;


            idxE = 1:2:length(angFreq_rads);
            idxV = 2:2:length(angFreq_rads);
            fdObj = idfrd(kernel,angFreq_rads,Ts); % Create a frequency response obj
            fdObjE = idfrd(kernel(idxE),angFreq_rads(idxE),Ts); % Create a frequency response obj
            fdObjV = idfrd(kernel(idxV),angFreq_rads(idxV),Ts); % Create a frequency response obj

            if isempty(obj.estimationSettings.modelOrder_num)
                nb = [1:10];    % Range for numerator order
            else
                nb = obj.estimationSettings.modelOrder_num;
            end
            if isempty(obj.estimationSettings.modelOrder_denom)
                na = [1:10];    % Range for denomenator order
            else
                na = obj.estimationSettings.modelOrder_denom;
            end

            optOrder = selstruc(arxstruc(fdObjE,fdObjV,struc(na,nb,0)),0); % Estimate optimum order

            % Update with optimium model orders
            obj.estimationSettings.modelOrder_num = optOrder(2);
            obj.estimationSettings.modelOrder_denom = optOrder(1);
            %
            %
            %                 % TF estimation settings
            %                 tfEstOpt = tfestOptions('EnforceStability',1);
            %                 tfEstOpt.InitializeOptions.MaxIterations = 7000;
            %
            %
            %                 % Estimate discrete time transfer function
            %                 optTF = tfest(fdObj,optOrder(1),optOrder(2),'Ts',Ts,tfEstOpt);
            %
            %                 numMean = sum(optTF.Numerator);
            %                 denomMean = sum(optTF.Denominator);
            %                 coeffSTD = sqrt(diag(optTF.Report.Parameters.FreeParCovariance));         % Standard deviations of numerator and denominator coefficients
            %                 numSTD = sum(coeffSTD(1:optOrder(2)));
            %                 denomSTD = sum(coeffSTD(optOrder(2)+1:end));
            %                 numFracUncertanity = abs(numSTD/numMean);
            %                 denomFracUncertanity = abs(denomSTD/denomMean);

            %                 obj.results.optimumTF = optTF;
            %                 obj.results.fitMetrics = optTF.Report.Fit;
            %                 obj.results.dUdT_mVpK = numMean/denomMean*1000;                                                   % dUdT term [mV/K]. Steady-state gain of discrete time transfer function
            %                 obj.results.dUdT_std = abs((numMean/denomMean))*(numFracUncertanity + denomFracUncertanity)*1000;  % Std of dUdT

            %%% Go back to previous estimation code %%%
            % Parameterise the frequency response.
            % Initial TF estimation using Levi method

             wNorm = angFreq_rads*Ts;
            leviTF = LeviAlgorithm(obj,kernel,wNorm,optOrder(2),optOrder(1),"fs",1/Ts);       % Transfer function fi using Levi method
            theta0 = [leviTF.B;leviTF.A];

            fcnTF = @(theta,angFreq_rads)freqs(theta(1:optOrder(2)+1),[theta(optOrder(2)+2:end);1],angFreq_rads);
            fcnGReIm = @(theta,angFreq_rads) [real(fcnTF(theta,angFreq_rads));imag(fcnTF(theta,angFreq_rads))];

            [optimumTf,tfResults] = LMAlgorithm(obj,fcnGReIm, theta0, angFreq_rads, [real(kernel);imag(kernel)],"s",sqrt([Cg; Cg]/2));

            GOpt = freqs(optimumTf(1:optOrder(2)+1),[optimumTf(optOrder(2)+2:end);1],angFreq_rads);
            SSE = sum(abs(GOpt-kernel).^2);
            SST = sum(abs(kernel-mean(kernel)).^2);
            tfResults.FitPercent = (1- SSE/SST)*100;

            obj.results.optimumTF = optimumTf;
            obj.results.optimumFRF = GOpt;
            obj.results.fitMetrics = tfResults;
            obj.results.dUdT_mVpK = optimumTf(optOrder(2)+1)*1000;                                                   % dUdT term [mV/K]. Steady-state gain of discrete time transfer function
            obj.results.dUdT_std = tfResults.stdTheta(optOrder(2)+1)*1000;  % Std of dUdT



            %                 if log10(obj.results.dUdT_std) >= 1 && obj.results.cntr < 1
            %                     obj.results.cntr = obj.results.cntr + 1;
            %                     obj.EstimateEntropyCoeff("modelOrder_denom",1,"modelOrder_num",1,"transientOnOff",obj.estimationSettings.transientOnOff);
            %                 end

        end

        function [G, T, Cv ,Cg, Ct, alpha] = LPMethod(obj,X,Y,F,nvp)
            %
            % Local polynomial method function to estimate a frequency response
            % function
            % Madotory inputs:
            %   X - FFT of input signal size lf x 1
            %   Y - FFT of output signal size lf x 1
            %   F - FFT lines size lf x 1
            %
            % Optional arguments:
            %   nvp.order - order of local polynomial default value 2, size 1 x 1
            %   nvp.transient - Set to 1 if transients are to be included in the output sprectrum (default method.transient = 1)
            %
            % Outputs:
            %   G - estimated FRF, size lf x 1
            %   T - Estimated Transients, size lf x 1
            %   Cv - Noise variance, size lf x 1
            %   Cg - FRF variance, size lf x 1
            %   Ct - Transient variance, size lf x 1
            %   alpha - polynomial variables size lf x 2*poly_order
            %
            % Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 28/06/11-20/10/2021 (~)
            % All Rights Reserved
            % Software may be used freely for non-comercial purposes only

            arguments
                obj
                X
                Y
                F
                nvp.transient = 1;
                nvp.polyOrder = 2;
            end

            X = X(:);
            Y = Y(:);
            F = F(:);


            % Initialise
            lf = length(F);
            G = zeros(lf,1);
            T = zeros(lf,1);
            Cv = zeros(lf,1);
            Cg = zeros(lf,1);
            Ct = zeros(lf,1);

            poly_order = nvp.polyOrder;
            transient = nvp.transient;

            if transient == 1
                alpha=zeros(lf,2*poly_order);
            else
                alpha=zeros(lf,poly_order);
            end

            if transient == 1
                ntheta= 2*poly_order+2;
            else
                ntheta= poly_order+1;
            end

            % Use 2*R+1 frequencies to estimate the ntheta parameters
            if mod(ntheta,2)==0
                R = ntheta/2;
            else
                R = (ntheta + 1)/2;
            end

            Rb = R; % use Rb frequencies before kk
            Ra = R; % use Ra frequencies after kk

            for kk = 1:length(F)
                %Slice input and output over frequency range
                if (kk<=Rb)
                    Fslice = F(1:2*R+1);
                    Xslice = X(1:2*R+1);
                    Yslice = Y(1:2*R+1);
                end
                if (kk>Rb) && (kk<=lf-Ra)
                    Fslice = F(kk-Rb:kk+Ra);
                    Xslice = X(kk-Rb:kk+Ra);
                    Yslice = Y(kk-Rb:kk+Ra);
                end
                if (kk>lf-Ra)
                    Fslice = F(lf-2*R:lf);
                    Xslice = X(lf-2*R:lf);
                    Yslice = Y(lf-2*R:lf);
                end

                % Regressor matrix
                Kn = zeros(2*R+1,ntheta);

                for rr = 1:2*R+1
                    %Taylor series variable r
                    r_Taylor = Fslice(rr)-F(kk);
                    for ss = 1:poly_order             % Taylor series model order
                        r_power(1,ss) = r_Taylor^ss;  % r_power for the transient taylor expansion
                        Xr = r_power*Xslice(rr);      % Xr for the frf taylor expansion
                    end
                    if transient == 1 % estimate with transients
                        Kn(rr,:)=[Xslice(rr), Xr, 1, r_power];
                    else
                        Kn(rr,:)=[Xslice(rr), Xr];
                    end
                end

                %Call numerically stable linear least squares function
                [Theta,lls_results] = LLs(obj,Kn,Yslice);

                %Extract FRF's and their variances
                G(kk,1) = Theta(1);                     % LTI branch dynamics
                Cg(kk,1) = lls_results.paraVar(1);          % Variance of G
                Cv(kk,1) = lls_results.noiseVar;            % Noise variance

                if transient == 1
                    T(kk,1) = Theta(poly_order+2);      % Transient spectrum
                    Theta([1,poly_order+2])=[];
                    alpha(kk,:) = Theta;                % polynomial variables
                    Ct(kk,1) = lls_results.paraVar(poly_order+2);    % Transient variance
                else
                    T(kk,1) = nan;                      % Transient spectrum
                    Theta(1) = [];
                    alpha(kk,:) = Theta;                % polynomial variables
                    Ct(kk,1) = nan;                     % Transient variance
                end
            end

        end

        function Model = LeviAlgorithm(obj,G,w,nb,na,nvp)

            arguments
                obj,
                G,
                w,
                nb,
                na,
                nvp.domain = "s",
                nvp.fs = 1,
            end

            fs = nvp.fs;

            if nvp.domain == 's'
                v = 1i*w*fs;
            else
                v = exp(-1i*w);
            end


            F = length(w);

            %Construct regressor matrix
            for kk = 1:F
                Reg(kk,:) = [-G(kk)*(v(kk).^(na:-1:1)) v(kk).^(nb:-1:0)];
            end

            Regt = [real(Reg);imag(Reg)];
            Gt = [real(G);imag(G)];

            [theta, results] = LLs(obj, Regt,Gt);

            Model.A = theta(1:na);
            Model.B = theta(na+1:end);
            Model.results = results;
        end

        function [theta,results] = LMAlgorithm(obj, fh, theta0, u, d, nvp)
            %
            % Perform nonlinear least squares optimisation using the
            % Levenberg-Marquardt method
            %
            % Minimised cost-function f(theta)
            %   f(theta) = sum ((G(theta)_i-d_i)/s_i)^2  i = 1...M
            %
            % Mandotory input argumetns
            %   fh: function handle of nonlinear model. fh = @(theta, u)fcn(theta,u,optArg1,...,optArgN).
            %       Nonlinear model function should return two outputs, the model output and Jacobian as second
            %   theta0: Initial starting point of model parameters, size N x 1
            %   u: Input signal to simulate nonlinear model, size Mu x 1
            %   d: Measured output data, size M x 1
            %
            % Optional arguments. Create a structure variable with the following fields:
            %   Jacobian: Specify as 'on' if model function returns the Jacobian or to 'off' to approximate by finte forward difference, default Jacobian ='off'.
            %   s: Residual weights, normally std of noise, default s = 1, size M x 1
            %   iterMax: Maximum number of iterations, default iterMax = 1000, double size 1 x 1
            %   TolDT:  Termination tolerance of parameter update, default TolDT = 1E-6, double size 1 x 1
            %   diagnosis: Set daignosis to 'on' to plot cost-function and lambda vs iterations, default diagnosis = 'off'
            %   epsilon: if Jacobian is 'off', use epsilon for parameter incerement to calculate approximate Jacobian, default epsilon = 1E-6, double size 1 x 1
            %   dispMsg: Specifiy as 'on' or 'off'. If 'on' messages will be printed
            %            else no messages are printed. Default 'on'.
            %
            % Output arguments:
            %   theta: Optimised parameter vector, size N x 1
            %   results: Sturcture variable with fields
            %            - covTheta: Covariance matrix of optimum parameters, size N x N
            %            - stdTheta: Standard deviation of optimum parameters, size N x 1
            %            - fracEr: Fractional error of optimum parameters, size N x 1
            %            - cF_iter: cost-function value at each iteration
            %            - L_iter: Lambda weight at each iteration
            %            - termMsg: Reason for Levenberg-Marquardt iteration termination
            %            - rankMsg: Message stating rank of Levenberg-Marquardt regressor for each iteration, size ~ x 1
            %            - LMRankFull: A flag at each iteration with value 0 or 1 if regressor is rank deficient or not, size ~ x 1
            %
            % Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 14/10/2015 (Highway to hell!!)
            % All Rights Reserved
            % Software may be used freely for non-comercial purposes only

            arguments
                obj,
                fh,
                theta0,
                u,
                d,
                nvp.Jacobian = "off",
                nvp.s = ones(size(d)),
                nvp.iterMax = 1000,
                nvp.TolDT = 1E-6,
                nvp.diagnosis = "off",
                nvp.epsilon = 1E-6,
                nvp.dispMsg = "on";
            end

            theta0 = theta0(:);
            nData = length(d);
            nPara = length(theta0);


            % Initialise
            theta_prev = theta0;
            JacStatus = lower(nvp.Jacobian);
            Jw = spdiags(1./nvp.s,0,nData,nData); % Residual weights (noise std) of cost-function to scale each row of the Jacobain


            % Evalaute model function and Jacobian for initial parameter values
            if ismember(JacStatus,'on')
                [y,J_prevTmp] = fh(theta_prev,u);
            elseif ismember(JacStatus ,'off')
                y = fh(theta_prev,u);
                J_prevTmp = JacApprox(obj,theta_prev,y,u,fh,nvp.epsilon);
            end
            J_prev = Jw*J_prevTmp; % Scale Jacobian with residual weights

            [cF_prev,F_prev] = costFunctionEval(obj,y,d,nvp.s);

            cF = cF_prev*10; % Induce that the present cost-function is worse than previous one with the assumed value of lambda
            iterUpdate = 1;
            deltaT = 1E6;
            lambda = 10;
            innerLoop = 1; % Used for de-bugging purposes

            cF_iter(iterUpdate,1) = cF;
            L_iter(iterUpdate,1) = lambda;

            % Start solution update
            while norm(deltaT) > nvp.TolDT && iterUpdate <= nvp.iterMax
                if isinf(y)
                    error('\nLM Model output is infinity at iteration %d\n',iterUpdate);
                end
                while cF > cF_prev % Increase lambda and re-evaluate parameter update
                    lambda = lambda*10;
                    deltaT = parameterUpdate(obj,J_prev,F_prev,lambda);     % Calculate parameter update
                    theta = theta_prev + deltaT;                        % Update parameter estimate

                    % Evalaute model function and Jacobian for updated parameter
                    if ismember(JacStatus ,'on')
                        [y, J_Tmp] = fh(theta,u);
                    elseif ismember(JacStatus ,'off')
                        y = fh(theta,u);
                        J_Tmp = JacApprox(obj,theta,y,u,fh,nvp.epsilon);
                    end
                    J = Jw*J_Tmp;                                        % Scale Jacobian with residual weights
                    [cF, F] = costFunctionEval(obj,y,d,nvp.s);
                    innerLoop = innerLoop + 1;
                end

                cF_prev = cF;
                theta_prev = theta;
                J_prev = J;
                F_prev = F;

                lambda = lambda/10;
                [deltaT,regRank] = parameterUpdate(obj,J,F,lambda);       % Calucuate parameter update
                theta = theta + deltaT;                               % Update parameter estimate

                % Evalaute model function and Jacobian for updated parameter
                if ismember(JacStatus ,'on')
                    [y,J_Tmp] = fh(theta,u);
                elseif ismember(JacStatus ,'off')
                    y = fh(theta,u);
                    J_Tmp = JacApprox(obj,theta,y,u,fh,nvp.epsilon);
                end
                J = Jw*J_Tmp;                                        % Scale Jacobian with residual weights
                [cF, F] = costFunctionEval(obj,y,d,nvp.s);

                cF_iter(iterUpdate,1) = cF;
                L_iter(iterUpdate,1) = lambda;

                if regRank.unique == 0
                    results.rankMsg{iterUpdate,1} =  sprintf(['LM Regressor rank deficient at iteration %d %s'],iterUpdate,regRank.msg);
                    if strcmp(nvp.dispMsg,'on')
                        fprintf(['\nLM Regressor rank deficient at iteration %d %s'],iterUpdate,regRank.msg);
                    end
                    results.LMRankFull(iterUpdate,1) = regRank.unique;
                else
                    results.rankMsg{iterUpdate,1} = sprintf(['LM Regressor preserves full rank at iteration %d %s'],iterUpdate,regRank.msg);
                    results.LMRankFull(iterUpdate,1) = regRank.unique;
                end

                iterUpdate = iterUpdate + 1;
            end % End of main while iterative loop
            iterUpdate = iterUpdate - 1; % Reduce iteration count by one when loop is exited

            % Estiamte parameter covariance matrix
            if all(nvp.s)                    % If measurement variance is not used in the cost function eistamte from residue for paramter variance scaling
                sCF = cF/(nData-nPara);
            else
                sCF = 1;                     % Else if measurement variance is used in the cost function, paramter variance estimate does not need scaling
            end
            covTheta = CovTheta(obj,sCF,J);      % Parameter variance

            if ismember(lower(nvp.diagnosis),'on')
                Diagnosis(obj,cF_iter,L_iter,iterUpdate)
            end

            idx = [norm(deltaT) < nvp.TolDT, iterUpdate == nvp.iterMax];
            termStr = {[' Parameter update is smaller than specified tolerance, TolDT = ', num2str(nvp.TolDT),'.'],...
                [' Maximum iteration reached, iterMax = ', num2str(nvp.iterMax),'.']};
            if strcmp(nvp.dispMsg,'on')
                fprintf('\n\nIteration terminated: %s\n',termStr{idx});
            end

            results.covTheta = covTheta;
            results.stdTheta = sqrt(diag(covTheta));
            results.cF_iter = cF_iter;
            results.L_iter = L_iter;
            results.Jacobian = J;
            results.termMsg = ['Iteration terminated: ',termStr{idx}];
            results.fracErr = results.stdTheta./abs(theta); % Fractional error
        end
        
        function [cF,F] = costFunctionEval(obj,y,d,s)
            arguments
                obj,y,d,s
            end
            F = (y-d)./s;       % Weighted residual
            cF = norm(F)^2;     % Cost-function
        end

        function [deltaT,regRank] = parameterUpdate(obj,J,F,lambda)
            arguments
                obj,J,F,lambda
            end
            K = ((J')*J + lambda*diag(diag((J')*J)));       % Create LM regressor matrix (cost-function Hessian + Steepest descent)
            Z = (-J')*F;                                    % Negative cost-function gradient
            [deltaT,regRank] = LLs(obj,K,Z);                    % Call numerically stable linear least squares method
        end

        function Diagnosis(obj,cF,L,I)
            arguments
                obj,cF,L,I
            end
            figure()
            semilogx([0:I-1],cF,'.-')
            xlabel('Iteration number')
            ylabel('Cost-fucntion (y-d/s)^2')

            figure()
            plot([0:I-1],L,'.-')
            xlabel('Iteration number')
            ylabel('Steepest descent lambda factor')
        end

        function [theta,results] = LLs(obj,K,Z,nvp)
            %
            % Computes a numerically stable (weighted) linear leaast squares estimate
            % and returns the optimal estimate, its variance and an estimate of the
            % noise variance. Parameters will be non-unique if K is rank deficient.
            %
            % Inputs (mandatory):
            %   K: Regresor matrix, size n x m
            %   Z: Output vector, size n x 1
            %
            % Inputs (optional):
            %   ny: Output noise vector, std of output noise for weighted least
            %   squares, size n x 1. default ny = ones(n,1);
            %   plotFit: Set plot as 1 or zero to look at the fitted vs data
            %
            % Outputs:
            %   theta: Optimum parameter estimate, size m x 1
            %   results: Structure variable with following fields
            %            - paraVar: Variance estimate of theta, size m x 1
            %            - paraCov: Parameter covariance matrix size m x m
            %            - noiseVar: Estimate of noise variance, size 1 x 1
            %            - resVec: Residual vector, size n x 1
            %            - resNorm: 2 norm of residual vector, size 1 x 1
            %            - regMsg: Message stating rank of regressor
            %            - regRankFull: A flag with value 0 or 1 if regressor is rank deficient or not.
            %
            % Copyright (C) W. D. Widanage -  WMG, University of Warwick, U.K. 11-02-2012-12/01/2016 (Through the never)
            % All Rights Reserved
            % Software may be used freely for non-comercial purposes only

            arguments
                obj
                K
                Z
                nvp.ny = ones(size(Z))
                nvp.plotFit = 0
            end

            ny = nvp.ny;
            plotFit = nvp.plotFit;

            % Weight matrix
            W = spdiags(1./ny,0,length(ny),length(ny));

            %Multiply by weights
            Z = W*Z;
            K = W*K;

            %Normalise K with the 1/norm of each column
            Knorm = sqrt(sum(abs(K).^2));
            idxZeros = Knorm<1E-14;
            Knorm(idxZeros) = 1;
            N = diag(1./Knorm);
            Kn = K*N;

            %Compute Lls via economic SVD decompostion
            [U, S, V] = svd(Kn,0);    % Perform SVD
            ss = diag(S);             % Singular values
            idxZeros = ss < 1E-14;    % Replace any small singular value with inf
            nCol = size(Kn,2);
            if sum(idxZeros)>0        % If there are zero singular values sum(idxZeros) > 0 and regressor is rank deficient
                results.msg = sprintf('Estimated parameters are NON-UNIQUE.\n Lls regressor is rank defficient. Rank = %d instead of %d. \nParameters estimated from a subspace.',nCol-sum(idxZeros),nCol);
                fprintf('\nEstimated parameters are NON-UNIQUE.\nRegressor rank defficient. Rank = %d instead of %d. \nParameters estimated from a subspace.\n\n',nCol-sum(idxZeros),nCol)
                results.regRankFull = 0;
                results.unique = 0;
            else
                results.msg = sprintf('Estimated parameters are unique.\nRegressor preserves full rank. Rank =  %d.',nCol);
                results.regRankFull = 1;
                results.unique = 1;
            end
            ss(idxZeros) = inf;
            Sv = diag(1./ss);               % Inverse singular value matrix

            % Least squares solution
            theta = N*V*Sv*U'*Z;
            % cond(N*V*Sv*U')

            %Projection matrix and residuals
            P = (eye(length(Z))-U*(U')); % Projection matrix
            R = P*Z;                     % Residuals

            %Estimate of noise and parameter variance
            if all(ny)                                                  % If measurement variance is not provided estiamte via residuals
                results.noiseVar = ((R')*R)/real(trace(P));             % Noise variance
                results.paraCov = (N')*V*Sv*Sv*(V')*N*results.noiseVar; % Parameter covariance
            else                                                        % Else no need to estimate noise variance
                results.noiseVar = [];                                  % Noise variance
                results.paraCov = (N')*V*Sv*Sv*(V')*N;                  % Parameter covariance
            end

            if plotFit
                figure
                Zm = K*theta;
                plot([1:length(Z)]',Zm,'. -',[1:length(Z)]',Z,'o-');
                xlabel('Index (-)'); ylabel('Model and Measured')
                legend('Model','Data')
            end
            results.paraVar = diag(results.paraCov);                % Parameter variance
            results.resVec = R;                                     % Residual vector
            results.resNorm = norm(R);                              % Residual norm
        end

        function ct = CovTheta(obj,s,J)
            % Numerically stable calculation of the parameter covaraince matrix.
            %
            % Inputs
            %   s: Sum of squared residuals/(nDataPts - nPara) or 1, if CF is not weighted or weighted respectively
            %   J: Jacobian matrix
            %
            % Outputs
            %   ct: Parameter covariance matrix

            arguments
                obj,s,J
            end

            % Normalise J with the 1/norm of each column
            Jnorm = sqrt(sum(abs(J).^2));
            idxZeros = Jnorm<1E-14;
            Jnorm(idxZeros) = 1;
            N = diag(1./Jnorm);
            Jn = J*N;

            [~, S, V] = svd(Jn,0);          % Perform SVD
            ss = diag(S);                   % Singular values
            idxZeros = ss < 1E-14;          % Replace any small singular value with inf
            nCol = size(Jn,2);
            if sum(idxZeros)>0
                regRank.msg = sprintf('\nRegressor is rank defficient. Rank = %d instead of %d. \nParameters estimated from a subspace.',nCol-sum(idxZeros),nCol);
                regRank.unique = 0;
            else
                regRank.msg = sprintf('Estimated parameters are unique.\nRegressor preserves full rank. Rank =  %d.', nCol);
                regRank.unique = 1;
            end

            ss(idxZeros) = inf;
            Sv = diag(1./ss);  %Inverse singular value matrix

            ct = (N')*V*(Sv)*Sv*(V')*N*s;    % Parameter covariance matrix
        end

        function J = JacApprox(obj,theta,y,u,fh,epsilon)
            % Approximate Jacobian with a first order finite difference
            %
            % Inputs:
            %   theta : Paramater vector. Size nTheta x 1
            %   y: Model output evaluated at theta. Size nData x 1
            %   fh: Function handle
            %   p: Structure to epsilon and input data to simulate model function
            %
            % Outputs:
            %   J: Approximated Jacobian size nData x nTheta

            arguments
                obj,theta,y,u,fh,epsilon
            end

            nTheta = length(theta);
            nData = length(y);
            J = zeros(nData,nTheta);

            deltaTheta = epsilon*theta;                     % Multiply vector of all the parameters by epsilon
            for nn = 1:nTheta
                eVec = zeros(nTheta,1);
                eVec(nn) = 1;
                theta_inc = theta + deltaTheta(nn)*eVec;    % Select the increment for each parameter
                yInc = fh(theta_inc,u);
                J(:,nn) = (yInc-y)./deltaTheta(nn);
            end
        end

        function obj = dataInterPolate(obj)
            arguments
                obj EntropyCoeffEstimator
            end
            try
                Ts = seconds(obj.refSig.Ts_s);
            catch
                error("Please import reference data before importing experimental data")
            end

           obj.measData.tempData = retime(obj.measData.tempData,'regular','linear',TimeStep = Ts);
            
        end
    end



end