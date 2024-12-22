%% A script to create the paper figures and tables
%
% W.D. Widanage 21/12/2024 (Somewhere in Germany while listening to Christmas songs)


clc
clear
close all

import ECEstimator.*

%% Figure 1: Reversible heat plot

load LGM50_RateTests.mat

cRate = '0p1C';
temp = 'T25';

figure()
plot(LGM50_5Ah_RateTest.(temp).(['cRate_',cRate]).timeVec/3600,LGM50_5Ah_RateTest.(temp).(['cRate_',cRate]).currVec,'. -')
xlabel("Time [H]"); ylabel("Applied current [A]"); grid on;
savefig(gcf,fullfile(pwd,'LGM50_Reversible_Heat_Current.fig'))


figure()
plot(LGM50_5Ah_RateTest.(temp).(['cRate_',cRate]).timeVec/3600,LGM50_5Ah_RateTest.(temp).(['cRate_',cRate]).cellTemp_mid,'. -')
xlabel("Time [H]"); ylabel({"Cell surface", "temperature [degC]"}); grid on;
savefig(gcf,fullfile(pwd,'LGM50_Reversible_Heat.fig'))


%% Figure 2: Time signal and frequency domain

% Generate a reference signal from 0degC to 50degC with 1 period and save 
obj = EntropyCoeffEstimator("Tmin_degC",10,"Tmax_degC",50,"numLevels",5,"startTemp",30,"fMax_Hz",1.2E-3,"Ts_s",2,"Periods",1,"fileName","refSig_1050_60per.csv");

% Frequency domain
refHarm = 1:floor(obj.refSig.Nref/2);
refFreqVec = refHarm/(obj.refSig.Nref*obj.refSig.Ts_s);
UTmp = fft(obj.refSig.refTempSig);

% Excited components
UrefHarm = UTmp(refHarm+1);
UrefExc = UTmp(obj.refSig.excHarmonics+1);

% Suppressed components
suppHarmMul = obj.refSig.Hsupp; % Suppressed harmonic multiples
suppHarmsTmp  = [];
for hh = 1:length(suppHarmMul)
    suppHarmsTmp = [suppHarmsTmp ,suppHarmMul(hh):suppHarmMul(hh): obj.refSig.excHarmonics(end)];
end
suppHarms = sort(suppHarmsTmp); % Sort the suppressed harmonics in ascending order
suppFreq_Hz = suppHarms/(obj.refSig.Nref*obj.refSig.Ts_s);
UrefSupp = UTmp(suppHarms+1);

figure()
stairs(obj.refSig.refTimeVec_s/3600,obj.refSig.refTempSig,'. -'); grid on;
xlabel('Time [H]'); ylabel({'Reference temperature';'signal [deg$^\circ$C]'})
savefig(gcf,fullfile(pwd,'Reference_Signal.fig'))


figure()
semilogx(obj.refSig.excFreq_Hz*1000,abs(UrefExc),'o r',refFreqVec*1000,abs(UrefHarm),'* b'); hold on;
xlabel('Frequency [mHz]'); ylabel({'FFT magnitude of';'temperature signal [-]'})
if ~isempty(suppHarms)
    semilogx(suppFreq_Hz*1000,abs(UrefSupp),'o g'); grid on, hold on; 
end
savefig(gcf,fullfile(pwd,'Reference_Signal_FFT.fig'))



%% Perform Kernel estimation and dUdT fitting 

dataPth = what('Measurement_Data/measurements_Aug2023').path;
hdrNames = ["time", "TEC1", "TEC2", "BoxTop", "TabAnode", "SurfaceBottomAnode", "SurfaceTopAnode", "SurfaceBottomCathode", "SurfaceTopCathode", "TabCathode", "SurfaceTopCenter", "SurfaceBottomCenter", "CoolingBlockTop", "Ambient", "U"];
freqTextFilesInfo = dir(fullfile(dataPth,"*Frequency.txt"));
z = (0:5:100)'; % SoC break points

% Estimate dUdT over each SoC with auto model order selection
cntr = 0;
for zz = 1:numel(z)
    kerObj(zz,1) = EntropyCoeffEstimator();
    kerObj(zz,1).ImportRefSig("filePth",fullfile(dataPth,"refSig","refSig_1050_July2022.mat")); % Import reference signal
    kerObj(zz,1).ImportExpData("filePth",fullfile(dataPth,freqTextFilesInfo(zz).name),"HdrNames",hdrNames);
    fprintf("Estimating for SoC %d\n",z(zz))
    kerObj(zz,1).EstimateEntropyCoeff("usePeriods",1,"transientOnOff","on");
    
    % Collect GoF, model order and full-rank status
    GoF(zz,1) = kerObj(zz).results.fitMetrics.FitPercent;
    RMSE(zz,1) = kerObj(zz).results.fitMetrics.RMSE;
    model_order(zz,:) = [kerObj(zz).estimationSettings.modelOrder_num,kerObj(zz).estimationSettings.modelOrder_denom]; 
    full_rank(zz,1) = kerObj(zz,1).results.fitMetrics.LMRankFull(end);

    % Collect dUdT and std
    dUdTK(zz,1) = kerObj(zz,1).results.dUdT_mVpK;
    dUdTK_std(zz,1) = kerObj(zz,1).results.dUdT_std;
end

kernel_results_table = table(z,GoF,RMSE,full_rank,model_order);

% Improve fit for low GoFs
soc_select = 10;
idx = find(soc_select == z);
kerObj(idx,1).EstimateEntropyCoeff("usePeriods",1,"transientOnOff","on","modelOrder_num",2,"modelOrder_denom",3,freqIdx_estimation=(1:5));

soc_select = 20;
idx = find(soc_select == z);
kerObj(idx,1).EstimateEntropyCoeff("usePeriods",1,"transientOnOff","on","modelOrder_num",2,"modelOrder_denom",3,freqIdx_estimation=(1:5));

soc_select = 25;
idx = find(soc_select == z);
kerObj(idx,1).EstimateEntropyCoeff("usePeriods",1,"transientOnOff","on","modelOrder_num",2,"modelOrder_denom",3,freqIdx_estimation=(1:5));
kerObj(idx).results.fitMetrics.FitPercent

soc_select = 30;
idx = find(soc_select == z);
kerObj(idx,1).EstimateEntropyCoeff("usePeriods",1,"transientOnOff","on","modelOrder_num",2,"modelOrder_denom",3,freqIdx_estimation=(1:5));

soc_select = 95;
idx = find(soc_select == z);
kerObj(idx,1).EstimateEntropyCoeff("usePeriods",1,"transientOnOff","on","modelOrder_num",2,"modelOrder_denom",3,freqIdx_estimation=(1:5));

% Re-collect kernel based dUdT and the standard deviation and fit metrics
for zz = 1:length(z)
    dUdTK(zz,1) = kerObj(zz,1).results.dUdT_mVpK;
    dUdTK_std(zz,1) = kerObj(zz,1).results.dUdT_std;

    % Collect GoF, model order and full-rank status
    GoF(zz,1) = kerObj(zz).results.fitMetrics.FitPercent;
    RMSE(zz,1) = kerObj(zz).results.fitMetrics.RMSE;
    full_rank(zz,1) = kerObj(zz,1).results.fitMetrics.LMRankFull(end);

    model_order(zz,:) = [kerObj(zz).estimationSettings.modelOrder_num,kerObj(zz).estimationSettings.modelOrder_denom]; 
   
end

kernel_results_table = table(z,GoF,RMSE,full_rank,model_order,dUdTK,dUdTK_std);
head(kernel_results_table)

%% Figure 3: Reference signal, processed temperature and OCV

close all
soc_select = 80;
idx = find(soc_select == z);

ref_time = kerObj(idx).refSig.refTimeVec_s;
ref_temperature_signal = kerObj(idx).refSig.refTempSig;
meas_mean_temp_signal = kerObj(idx).processedData.CaloricPeriod_degC.Caloric_Period;
meas_ocv_signal = kerObj(idx).processedData.ocvPeriod_V.OCV_Period;

figure();
stairs(ref_time/3600,ref_temperature_signal); hold on
stairs(ref_time/3600,meas_mean_temp_signal,'. -')
xlabel("Time [H]"); ylabel("Temperature [degC]"); grid on;
savefig(gcf,fullfile(pwd,'Reference_Processed_Temperature.fig'))

figure
plot(ref_time/3600,meas_ocv_signal,'. -');
xlabel("Time [H]"); ylabel("OCV [V]"); grid on;
savefig(gcf,fullfile(pwd,'Measured_OCV.fig'))

%% Figure 4: Temperature and OCV frequency content
close all
soc_select = 80;
idx = find(soc_select == z);
kerObj(idx).PlotProcessedSig;

% FFT plot of measured mean temperature 
figure(1)
p = gcf().Children(2).Children; % FFT of measured mean temperature
figure
plot(p(1).XData,p(1).YData,'o')
semilogx(p(3).XData,p(3).YData,'o r',p(2).XData,p(2).YData,'* b',p(1).XData,p(1).YData,'o g'); grid on;
xlabel('Freq (mHz)'); ylabel({' Temperature FFT'; 'magnitude [-]'}); 
legend('Excited frequencies','All frequencies', 'Suppressed frequencies')
savefig(gcf,fullfile(pwd,'Measured_Temperature_FFT.fig'))


% FFT plot of measured mean temperature 
figure(2)
p = gcf().Children(2).Children; % FFT of measured mean temperature
figure
plot(p(1).XData,p(1).YData,'o')
semilogx(p(3).XData,p(3).YData,'o r',p(2).XData,p(2).YData,'* b',p(1).XData,p(1).YData,'o g'); grid on;
xlabel('Freq (mHz)'); ylabel({' Temperature FFT'; 'magnitude [-]'}); 
legend('Excited frequencies','All frequencies', 'Suppressed frequencies')
savefig(gcf,fullfile(pwd,'Measured_OCV_FFT.fig'))


%% Figure 5: FRF and TF Fit

close all
soc_select = 80;
idx = find(soc_select == z);

freq_mHz = kerObj(idx).refSig.excFreq_Hz*1000;
magMeas = 20*log10(abs(kerObj(idx).results.kernel));
phaseMeas = unwrap(angle(kerObj(idx).results.kernel));
stdMeas = 10*log10(kerObj(idx).results.varKernel);
magFit = 20*log10(abs(kerObj(idx).results.optimumFRF));
phaseFit = unwrap(angle(kerObj(idx).results.optimumFRF));

figure
yyaxis left
plot(freq_mHz,magMeas,'. -'); grid on
xlabel("Frequency [mHz]"); ylabel("Magnitude [dB]");
yyaxis right
plot(freq_mHz,stdMeas,'. -'); grid on
xlabel("Frequency [mHz]"); ylabel("Magnitude [dB]"); legend(["$|\hat{G}(\omega_k)|$";"$\sigma_{\hat{G}}(\omega_k)$"],Interpreter="latex")
savefig(gcf,fullfile(pwd,'Kernel_Estimate_Mag.fig'))

figure
plot(freq_mHz,phaseMeas,'. -'); grid on
xlabel("Frequency [mHz]"); ylabel("Phase [rad]");
savefig(gcf,fullfile(pwd,'Kernel_Estimate_Phase.fig'))

figure
plot(freq_mHz,magMeas,'. -',freq_mHz,magFit); grid on
xlabel("Frequency [mHz]"); ylabel("Magnitude [dB]"); legend(["$|\hat{G}(\omega_k)|$";"$|G(\omega_k,\theta)|$"],Interpreter="latex")
savefig(gcf,fullfile(pwd,'Kernel_Fit_Mag.fig'))

figure
plot(freq_mHz,phaseMeas,'. -',freq_mHz,phaseFit); grid on
xlabel("Frequency [mHz]"); ylabel("Phase [rad]"); legend(["$\angle \hat{G}(\omega_k)$";"$\angle G(\omega_k,\theta)$"],Interpreter="latex")
savefig(gcf,fullfile(pwd,'Kernel_Fit_Phase.fig'))


%% Figure 6: Potentiometric based method

close all
potTextFilesInfo = dir(fullfile(dataPth,"*Potentiometric.txt"));
hdrNames = ["time", "TEC1", "TEC2", "BoxTop", "TabAnode", "SurfaceBottomAnode", "SurfaceTopAnode", "SurfaceBottomCathode", "SurfaceTopCathode", "TabCathode", "SurfaceTopCenter", "SurfaceBottomCenter", "CoolingBlockTop", "Ambient", "U"];
z = (0:5:100)'; % SoC break points
plot_soc = 80;

for zz = 1:numel(z)
    filePath = fullfile(dataPth,potTextFilesInfo(zz).name);
    imOpts = detectImportOptions(filePath);
    imOpts.VariableNames = hdrNames;
    potData = readtable(filePath,imOpts);
    Time_s = seconds(potData.(1) - potData.(1)(1)); % Reset time to 0 seconds
    potData.Time_s = Time_s;
    potData = table2timetable(potData,"RowTimes",Time_s);
    potDuration(zz,1) = hours(potData.Time(end));
    caloricTemp = CaloricCellTemperature(potData);
    
    % Signals for potentiometric method
    time = caloricTemp.Time;
    temp = caloricTemp.mean_temperature_degC;
    ocv = potData.U;

    % Index set
    idx50 = find(temp > 49.5 & temp < 50.5, 1, 'last');
    idx40 = find(temp > 39.3 & temp < 40.3, 1,'last');
    idx30 = find(temp > 29.01 & temp < 30.01, 1,'last');
    idx20 = find(temp > 19.01 & temp < 20.01, 1,'last');
    idx10 = find(temp > 9.05 & temp < 10.5, 1,'last');

    idxSS = [idx50,idx40,idx30,idx20,idx10];

    ssTemp = temp(idxSS);
    ssOCV = ocv(idxSS);

    % Best fit
    K = [ones(size(ssTemp)),ssTemp];                  % Regressor matrix;
    [dUdTP_Tmp,dudTP_info] = Lls(K,ssOCV);            % Straight line fit
    dUdTFit = polyval(flipud(dUdTP_Tmp),ssTemp);      % [V]

    dUdTP(zz,1) = dUdTP_Tmp(2)*1000;                        % [mV/K]
    dUdTP_std(zz,1) = sqrt(dudTP_info.paraVar(2))*1000;     % [mV/K]

    % Error metrics
    [RMSE_tmp, GoF_tmp] = CalculateErrorMetrics(dUdTFit,ssOCV);
    RMSE(zz,1) = RMSE_tmp;
    GoF(zz,1) = GoF_tmp;

    if (z(zz) == plot_soc)
        figure
        subplot(2,1,1);
        plot(hours(time),potData.U); grid on;
        xlabel("Time [H]"); ylabel("OCV [V]"); title(['SoC: ' num2str(z(zz)) '%'])

        subplot(2,1,2);
        plot(hours(time),temp,'. -'); grid on;
        xlabel("Time [H]"); ylabel("Temperature [degC]")
        savefig(gcf,fullfile(pwd,sprintf('Potentiometric_Signal_%d.fig',z(zz))))

        figure
        plot(ssTemp,ssOCV,'x',ssTemp,dUdTFit,'-'); grid on;
        xlabel("Temperature steady-state [degC]"); ylabel("OCV steadystate [V]"); title(['SoC: ' num2str(z(zz)) '%'])
        savefig(gcf,fullfile(pwd,sprintf('Potentiometric_Fit_%d.fig',z(zz))))

    end
end

pot_results_table = table(z,GoF,RMSE,potDuration,dUdTP,dUdTP_std);
head(pot_results_table)

%% Figure 7: Plot dUdT for Kernel and potentiometric method

close all
figure
[z_sort,idx_sort] = sort(z);
errorbar(z_sort,dUdTK(idx_sort),dUdTK_std(idx_sort))
hold on
errorbar(z_sort,dUdTP(idx_sort),dUdTP_std(idx_sort),[],[],[],LineStyle="- ."); grid on;
xlabel("SoC [%]"); ylabel("dUdT [mV/K]")
legend(["Kernel based","Potentiometirc"],"Location","best")

savefig(fullfile(pwd,'Kernal_Potentiometric_dUdT.fig'))

%% Helper functions
% Function to calculate mean temperature 

function [caloricTemp] = CaloricCellTemperature(potData)
MeanTemp = (potData.SurfaceTopCenter + potData.SurfaceTopCathode + potData.SurfaceTopAnode + potData.SurfaceBottomCenter + potData.SurfaceBottomCathode)/5;
caloricTemp = timetable(potData.Time_s, MeanTemp,'VariableNames',"mean_temperature_degC");
end

% Error metrics
function [RMSE,GoF] = CalculateErrorMetrics(model,meas)
SSE = sum((model - meas).^2);
SST = sum(abs(meas-mean(meas)).^2);
RMSE = sqrt(SSE);
GoF = (1 - SSE/SST)*100;
end