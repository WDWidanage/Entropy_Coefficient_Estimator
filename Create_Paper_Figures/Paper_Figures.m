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
    kerObj(zz,1).EstimateEntropyCoeff("usePeriods",[1,2],"transientOnOff","on");
    
    % Collect GoF, model order and full-rank status
    GoF(zz,1) = kerObj(zz).results.fitMetrics.FitPercent;
    RMSE(zz,1) = kerObj(zz).results.fitMetrics.RMSE;
    model_order(zz,:) = [kerObj(zz).estimationSettings.modelOrder_num,kerObj(zz).estimationSettings.modelOrder_denom]; 
    full_rank(zz,1) = kerObj(zz,1).results.fitMetrics.LMRankFull(end);

    % Collect dUdT and std
    dUdTK(zz,1) = kerObj(zz,1).results.dUdT_mVpK;
    dUdTK_std(zz,1) = kerObj(zz,1).results.dUdT_std;
end

results_table = table(z,GoF,RMSE,full_rank,model_order);

% Improve fit for low GoFs
close all
soc_select = 20;
idx = find(soc_select == z);
kerObj(idx,1).EstimateEntropyCoeff("usePeriods",1,"transientOnOff","on","modelOrder_num",6,"modelOrder_denom",5);

soc_select = 30;
idx = find(soc_select == z);
kerObj(idx,1).EstimateEntropyCoeff("usePeriods",1,"transientOnOff","on","modelOrder_num",6,"modelOrder_denom",5);


soc_select = 95;
idx = find(soc_select == z);
kerObj(idx,1).EstimateEntropyCoeff("usePeriods",1,"transientOnOff","on","modelOrder_num",2,"modelOrder_denom",4);

% Re-collect kernel based dUdT and the standard deviation and fit metrics
for zz = 1:length(z)
    dUdTK(zz,1) = kerObj(zz,1).results.dUdT_mVpK;
    dUdTK_std(zz,1) = kerObj(zz,1).results.dUdT_std;

    % Collect GoF, model order and full-rank status
    GoF(zz,1) = kerObj(zz).results.fitMetrics.FitPercent;
    RMSE(zz,1) = kerObj(zz).results.fitMetrics.RMSE;
    model_order(zz,:) = [kerObj(zz).estimationSettings.modelOrder_num,kerObj(zz).estimationSettings.modelOrder_denom]; 
   
end

results_table = table(z,GoF,RMSE,model_order);


%% Figure 3: FRF and TF Fit

close all
soc_select = 45;
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


%% Potentiometric based method

potTextFilesInfo = dir(fullfile(dataPth,"*Potentiometric.txt"));
hdrNames = ["time", "TEC1", "TEC2", "BoxTop", "TabAnode", "SurfaceBottomAnode", "SurfaceTopAnode", "SurfaceBottomCathode", "SurfaceTopCathode", "TabCathode", "SurfaceTopCenter", "SurfaceBottomCenter", "CoolingBlockTop", "Ambient", "U"];
z = 0:5:100; % SoC break points
plot_intermediate = false;

for zz = 1:numel(z)
    filePath = fullfile(dataPth,potTextFilesInfo(zz).name);
    imOpts = detectImportOptions(filePath);
    imOpts.VariableNames = hdrNames;
    potData = readtable(filePath,imOpts);
    Time_s = seconds(potData.(1) - potData.(1)(1)); % Reset time to 0 seconds
    potData.Time_s = Time_s;
    potData = table2timetable(potData,"RowTimes",Time_s);
    potDuration = hours(potData.Time(end))
    caloricTemp = CaloricCellTemperature(potData);
    time = caloricTemp.Time;
    temp = caloricTemp.meanCaloric;
    ocv = potData.U;

    idx50 = find(temp > 49.5 & temp < 50.5, 1, 'last');
    idx40 = find(temp > 39.3 & temp < 40.3, 1,'last');
    idx30 = find(temp > 29.01 & temp < 30.01, 1,'last');
    idx20 = find(temp > 19.01 & temp < 20.01, 1,'last');
    idx10 = find(temp > 9.05 & temp < 10.5, 1,'last');

    idxSS = [idx50,idx40,idx30,idx20,idx10];

    ssOCV = ocv(idxSS);
    ssTemp = temp(idxSS);

    tempSelected = ssTemp;
    ocvSelected = ssOCV*1000; % [mV]

    % Best fit
    K = [ones(size(tempSelected)),tempSelected];    % Regressor matrix;
    [dUdTP_Tmp,dudTP_info] = Lls(K,ocvSelected);
    dUdTFit = polyval(flipud(dUdTP_Tmp),ssTemp)/1000;       % [V]

    dUdTP(zz,1) = dUdTP_Tmp(2);
    dUdTP_std(zz,1) = sqrt(dudTP_info.paraVar(2));

    if (plot_intermediate)
        figure
        subplot(2,1,1);
        plot(hours(time),potData.U);
        xlabel("Time [H]"); ylabel("OCV [V]"); title(['SoC: ' num2str(z(zz)) '%'])

        subplot(2,1,2);
        plot(hours(time),temp,'. -',...
            hours(time(idx50)),temp(idx50),'o',...
            hours(time(idx40)),temp(idx40),'o',...
            hours(time(idx30)),temp(idx30),'o',...
            hours(time(idx20)),temp(idx20),'o',...
            hours(time(idx10)),temp(idx10),'o')
        xlabel("Time [H]"); ylabel("Temperature [degC]")

        figure
        plot(ssTemp,ssOCV,'x',ssTemp,dUdTFit,'-')
        xlabel("Temperature steady-state [degC]"); ylabel("OCV steadystate [V]"); title(['SoC: ' num2str(z(zz)) '%'])
    end
end
%% Plot dUdT for Kernel, potentiometric

figure
[z_sort,idx_sort] = sort(z);
errorbar(z_sort,dUdTK(idx_sort),dUdTK_std(idx_sort))
hold on
errorbar(z_sort,dUdTP(idx_sort),dUdTP_std(idx_sort),[],[],[],LineStyle="- ."); grid on;
xlabel("SoC [%]"); ylabel("dUdT [mV/K]")
legend(["Kernel based","Potentiometirc"],"Location","best")

savefig(fullfile(pwd,'Kernal_Potentiometric_dUdT.fig'))
%% Function to calculate caloric temperature 

function [caloricTemp] = CaloricCellTemperature(potData)
p = EntropyCoeffEstimator();

alphaCuCell = p.cellThermalProperties.copperConductivity_wpmk/p.cellThermalProperties.copperThickness_m; % Approximation heat conductivity Cu-plate to Cell, W/(m^2*K)
KeyFigures.Nu2Infinity = (pi^2)/2;
KeyFigures.Bi = (alphaCuCell * p.cellThermalProperties.cellThickness_m)/ p.cellThermalProperties.cellConductivity_wpmk; % Biot-Number
KeyFigures.NuiInfinity = (4+p.cellThermalProperties.GeoFactor + KeyFigures.Bi)/(1+(KeyFigures.Bi/KeyFigures.Nu2Infinity));


MeanTemp.TECRef = (potData.TEC1 + potData.TEC2)/2;
MeanTemp.Caloric = (potData.SurfaceTopCenter + potData.SurfaceTopCathode + potData.SurfaceTopAnode + potData.SurfaceBottomCenter + potData.SurfaceBottomCathode)/5; %for erste Zeile kann wegen t=0 keine kalorische Mitteltemperatur berechnet werden
setTemp = MeanTemp.TECRef;

Zeilen = size(MeanTemp.TECRef,1);
t = 2;
Caloric(1,1) = 0;
for ii = 2:Zeilen
    if isequal(setTemp(ii), setTemp(ii-1))
        %Kalorische Mitteltemperatur
        KeyFigures.Fo(ii) = (p.cellThermalProperties.cellConductivity_wpmk*t)/(p.cellThermalProperties.cellThickness_m^2 * p.cellThermalProperties.cellDenisty_kgpm3 * p.cellThermalProperties.cellSpecificHeatCapacity_Jp);
        KeyFigures.Nui0(ii) = (sqrt(pi)+ 10*KeyFigures.Bi*sqrt(KeyFigures.Fo(ii)))/(1+5*KeyFigures.Bi*sqrt(pi*KeyFigures.Fo(ii)))*1/sqrt(KeyFigures.Fo(ii));
        KeyFigures.Nui(ii) = sqrt((KeyFigures.NuiInfinity^2)-(0.4^2)+(KeyFigures.Nui0(ii) + 0.4)^2);
        KeyFigures.NTUi(ii) = (p.cellThermalProperties.GeoFactor * KeyFigures.Fo(ii))/((1/KeyFigures.Bi) + (1/KeyFigures.Nui(ii)));
        Caloric(ii,1) = (MeanTemp.TECRef(ii) + (MeanTemp.Caloric(ii-(t/2))-MeanTemp.TECRef(ii))*exp(-KeyFigures.NTUi(ii)));
        t = t+2;
    else
        t = 2;
        %Kalorische Mitteltemperatur
        KeyFigures.Fo(ii) = (p.cellThermalProperties.cellConductivity_wpmk*t)/(p.cellThermalProperties.cellThickness_m^2 * p.cellThermalProperties.cellDenisty_kgpm3 * p.cellThermalProperties.cellSpecificHeatCapacity_Jp);
        KeyFigures.Nui0(ii) = ((sqrt(pi)+ 10*KeyFigures.Bi*sqrt(KeyFigures.Fo(ii)))/(1+5*KeyFigures.Bi*sqrt(pi*KeyFigures.Fo(ii))))*(1./sqrt(KeyFigures.Fo(ii)));
        KeyFigures.Nui(ii) = sqrt(KeyFigures.NuiInfinity^2-0.4^2.+(KeyFigures.Nui0(ii) + 0.4)^2);
        KeyFigures.NTUi(ii) = (p.cellThermalProperties.GeoFactor * KeyFigures.Fo(ii))/((1/KeyFigures.Bi) + (1/KeyFigures.Nui(ii)));
        Caloric(ii,1) =(MeanTemp.TECRef(ii) + (MeanTemp.Caloric(ii-1) - MeanTemp.TECRef(ii))*exp(-KeyFigures.NTUi(ii)));
    end

end
caloricTemp = timetable(potData.Time_s,Caloric, MeanTemp.Caloric, MeanTemp.TECRef,'VariableNames',["dynamicCaloric","meanCaloric","meanTecRef"]);
end