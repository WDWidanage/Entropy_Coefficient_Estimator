% Prepare figures for paper
close all
clc
%% Signals

obj = load("ecObj_1050.mat").ecObj_1050;

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

figure(1)
axes(ColorOrder=brewermap([],'Set2'),NextPlot = 'replacechildren');
stairs(obj.refSig.refTimeVec_s/3600,obj.refSig.refTempSig,'. -'); grid on;
xlabel('Time [H]'); ylabel({'Reference temperature';'signal [deg$^\circ$C]'})
PrepareFigure(1,'LineWidth',2,'pdf','fileName','signal-time','axisFontSize',24,'labelFontSize',30,'annotationFontSize',30')

figure(2)
semilogx(obj.refSig.excFreq_Hz*1000,abs(UrefExc),'o r',refFreqVec*1000,abs(UrefHarm),'* b'); hold on;
xlabel('Frequency [mHz]'); ylabel({'FFT magnitude of';'temperature signal [-]'})
% legend('Excited frequencies','All frequencies')
if ~isempty(suppHarms)
    semilogx(suppFreq_Hz*1000,abs(UrefSupp),'o g'); hold on; 
end
annotation(gcf,'textarrow',[0.331818181818182 0.307272727272727],...
    [0.288333333333333 0.223333333333333],...
    'String',{'Suppressed even','frequencies'});
annotation(gcf,'textarrow',[0.464545454545454 0.44],...
    [0.728333333333337 0.663333333333337],...
    'String',{'Excited odd','frequencies'});
annotation(gcf,'textarrow',[0.571818181818182 0.547272727272727],...
    [0.398333333333337 0.333333333333337],'String',{'All other','frequencies'});

PrepareFigure(2,'LineWidth',2,'pdf','fileName','signal-fft','axisFontSize',24,'labelFontSize',30,'annotationFontSize',30')


%% Reversible heat
clear
close all
clc

load LGM50_RateTests.mat

cRate = '0p1C';
temp = 'T25';

figure(3)
plot(LGM50_5Ah_RateTest.(temp).(['cRate_',cRate]).timeVec/3600,LGM50_5Ah_RateTest.(temp).(['cRate_',cRate]).currVec,'. -')
xlabel("Time [H]"); ylabel("Applied current [A]")
annotation(gcf,'textarrow',[0.294545454545455 0.316363636363636],...
    [0.424 0.328333333333333],'String',{'0.1C discharge'});
annotation(gcf,'textarrow',[0.656363636363636 0.678181818181818],...
    [0.892333333333335 0.796666666666668],'String',{'C/3 CCCV charge'});
PrepareFigure(3,'LineWidth',2,'pdf','fileName','discharge-charge-current','axisFontSize',24,'labelFontSize',30,'annotationFontSize',30)

figure(4)
plot(LGM50_5Ah_RateTest.(temp).(['cRate_',cRate]).timeVec/3600,LGM50_5Ah_RateTest.(temp).(['cRate_',cRate]).cellTemp_mid,'. -')
xlabel("Time [H]"); ylabel({"Cell surface", "temperature [degC]"})
annotation(gcf,'textarrow',[0.38 0.349090909090909],...
    [0.348333333333333 0.245],'String',{'Reversible heat'});
annotation(gcf,'textarrow',[0.627272727272727 0.671818181818182],...
    [0.745 0.59],'TextLineWidth',1,'String',{'Reversible heat'});
PrepareFigure(4,'LineWidth',2,'pdf','fileName','reversible-heat','axisFontSize',24,'labelFontSize',30,'annotationFontSize',30)

%% Kernel FRF estimate
clear
close all

open("Kernel_Estimate_Mag.fig")
legend("off")
PrepareFigure(gcf,'LineWidth',2,'pdf','fileName','Kernel_Estimate_Mag','axisFontSize',20,'labelFontSize',20,'markerSize',24,'noSave')
legend(Interpreter="latex")
exportgraphics(gcf,'Kernel_Estimate_Mag.pdf')

open("Kernel_Estimate_Phase.fig")
legend("off")
PrepareFigure(gcf,'LineWidth',2,'pdf','fileName','Kernel_Estimate_Phase','axisFontSize',20,'labelFontSize',20,'markerSize',24)



open("Kernel_Fit_Mag.fig")
legend("off")
annotation(gcf,'textarrow',[0.446363636363636 0.492727272727273],...
    [0.316666666666668 0.441666666666668],'String',{'Estimated FRF'});
annotation(gcf,'textarrow',[0.569090909090909 0.55],...
    [0.500000000000001 0.423333333333335],'String',{'Transfer function','fit'});
PrepareFigure(3,'LineWidth',2,'pdf','fileName','Kernel_Fit_Mag','axisFontSize',20,'labelFontSize',20,'markerSize',24,'annotationFontSize',20)
legend(Interpreter="latex")
exportgraphics(gcf,'Kernal_Fit_Mag.pdf')


open("Kernel_Fit_Phase.fig")
legend("off")
annotation(gcf,'textarrow',[0.565454545454546 0.541818181818182],...
    [0.383333333333337 0.308333333333341],'String',{'Estimated FRF phase'});
annotation(gcf,'textarrow',[0.4 0.380909090909091],...
    [0.445000000000004 0.368333333333338],'String',{'Transfer function','fit'});
PrepareFigure(gcf,'LineWidth',2,'pdf','fileName','Kernel_Fit_Phase','axisFontSize',20,'labelFontSize',20,'markerSize',24,'annotationFontSize',20)
legend(Interpreter="latex")
exportgraphics(gcf,'Kernel_Fit_Phase.pdf')
%% dUdT
clear
close all
open("Kernel_Potentio_dUdT.fig")
ylabel("Entropy coefficient [mV/K]")
xlabel("SoC [\%]")

legend("off")
annotation(gcf,'textarrow',[0.742727272727273 0.717272727272727],...
    [0.623333333333334 0.736666666666667],'String',{'Kernel based method'});
annotation(gcf,'textarrow',[0.591818181818182 0.614545454545455],...
    [0.771666666666667 0.685000000000003],'String',{'Potentiometric method'});
PrepareFigure(gcf,'LineWidth',2,'pdf','fileName','Kernel_Potentiometric_dUdT','axisFontSize',20,'labelFontSize',20,'annotationFontSize',20)
