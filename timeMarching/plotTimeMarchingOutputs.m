%*******************************************************************************
%                                                                              *
%                    _   _            _     ____ ___                           *
%                   | \ | | ___  ___ | |   / ___/ _ \                          *
%                   |  \| |/ _ \/ _ \| |  | |  | | | |                         *
%                   | |\  |  __/ (_) | |__| |__| |_| |                         *
%                   |_| \_|\___|\___/|_____\____\___/                          *
%                                                                              *
%                                                                              *
% Copyright (C) 2020 - 2024                                                    *
%                                                                              *
% Nicola Fonzi (nicola.fonzi@outlook.com)                                      *
%                                                                              *
%                                                                              *
% This file is part of NeoLCO Software (github.com/Nicola-Fonzi/NeoLCO).       *
% The use of this software is licensed based on the licence distributed        *
% together with the source code. If you have not received the license please   *
% contact the copywright owner before using the software.                      *
%                                                                              *
%*******************************************************************************
function plotTimeMarchingOutputs(time,rotation,torque,modes,pdyn,processedFrequency,modelIntegration,options)

figure('name','Amplitude at gap points')
plot(time,rotation,'LineWidth',1.5)
xlabel("Time [s]",'FontSize',16)
ylabel("Amplitude at the nonlinear elements",'FontSize',16)
legend(options.gapPoints{:,4},'FontSize',16)
saveas(gcf,"amplitude.fig")
close

figure('name','Forces at the gap points')
plot(time,torque,'LineWidth',1.5)
xlabel("Time [s]",'FontSize',16)
ylabel("Torque/force at the nonlinear elements",'FontSize',16)
saveas(gcf,"torque.fig")
close

figure('name','Wind speed')
plot(time,sqrt(pdyn*2/modelIntegration.rho),'LineWidth',1.5)
xlabel("Time [s]",'FontSize',16)
ylabel("Speed [m/s]",'FontSize',16)
saveas(gcf,"speed.fig")
close

% Custom graphical outputs: monitor points
if ~isempty(options.monitorPoints)
    figure('name','Amplitude at monitor elements')
    plot(time,(modelIntegration.Umonitor*modes.').','LineWidth',1.5)
    xlabel("Time [s]",'FontSize',16)
    ylabel("Amplitude at the monitor elements",'FontSize',16)
    legend(options.monitorPoints{:,4},'FontSize',16)
    saveas(gcf,"amplitudeMonitor.fig")
    close
end

figure('name','Main frequency of LCO')
plot(processedFrequency.fVect,processedFrequency.pVect,'LineWidth',1.5)
xlabel("Frequency [Hz]","FontSize",24)
ylabel("Amplitude","FontSize",24)
saveas(gcf,"FFT.fig")
close

return