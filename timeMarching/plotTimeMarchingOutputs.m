%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2008 - 2018
%
% Sergio Ricci (sergio.ricci@polimi.it)
%
% Politecnico di Milano, Dipartimento di Ingegneria Aerospaziale
% Via La Masa 34, 20156 Milano - ITALY
%
% This file is part of NeoCASS Software (www.neocass.org)
%
% NeoCASS is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation;
% either version 2, or (at your option) any later version.
%
% NeoCASS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public
% License along with NeoCASS; see the file GNU GENERAL
% PUBLIC LICENSE.TXT.  If not, write to the Free Software
% Foundation, 59 Temple Place -Suite 330, Boston, MA
% 02111-1307, USA.
%
%
%
%*******************************************************************************
%
%  NeoCASS
%  Next generation Conceptual Aero Structural Sizing
%
%  SMARTCAD
%  Simplified Models for Aeroelasticity in Conceptual Aircraft Design
%
%                      Sergio Ricci         <ricci@aero.polimi.it>
%                      Luca Cavagna         <cavagna@aero.polimi.it>
%                      Alessandro Degaspari <degaspari@aero.polimi.it>
%                      Luca Riccobene       <riccobene@aero.polimi.it>
%                      Federico Fonte       <federico.fonte@polimi.it>
%                      Francesco Toffol     <francesco.toffol@polimi.it>
%                      Nicola Fonzi         <nicola.fonzi@polimi.it>
%
%
%  Department of Aerospace Sciences and Technologies (DAST)
%  Politecnico di Milano
%
%*******************************************************************************
%
% MODIFICATIONS:
%     DATE        VERS    PROGRAMMER       DESCRIPTION
%     30-09-2021          Nicola Fonzi     Creation
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
ylabel("Amplitude [deg]","FontSize",24)
saveas(gcf,"FFT.fig")
close

return