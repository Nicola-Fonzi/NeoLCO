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
%     29-06-2021          Nicola Fonzi     Creation
%*******************************************************************************
function [] = plotMaxRealEigSS(modelForIntegration)


V_vect = linspace(0.001,160,300);
Pdyn = V_vect.^2*modelForIntegration.rho*0.5;

for i = 1:length(Pdyn)
    reducedFreqRatio = V_vect(i)/(modelForIntegration.lref);
    Mtil=modelForIntegration.Mhh-modelForIntegration.D2/reducedFreqRatio^2*Pdyn(i);
    Ktil=modelForIntegration.Khh-modelForIntegration.D0*Pdyn(i);
    Ktil_stiff=modelForIntegration.Khh_stiff-modelForIntegration.D0*Pdyn(i);
    Ctil=modelForIntegration.Chh-modelForIntegration.D1/reducedFreqRatio*Pdyn(i);
    invMKtil = inv(Mtil)*Ktil;
    invMKtil_stiff = inv(Mtil)*Ktil_stiff;
    invMCtil = inv(Mtil)*Ctil;
    invMctil = inv(Mtil)*modelForIntegration.Ca*Pdyn(i);
    Amatrix = [zeros(modelForIntegration.nStru) eye(modelForIntegration.nStru) ...
        zeros(modelForIntegration.nStru,size(modelForIntegration.A,1));...
        -invMKtil -invMCtil invMctil;...
        reducedFreqRatio*modelForIntegration.B0-modelForIntegration.B2/reducedFreqRatio...
        *invMKtil modelForIntegration.B1-modelForIntegration.B2/reducedFreqRatio*invMCtil...
        reducedFreqRatio*modelForIntegration.A+modelForIntegration.B2/reducedFreqRatio*invMctil];
    Amatrix_stiff = [zeros(modelForIntegration.nStru) eye(modelForIntegration.nStru) ...
        zeros(modelForIntegration.nStru,size(modelForIntegration.A,1));...
        -invMKtil_stiff -invMCtil invMctil;...
        reducedFreqRatio*modelForIntegration.B0-modelForIntegration.B2/reducedFreqRatio...
        *invMKtil_stiff modelForIntegration.B1-modelForIntegration.B2/reducedFreqRatio*invMCtil...
        reducedFreqRatio*modelForIntegration.A+modelForIntegration.B2/reducedFreqRatio*invMctil];
    lambda = eig(Amatrix);
    lambda_stiff = eig(Amatrix_stiff);
    maxEig(i) = max(real(lambda));
    maxEig_stiff(i) = max(real(lambda_stiff));
end


figure
plot(V_vect,maxEig)
xlabel("Flow speed [m/s]")
ylabel("Max real part of eigenvalues")
title("Free system")
h = gcf;
saveas(h,"maxRealEig.fig")
close

figure
plot(V_vect,maxEig_stiff)
xlabel("Flow speed [m/s]")
ylabel("Max real part of eigenvalues")
title("Stiff system")
h = gcf;
saveas(h,"maxRealEig_stiff.fig")
close

return