%*******************************************************************************
% Copyright (C) 2020 - 2021                                                    *
%                                                                              *
% Nicola Fonzi (nicola.fonzi@polimi.it)                                        *
%                                                                              *
% Politecnico di Milano, Dipartimento di Ingegneria Aerospaziale               *
% Via La Masa 34, 20156 Milano - ITALY                                         *
%                                                                              *
% This file is part of NeoLCO Software (github.com/Nicola-Fonzi/NeoLCO)        *
%                                                                              *
%*******************************************************************************
%                                                                              *
%                                                                              *
%                                                                              *
% Version: 2.0.0                                                               *
%                                                                              *
%                                                                              *
%                                                                              *
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