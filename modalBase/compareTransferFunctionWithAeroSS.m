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
function [] = compareTransferFunctionWithAeroSS(modelForIntegration,Ha,k)

syms s
reducedFreqRatio = 1/(modelForIntegration.lref);

s = 1i*k*reducedFreqRatio;
for indexs = 1:length(s)
    H(:,:,indexs) = modelForIntegration.Ca*((s(indexs)/reducedFreqRatio*eye(size(modelForIntegration.A))...
        -modelForIntegration.A)\(modelForIntegration.B0 + modelForIntegration.B1*s(indexs)/reducedFreqRatio ...
        + modelForIntegration.B2*s(indexs)^2/reducedFreqRatio^2 ))...
        + modelForIntegration.D2*s(indexs)^2/reducedFreqRatio^2 + modelForIntegration.D1*s(indexs)...
        /reducedFreqRatio + modelForIntegration.D0;
end
for indexinput = 1:size(H,2)
    for indexoutput = 1:size(H,1)
        for i =1:length(s)
            toPlot(i) = H(indexoutput,indexinput,i);
            temp(i) = Ha(indexoutput,indexinput,i);
        end
        figure
        plot(imag(s),imag(toPlot))
        hold on
        plot(imag(s),imag(temp))
        figure
        plot(imag(s),real(toPlot))
        hold on
        plot(imag(s),real(temp))
    end
end

return