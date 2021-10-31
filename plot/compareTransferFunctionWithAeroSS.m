%*******************************************************************************
%                                                                              *
%                    _   _            _     ____ ___                           *
%                   | \ | | ___  ___ | |   / ___/ _ \                          *
%                   |  \| |/ _ \/ _ \| |  | |  | | | |                         *
%                   | |\  |  __/ (_) | |__| |__| |_| |                         *
%                   |_| \_|\___|\___/|_____\____\___/                          *
%                                                                              *
%                                                                              *
% Copyright (C) 2020 - 2021                                                    *
%                                                                              *
% Nicola Fonzi (nicola.fonzi@polimi.it)                                        *
%                                                                              *
% Politecnico di Milano, Dipartimento di Ingegneria Aerospaziale               *
% Via La Masa 34, 20156 Milano - ITALY                                         *
%                                                                              *
% This file is part of NeoLCO Software (github.com/Nicola-Fonzi/NeoLCO).       *
% You are not authorized to use, distribute, or modify this file in any way,   *
% unless explicitly decleared otherwise by the copyright owner.                *
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
        set(gcf,"name",strcat("Transfer function with input ",num2str(indexinput)," and output ",num2str(indexoutput)))
        plot(imag(s),imag(toPlot))
        hold on
        plot(imag(s),imag(temp))
        xlabel("Reduced frequency divided by reference length")
        ylabel("Imaginary part")
        figure
        set(gcf,"name",strcat("Transfer function with input ",num2str(indexinput)," and output ",num2str(indexoutput)))
        plot(imag(s),real(toPlot))
        hold on
        plot(imag(s),real(temp))
        xlabel("Reduced frequency divided by reference length")
        ylabel("Real part")
    end
end

return