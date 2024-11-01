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