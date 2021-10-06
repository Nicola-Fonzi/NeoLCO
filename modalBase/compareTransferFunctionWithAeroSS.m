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