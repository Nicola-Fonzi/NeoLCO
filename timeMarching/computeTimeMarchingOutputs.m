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
function [outputGap, outputFrequency, outputTorque, outputMonitor] = computeTimeMarchingOutputs(rotation, time, torque, modes, Umonitor, options)

type = options.amplitudeDefinition;

outputGap = switchOnType(type,rotation.');

outputTorque = switchOnType(type,torque.');

[outputFrequency.fVect, outputFrequency.pVect] = createFFTfromVariableTimeStepVector(time, rotation, options.timeStepForFFT, options.maxFrequencyInterest);
% We want everything in the format so that the row is the
% point, the column the value
outputFrequency.fVect = outputFrequency.fVect.';
outputFrequency.pVect = outputFrequency.pVect.';


outputMonitor = switchOnType(type,Umonitor*modes.');

return

function output = switchOnType(type,values)

if strcmp(type,'maxPeak')
    output = findPeaksMatrix(values);
elseif strcmp(type,'rms')
    output(1) = mean(values(values>0),2);
    output(2) = mean(values(values>=0),2);
    output(3) = rms(values,2);
elseif strcmp(type,'std')
    output(1) = std(values,[],2);
    output(2) = output(1);
    output(3) = output(1);
else
    error("No type of amplitude definition specified")
end

return

function peaks = findPeaksMatrix(values)

% The matlab function findPeaks is not able to work with matrices, this is
% just a workaround. The matrix will have the format:
% matrix = nonlinearity1 [value1 value2 value3;
%          nonlinearity2 [value1 value2 value3];
peaks=zeros(size(values,1),3);
for i=1:size(values,1)
    [maxPeaks, ~] = findpeaks(values(i,:));
    [minPeaks, ~] = findpeaks(-values(i,:));
    meanMax = mean(maxPeaks);
    meanMin = -mean(minPeaks);
    peaks(i,1) = meanMax;
    peaks(i,2) = meanMin;
    peaks(i,3) = (meanMax-meanMin)/2;
end

return
