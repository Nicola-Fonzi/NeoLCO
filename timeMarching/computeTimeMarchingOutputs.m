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
