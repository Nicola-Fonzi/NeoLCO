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
function [outputGap, outputFrequency, outputTorque, outputMonitor] = computeTimeMarchingOutputs(rotation, time, torque, modes, Umonitor, options)

type = options.amplitudeDefinition;

outputGap = switchOnType(type,rotation.');

outputTorque = switchOnType(type,torque.');

[outputFrequency.fVect, outputFrequency.pVect] = createFFTfromVariableTimeStepVector(time, rotation, options.timeStepForFFT, options.maxFrequencyInterest);
% We want everything in the format so that the row is the
% point, the column the value
outputFrequency.fVect = outputFrequency.fVect.';
outputFrequency.pVect = outputFrequency.pVect.';

if isempty(Umonitor)
    outputMonitor = [];
else
    outputMonitor = switchOnType(type,Umonitor*modes.');
end

return

function output = switchOnType(type,values)

if strcmp(type,'maxPeak')
    output = findPeaksMatrix(values);
elseif strcmp(type,'rms')
    output = findRmsMatrix(values);
elseif strcmp(type,'std')
    output(:,1) = std(values,[],2);
    output(:,2) = -output(:,1);
    output(:,3) = output(:,1);
else
    error("No type of amplitude definition specified")
end

return

function peaks = findPeaksMatrix(values)

% The matlab function findPeaks is not able to work with matrices, this is
% just a workaround. The matrix will have the format:
% matrix = nonlinearity1 [value1 value2 value3;
%          nonlinearity2  value1 value2 value3];
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

function Rms = findRmsMatrix(values)

Rms=zeros(size(values,1),3);
for i=1:size(values,1)
    bias = mean(values(i,:),2);
    Rms(i,1) = bias + std(values(i,:),[],2);
    Rms(i,2) = bias - std(values(i,:),[],2);
    Rms(i,3) = std(values(i,:),[],2);
end
return