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
function [f,P1] = createFFTfromVariableTimeStepVector(time,value,dt,fMax)

time = seconds(time-time(1));
TT1 = timetable(time,value);
TT2 = retime(TT1,'regular','linear','TimeStep',seconds(dt));


Fs = 1/dt;
L = size(TT2.value,1);

Y = fft(TT2.value);
P2 = abs(Y/L);
P1 = P2(1:L/2+1,:);
P1(2:end-1,:) = 2*P1(2:end-1,:);
f = Fs*(0:(L/2))/L;

fMaxIndex = find(f>=fMax,1);

f=f(1:fMaxIndex);
P1=P1(1:fMaxIndex,:);

return