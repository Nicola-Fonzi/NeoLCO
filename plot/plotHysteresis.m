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
function plotHysteresis(x,y,varargin)
% Very simple function that plots the hysteresis
% It assumes that, if more that one set of data is to be provided, they are
% organised so that the columns of y are in the same number as the columns
% of x

for i = 1:size(y,1)
    quiver(x(1:end-1),y(i,1:end-1),diff(x),diff(y(i,:)),0,varargin{:},'MaxHeadSize',0.01*max(y(i,:)));
end

return