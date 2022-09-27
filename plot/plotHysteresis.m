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
% You are not entitled to use, distribute, or modify this file in any way,     *
% unless explicitly authorized by the copyright owner.                         *
%                                                                              *
%*******************************************************************************
function plotHysteresis(x,y,initialColor,dashed)
% Very simple function that plots the hysteresis
% It assumes that, if more that one set of data is to be provided, they are
% organised so that the columns of y are in the same number as the columns
% of x

if nargin==2
    initialColor = 1;
    dashed = 0;
elseif nargin==3
    dashed = 0;
end

while initialColor>15
    initialColor = initialColor - 15;
end

colorTable = [ ...
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0.6350    0.0780    0.1840
    0.0000    0.0000    1.0000
    1.0000    0.0000    0.0000
    0.0000    1.0000    0.0000
    0.6000    0.6000    0.6000
    0.0000    0.0000    0.0000
    0.5000    0.0000    0.8000
    0.0000    0.4000    0.4000
    ];

% Treat row vectors
if size(x,1)==1
    x = x.';
end
if size(y,1)==1
    y = y.';
end

if any(y)
    minimum = min(min(y));
    maximum = max(max(y));
    ylim([minimum - 0.1*abs(minimum),maximum + 0.1*abs(maximum)])
end
if any(x)
    minimum = min(min(x));
    maximum = max(max(x));
    xlim([minimum - 0.1*abs(minimum),maximum + 0.1*abs(maximum)])
end

for j = 1:size(y,2)
    for i = 1:size(y,1)-1
        if isnan(x(i))==0 && isnan(y(i,j))==0 && isnan(x(i+1))==0 && isnan(y(i+1,j))==0
            ah = annotation('arrow','headStyle','cback1','HeadLength',8,'HeadWidth',5);
            set(ah,'parent',gca)
            set(ah,'position',[x(i),y(i,j),x(i+1)-x(i),y(i+1,j)-y(i,j)])
            set(ah,'color',colorTable(j+(initialColor-1),:))
            set(ah,'LineWidth',1.5)
            if dashed
                set(ah,'LineStyle','--')
            end
        end
    end
    hold on
    plot([0 0],[0 0],'color',colorTable(j+(initialColor-1),:),'LineWidth',1.5)
end

return