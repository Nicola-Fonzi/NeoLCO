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
function [t,q,qdot,qddot] = readHistoryModal(filename,nmodes,n,display)

if nargin == 2
    n = false;
end

if nargin == 3
    display = false;
end

if ~islogical(display)
    error('display must be logical');
end

fid = fopen(filename);

% Read the formatted file StructHistoryModal.dat. The first three columns
% are always present. Other 3 are added per each mode
formatSpec = '%f%f%f';
for i = 1:nmodes
    formatSpec = [formatSpec, '%f%f%f'];
end
data = textscan(fid,formatSpec,'HeaderLines',1);
data = cell2mat(data);

fclose(fid);

t = data(:,1);


if ~n
    n = str2double(input('n: ','s'));
end
if n > nmodes
    error(['Mode not found. NMODES = ',num2str(nmodes)]);
end

i = 3+(n-1)*3;
q = data(:,i+1);
qdot = data(:,i+2);
qddot = data(:,i+3);

if display
    figure();
    hold on;
    subplot(3,1,1);
    plot(t,q);
    xlabel('t');
    ylabel('$q$','interpreter','latex','FontSize',14);
    title(['Mode n. ', num2str(n)]);
    subplot(3,1,2);
    plot(t,qdot);
    xlabel('t');
    ylabel('$\dot{q}$','interpreter','latex','FontSize',14);
    subplot(3,1,3);
    plot(t,qddot);
    xlabel('t');
    ylabel('$\ddot{q}$','interpreter','latex','FontSize',14);
end

end

