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
function [ID,GridType,U,Ux,Uy,Uz,K,M,Uxr,Uyr,Uzr,Usp] = readPunchShapes(filename)

fid = fopen(filename);

% Read a punch file and extract the mode shapes
% A new mode starts with the keywords EIGENVALUE and MODE
% Treat differently the scalar points, as they only have 1 degree of
% freedom
justRead = 0;
while ~feof(fid)
    if ~justRead
        card = fgetl(fid);
    else
        justRead = 0;
    end
    if contains(card,'EIGENVALUE') && contains(card,'MODE')
        j = str2double(card(38:42));
        k = str2double(card(16:28));
        M(j,j) = 1;
        K(j,j) = k;
        card = fgetl(fid);
        keepGoing = 1;
        i = 1;
        isp = 1 ;
        index = 1;
        while keepGoing
            type = card(16:20);
            GID = str2double(card(1:10));
            if type == '  G  '
                GridID(index:index+5) = GID;
                GridType(index:index+5) = type(3);
                v = str2num(card(24:72));
                ux = v(1);
                uy = v(2);
                uz = v(3);
                Ux(i,j) = ux;
                Uy(i,j) = uy;
                Uz(i,j) = uz;
                U(index,j) = ux;
                index = index+1;
                U(index,j) = uy;
                index = index+1;
                U(index,j) = uz;
                index = index+1;
            else
                GridID(index) = GID;
                GridType(index) = type(3);
                v = str2num(card(24:72));
                usp = v(1);
                if nargout == 12
                    Usp(isp,j) = usp;
                end
                U(index,j) = usp;
                index = index+1;
                isp = isp+1;
            end
            card = fgetl(fid);
            if type == '  G  '
                v = str2num(card(24:72));
                uxr = v(1);
                uyr = v(2);
                uzr = v(3);
                if nargout > 8
                    Uxr(i,j) = uxr;
                    Uyr(i,j) = uyr;
                    Uzr(i,j) = uzr;
                end
                U(index,j) = uxr;
                index = index+1;
                U(index,j) = uyr;
                index = index+1;
                U(index,j) = uzr;
                index = index+1;
                i=i+1;
            end
            card = fgetl(fid);
            if card ~= -1
                keepGoing = strcmp(card(1),' ') || contains(card,'-CONT-');
                if contains(card,'EIGENVALUE') && contains(card,'MODE')
                    justRead = 1;
                end
            else
                keepGoing = 0;
            end
        end
    end
end
fclose(fid);

ID = [];
k = 1;
for i = 1:length(GridID)
    if strcmp(GridType(i),'S')
        ID = [ID; GridID(i), 0];
    elseif strcmp(GridType(i),'G')
        ID = [ID; GridID(i), k];
        k = k+1;
        if k == 7
            k = 1;
        end
    end
end

if isp == 1
    Usp = [];
end

