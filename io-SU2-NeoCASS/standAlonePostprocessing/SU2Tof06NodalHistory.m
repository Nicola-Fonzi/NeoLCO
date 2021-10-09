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
function SU2Tof06NodalHistory(filename_historymodal,filename_pch,filename_out,step_increment,analysed_modes)

% Read the time history of the modes and a punch file
% to write the time history of the nodal displacements
% in Nastran-like format in a .f06 file

if nargin < 4
    step_increment = 1;
end

if nargin < 3 || isempty(filename_out)
    filename_out = 'nodaldisp.f06';
else
    [filepath,name,ext] = fileparts(filename_out);
    if isempty(ext)
        filename_out = fullfile(filepath, strcat(name,'.f06'));
    elseif ~strcmp(ext,'.f06')
        error('filename_out extension must be ''.f06''');
    end
end


[ID,GridType,U,Ux,Uy,Uz,~,~,Uxr,Uyr,Uzr,Usp] = readPunchShapes(filename_pch);

node_list = unique(ID(:,1));
npoints = length(node_list);

i_G = unique(ID(GridType == 'G',1));
i_S = ID(GridType == 'S',1);

if nargin < 5
    analysed_modes = size(U,2);
end

for n = 1:analysed_modes
    [t,q] = readHistoryModal(filename_historymodal,analysed_modes,n,false);
    if n == 1
        jj = 1:step_increment:length(t);
        q_mat = zeros(n,length(jj));
    end
    q_mat(n,:) = q(jj)';
end
t = t(jj);

ux = Ux(:,1:analysed_modes)*q_mat;
uy = Uy(:,1:analysed_modes)*q_mat;
uz = Uz(:,1:analysed_modes)*q_mat;
uxr = Uxr(:,1:analysed_modes)*q_mat;
uyr = Uyr(:,1:analysed_modes)*q_mat;
uzr = Uzr(:,1:analysed_modes)*q_mat;
if ~isempty(Usp)
    usp = Usp(:,1:analysed_modes)*q_mat;
end


fido = fopen(filename_out,'w');
fprintf(fido,'1\n');
fprintf(fido,'\n');

for i = 1:npoints

    id = node_list(i);
    indexS = find(i_S==id,1);
    indexG = find(i_G==id,1);

    fprintf(fido,'1                                                       **STUDENT EDITION*      MAY  30, 2018  MSC Nastran  7/13/17   PAGE    1\n');
    fprintf(fido,'\n');
    fprintf(fido,'0\n');
    fprintf(fido,'      POINT-ID = %9d\n', id);
    fprintf(fido,'                                             D I S P L A C E M E N T   V E C T O R\n');
    fprintf(fido,'\n');
    fprintf(fido,'       TIME       TYPE          T1             T2             T3             R1             R2             R3\n');

    for j = 1:length(t)

        if indexG
            fprintf(fido,'%15.6e     G   %15.6e%15.6e%15.6e%15.6e%15.6e%15.6e\n',...
                t(j), ux(indexG,j), uy(indexG,j), uz(indexG,j), uxr(indexG,j), uyr(indexG,j), uzr(indexG,j));
        elseif indexS
            fprintf(fido,'%15.6e     S   %15.6e\n',...
                t(j), usp(indexS,j));
        end

    end

end

fclose(fido);


end
