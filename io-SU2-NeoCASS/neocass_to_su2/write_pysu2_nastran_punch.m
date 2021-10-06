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
function write_pysu2_nastran_punch(filename,K,U,ID_table)

% Read matrices and identifiers
% to write modal model (eigenvalues and mode shapes)
% in Nastran-like format in a .pch file

if isempty(filename)
    filename = 'model.pch';
else
    [filepath,name,ext] = fileparts(filename);
    if isempty(ext)
        filename = fullfile(filepath, strcat(name,'.pch'));
    elseif ~strcmp(ext,'.pch')
        error('filename extension must be ''.pch''');
    end
end

fid = fopen(filename,'w');

node_list = ID_table(:,1);
ngrid = sum(ID_table(:,2)==6);
nsp = sum(ID_table(:,2)==1);
GridType_list = strcat(repmat('G',1,ngrid), repmat('S',1,nsp));
npoints = length(node_list);
nmodes = size(U,2);

l = 1;

for j = 1:nmodes
    fprintf(fid,'$EIGENVALUE =  %9.7E  MODE = %5d%38d\n',K(j,j),j,l);
    l = l+1;
    index = 1;
    for i = 1:npoints
        id = node_list(i);
        if GridType_list(i) == 'G'
            fprintf(fid,'%10d       G     %13.6E     %13.6E     %13.6E%8d\n',id,...
                U(index,j),U(index+1,j),U(index+2,j),l);
            l = l+1;
            fprintf(fid,'-CONT-                 %13.6E     %13.6E     %13.6E%8d\n',...
                U(index+3,j),U(index+4,j),U(index+5,j),l);
        else
            fprintf(fid,'%10d       S     %13.6E     %13.6E     %13.6E%8d\n',id,...
                U(index,j),0,0,l);
            l = l+1;
            fprintf(fid,'-CONT-                 %13.6E     %13.6E     %13.6E%8d\n',...
                0,0,0,l);
        end
        l = l+1;
        if GridType_list(i) == 'G'
            index = index+6;
        else
            index = index+1;
        end
    end
end

fclose(fid);

end
