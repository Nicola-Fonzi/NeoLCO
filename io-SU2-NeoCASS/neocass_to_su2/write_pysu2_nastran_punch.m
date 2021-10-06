% \file write_pysu2_nastran_punch.m
%  \brief Writes modal model (eigenvalues and mode shapes) in Nastran-like format (.pch)
%  \authors Vittorio Cavalieri, Nicola Fonzi
%  \version 7.0.8 "Blackbird"
%
% SU2 Project Website: https://su2code.github.io
%
% The SU2 Project is maintained by the SU2 Foundation
% (http://su2foundation.org)
%
% Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
%
% SU2 is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% SU2 is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with SU2. If not, see <http://www.gnu.org/licenses/>.

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
