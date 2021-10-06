% \file writef06_nodalhistory.m
%  \brief Append nondiagonal matrices to common punch format
%  \authors Nicola Fonzi, Vittorio Cavalieri
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

function [] = append_matrices_to_punch(filename,M,K,C)


fid = fopen(filename,'a');
writeMatrixElements(fid,K,'NDK');   % NON-DIAGONAL STIFFNESS MATRIX
writeMatrixElements(fid,M,'NDM');   % NON-DIAGONAL MASS MATRIX

if nargin == 4
    writeMatrixElements(fid,C,'NDC');   % NON-DIAGONAL DAMPING MATRIX
end

fclose(fid);

function writeMatrixElements(fid,A,keyword)
    fprintf(fid,'%s\n',keyword);
    n = size(A,1);
    for i = 1:n
        j = 1;
        fprintf(fid,'%-8d',i);
        while j <= n
            k = 1;
            while k <= 5 && j <= n
                fprintf(fid,'%14.6E',A(i,j));
                k = k+1;
                j = j+1;
            end
            fprintf(fid,'\n');
            if j <= n
                fprintf(fid,'%-8s','-CONT-');
            end
        end
    end
end

end
