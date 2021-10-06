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
