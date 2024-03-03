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
