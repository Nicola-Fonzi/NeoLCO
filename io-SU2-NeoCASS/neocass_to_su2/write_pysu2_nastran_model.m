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
function write_pysu2_nastran_model(model_data,filename,SID)

% Read the model_data from NeoCASS model
% to write model info (cord2r, grid, set1)
% in Nastran-like format in a .f06 file (ECHO)

if nargin < 3
    SID = 1;
end

if nargin < 2 || isempty(filename)
    filename = 'model.f06';
else
    [filepath,name,ext] = fileparts(filename);
    if isempty(ext)
        filename = fullfile(filepath, strcat(name,'.f06'));
    elseif ~strcmp(ext,'.f06')
        error('filename extension must be ''.f06''');
    end
end

if isstruct(model_data)
    model = model_data;
    maxID = 1e8-1;
    ID = model.Node.ID;
    index_grid = ID<=maxID;
    ID = ID(index_grid);
    coord = model.Node.Coord(:,index_grid);
    CD = model.Node.CD(index_grid);
    CID = model.Coord.ID;
    X0 = model.Coord.X0;
    R = model.Coord.R;
    set_list = model.Set.grids{model.Set.ID==SID};
    npoints = size(coord,2);
else
    coord = model_data;
    npoints = size(coord,2);
    ID = 1:npoints;
    CD = zeros(1,npoints);
    CID = [];
    set_list = ID;
end

if size(coord,1) ~= 3
    error('coord must have 3 rows');   
end



fido = fopen(filename,'w');
fprintf(fido,'1\n');
fprintf(fido,'\n');

fprintf(fido,'1                                                       **STUDENT EDITION*      MAY  30, 2018  MSC Nastran  7/13/17   PAGE    1\n');
fprintf(fido,'\n');
fprintf(fido,'0\n');
fprintf(fido,'                                                  S O R T E D   B U L K   D A T A   E C H O                                         \n');
fprintf(fido,'                 ENTRY                                                                                                              \n'); 
fprintf(fido,'                 COUNT        .   1  ..   2  ..   3  ..   4  ..   5  ..   6  ..   7  ..   8  ..   9  ..  10  . \n'); 

k = 1;
RID = 0;
for j = 1:length(CID)
    A = X0(:,j);
    sA = [num2Nastranfield(A(1)),num2Nastranfield(A(2)),num2Nastranfield(A(3))];
    B = A + R(:,:,j) * [0;0;1];
    sB = [num2Nastranfield(B(1)),num2Nastranfield(B(2)),num2Nastranfield(B(3))];
    C = A + R(:,:,j) * [1;0;0];
    sC = [num2Nastranfield(C(1)),num2Nastranfield(C(2)),num2Nastranfield(C(3))];
    fprintf(fido,'%21d-        CORD2R  %-8d%-8d%24s%24s+\n',k,CID(j),RID,sA,sB);
    k = k+1;
    fprintf(fido,'%21d-        +       %24s\n',k,sC);
    k = k+1;
end

CP = 0;
for j = 1:npoints
    if CD(j) == 0
        cid = 0;
    else
        cid = CID(CD(j));
    end
    coord_j = [num2Nastranfield(coord(1,j)),num2Nastranfield(coord(2,j)),num2Nastranfield(coord(3,j))];
    fprintf(fido,'%21d-        GRID    %-8d%-8d%24s%d\n',k,ID(j),CP,coord_j,cid);
    k = k+1;
end

n = length(set_list);

if n <= 7
    fmt = ['%21d-        SET1    %-8d',repmat('%-8d',1,n),'\n'];
    fprintf(fido,fmt,k,SID,set_list);
else
    rows = floor(n/8)+1;
    fprintf(fido,'%21d-        SET1    %-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d+\n',k,SID,set_list(1:7));
    k = k+1;
    for i = 2:rows-1
        ids = set_list(8*(i-1):8*i-1);
        fprintf(fido,'%21d-        +       %-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d+\n',k,ids);
        k = k+1;
    end
    if isempty(i)
        i = 2;
    else
        i = i+1;
    end
    ids = set_list(8*(i-1):n);
    fmt = ['%21d-        +       ',repmat('%-8d',1,length(ids)),'\n'];
    fprintf(fido,fmt,k,ids);
end

fclose(fido);



function s = num2Nastranfield(a)
    
    a_max = 1e10-1;
    if abs(a) > a_max
        error('%d is the maximum coordinate value supported for conversion',a_max);
    end
    
    a_min = 1e-9;
    if a ~= 0 && abs(a) < a_min
        a = 0;
        warning('Small coordinate value (less than %.1e) is rounded to 0.',a_min);
    end
    
    if a == 0
        s = '0.      ';
    else
        order = floor(log10(abs(a)));
        if order <= -3
            temp = num2str(a,'%.3e');
            index = strfind(temp,'e-0');
            temp([index,index+2]) = '';
        elseif order == -2
            temp = num2str(a,'%.5g');
        elseif order >= -1 && order <= 4
            temp = num2str(a,'%.6g');
        elseif order == 5
            temp = num2str(a,'%.6g');
        else
            temp = num2str(a,'%.3e');
            index = strfind(temp,'e+0');
            temp([index,index+2]) = '';
        end
        if ~contains(temp,'.')
            temp = [temp, '.'];
        end
        if strcmp(temp(1:2),'0.')
            temp(1) = '';
        elseif contains(temp,'-0.')
            temp(2) = '';
        end
        l = length(temp);
        if l < 8
            s = [temp,repmat(' ',1,8-l)];
        elseif l == 8
            s = temp;
        else
            error('Maximum of 8 characters for field is exceeded');
        end
    end
    
end

end