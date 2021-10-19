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
function [DOF, pointIndex] = obtainDOF(gridCellArray,model)
%
% This function obtain the global degree of freedom of a certain
% displacement component of a node in the model. It uses the ID of the node
% (or scalar point) and the component number.
%
% These inputs must be provided via a cell array with the format {ID1,"s"
% or "g",component1;ID2 ...}. Where "s" is used for scalar points, "g" for
% physical grid nodes.

[~, gridDOF, spointDOF] = getIDtable(model.Node.ID, model.Spoint.ID);

DOF = zeros(size(gridCellArray,1),1);
if nargout>1
    pointIndex= zeros(size(gridCellArray,1),1);
end

for i = 1:length(DOF)
    if strcmp(gridCellArray{i,2},'s')
        Pos = find(model.Spoint.ID == gridCellArray{i,1},1);
        DOF(i) = spointDOF(Pos);
    else
        Pos = find(model.Node.ID == gridCellArray{i,1},1);
        DOF(i) = gridDOF(Pos,gridCellArray{i,3});
    end
    if nargout>1
        pointIndex=Pos;
    end
end

return
