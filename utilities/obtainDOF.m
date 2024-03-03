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
        pointIndex(i)=Pos;
    end
end

return
