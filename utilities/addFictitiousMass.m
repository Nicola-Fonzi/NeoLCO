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
function [model_loaded, CmassIndex] = addFictitiousMass(model, gapPoints, masses)

% This function takes as an input the free system (i.e. the one with no
% stiffness at the nonlinearity points) and add the stiffness contained in
% the vector "stiffnesses". The points where to add it are defined in the
% cell array "gapPoints" which has the format {ID1,"s" or "g", component1
% ;ID2,"s" or "g", component2 ...}. Other columns can be present, and
% create no issues.
%
% Please note that at the moment the code can only add a stiffness for a
% single dof. This implies that it can be either applied to a relative
% motion via a scalar point, or to an absolute motion via a single grid
% point. A relative motion between two grid points is not yet supported.

[~, pointIndex] = obtainDOF(gapPoints, model);

% The index of the new elements are large enough not to conflict with
% anything else in the model
CmassIndex = 4e9+(1:size(gapPoints,1));
cFictmassIndex = 5e9;

% Check that the components for the scalar points are 0
for i = 1:size(gapPoints,1)
    if gapPoints{i,2}=='s'
        gapPoints{i,3} = 0;
    end
end

model_loaded = model;

model_loaded.Fictmass.ID = [model_loaded.Fictmass.ID, cFictmassIndex];
model_loaded.Fictmass.elemList = cat(2, model_loaded.Fictmass.elemList, {CmassIndex});
for i = 1:size(gapPoints,1)
    model_loaded.Fictmass.Cmass.ID = [model_loaded.Fictmass.Cmass.ID, CmassIndex(i)];
    model_loaded.Fictmass.Cmass.point = cat(2,model_loaded.Fictmass.Cmass.point,[pointIndex(i);0]);
    model_loaded.Fictmass.Cmass.component = cat(2,model_loaded.Fictmass.Cmass.component,[gapPoints{i,3};0]);
    model_loaded.Fictmass.Cmass.value = [model_loaded.Fictmass.Cmass.value, masses(i)];
end

return
