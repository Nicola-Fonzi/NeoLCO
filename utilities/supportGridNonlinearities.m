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
% You are not entitled to use, distribute, or modify this file in any way,     *
% unless explicitly authorized by the copyright owner.                         *
%                                                                              *
%*******************************************************************************
function [gapPoints, model] = supportGridNonlinearities(gapPoints, model)

% This function is used to support the definition of nonlinearities between grid points
% In practise, we do what we would do in the linear FEM: we add a scalar point and a MPC between the two nodes

% Take the max ID to avoid overlap between the new scalar points and other nodes
maxs = max(model.Spoint.ID);
maxg = max(model.Node.ID);

maxID = max([maxs, maxg])+1;

% Extract the nonlinearities at grid points
labels = cellfun(@(x) x, gapPoints(:, 4));

duplicateLabel = [];
while ~isempty(labels)
    label = labels(1);
    I = find(strcmp(label,labels));
    if length(I)==1
        % This is a grid point fixed to the ground
        labels(1) = [];
        continue
    elseif length(I)==3
        % This is not possible
        error("Nonlinearities can be applied to scalar points, grid points, or between TWO (not more) grid points");
    end
    duplicateLabel = labels(1);
    labels(I) = [];
end

while ~isempty(duplicateLabel)
    % Find connected grid nodes in array
    indices =  find(cellfun(@(x) strcmp(x, duplicateLabel(1)), gapPoints(:,4)));

    % Add a new scalar point
    model.Spoint.ID = [model.Spoint.ID, maxID];

    % Add a new MPC connecting that scalar point
    model.MPC.SID = [model.MPC.SID, model.MPC.SID(1)];
    newMPC.G = [gapPoints{indices(1),1}, gapPoints{indices(2),1}, maxID];
    newMPC.DOF = obtainDOF({gapPoints{indices(1),:}; gapPoints{indices(2),:}; maxID, "s", 0, duplicateLabel(1)}, model).';
    newMPC.C = [gapPoints{indices(1),3}, gapPoints{indices(2),3}, 0];
    newMPC.A = [1, -1, -1];
    model.MPC.data = [model.MPC.data, newMPC];

    % Delete connected grid nodes from array
    gapPoints(indices,:) = [];

    % Add newly created scalar point
    gapPoints = {gapPoints{:,:}; maxID, "s", 0, duplicateLabel(1)};

    % Prepare for next iteration
    duplicateLabel(1) = [];
    maxID = maxID + 1;
end

end