
function [model_stiff, celasIndex] = addNonlinearityStiffness(model, gapPoints, stiffnesses)

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
celasIndex = 3e9+(1:size(gapPoints,1));

% Check that the components for the scalar points are 0
gapPoints{cellfun(@(x) x=="s",gapPoints(:,2)),3}=0;

model_stiff = model;
for i = 1:size(gapPoints,1)
    model_stiff.Celas.ID = [model_stiff.Celas.ID, celasIndex(i)];
    model_stiff.Celas.point = cat(2,model_stiff.Celas.point,[pointIndex(i);0]);
    model_stiff.Celas.component = cat(2,model_stiff.Celas.component,[gapPoints{i,3};0]);
    model_stiff.Celas.value = [model_stiff.Celas.value, stiffnesses(i)];
end

return
