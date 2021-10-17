function plotFixedPoints(modelForIntegration, stiffness, gapPoints, resultsTrim, options)

% This function only has a proper meaning if freeplay nonlinearities are
% concerned.
%
% In this case, each nonlinearity can be divided into three subspaces
%
%               force ^
%                     |     /
%                     |    / subspace 3
%       ______________|___/______________>
%                /                       displacement
%  subspace 1   /
%              /subspace 2
%
%
% The total number of subspaces depends on the number of nonlinearities
%
% Please note that this plot is made for unitary (in radians) gap HALF size

nNonlinearities = length(stiffness);

gapRatios = [0.5 2 6];
nGapRatios = length(gapRatios);

subspaces = repmat({[1,2,3]},nNonlinearities,1);
subspaces = combvec(subspaces{:});
nSubspaces = size(subspaces,2);

model = modelForIntegration.model;
struData = modelForIntegration.struData;
reducedBasis = modelForIntegration.reducedBasis;
aeroData = modelForIntegration.aeroData;

V_vect = linspace(0.001,max(options.speedVector),300);
Pdyn = V_vect.^2*options.rho*0.5;

nr = size(struData.suportTable,1);
ns = length(struData.surfNames);
nd = size(struData.Tgz,2) - nr - ns;

posd = 1:nd;

% We exclude rigid motion as it will never have a fixed point
struStiffness = struData.Kzz(posd,posd);

fixPoint = zeros(size(struStiffness,1),length(Pdyn),nSubspaces);
fixpointAtNonlinearity = zeros(nNonlinearities,length(Pdyn),nSubspaces,nGapRatios);

DOF = obtainDOF(gapPoints,model);

for i = 1:nSubspaces
    offsetForce = zeros(size(struData.Tgz,1),1);
    offsetForce(DOF(subspaces(:,i)==1)) = stiffness(subspaces(:,i)==1);
    offsetForce(DOF(subspaces(:,i)==3)) = -stiffness(subspaces(:,i)==3);
    offsetForce = struData.Tgz'*offsetForce;
    offsetStiffness = zeros(size(struData.Tgz,1),size(struData.Tgz,1));
    offsetStiffness(DOF(subspaces(:,i)==1 | subspaces(:,i)==3),DOF(subspaces(:,i)==1 | subspaces(:,i)==3)) = stiffness(subspaces(:,i)==1 | subspaces(:,i)==3);
    offsetStiffness = struData.Tgz'*offsetStiffness*struData.Tgz;
    for k = 1:nGapRatios
        for j = 1:length(Pdyn)
            fixPoint(:,j,i,k) = (struStiffness+offsetStiffness-Pdyn(j)*aeroStiffness)\(Pdyn(j)*aeroForce*gapRatios(k) + offsetForce);
            fixpointAtNonlinearity(:,j,i,k) = struData.Tgz(DOF,:)*squeeze(fixPoint(:,j,i));
        end
    end
end

for i = 1:nNonlinearities
    figure
    hold on
    for k = 1:nGapRatios
        if k == 1
            plot(V_vect,squeeze(fixpointAtNonlinearity(i,:,:,k)),'bo')
        elseif k == 2
            plot(V_vect,squeeze(fixpointAtNonlinearity(i,:,:,k)),'ro')
        else
            plot(V_vect,squeeze(fixpointAtNonlinearity(i,:,:,k)),'go')
        end
    end
    xlabel("Flow speed [m/s]")
    ylabel("Fixed points per unit gap size")
    title(strcat("Nonlinearity ",gapPoints{i,4}))
    h = gcf;
    saveas(h,strcat("Nonlinearity n°",num2str(i),".fig"))
    close
end

return
