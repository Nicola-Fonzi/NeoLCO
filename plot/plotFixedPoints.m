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

gapRatios = [0.5 2 8];
nGapRatios = length(gapRatios);

subspaces = repmat({[1,2,3]},nNonlinearities,1);
subspaces = combvec(subspaces{:});
nSubspaces = size(subspaces,2);

aeroMatrixUsed = modelForIntegration.aeroData.aeroMatrix_vlm.aeroMatrix;
aeroForce = aeroMatrixUsed.Kaero(aeroMatrixUsed.outputGroup.d,[aeroMatrixUsed.inputGroup.f0,aeroMatrixUsed.inputGroup.u])*[1; squeeze(resultsTrim.Vbody(:,resultsTrim.iLoadRec,1))];
aeroStiffness = aeroMatrixUsed.Kaero(aeroMatrixUsed.outputGroup.d,aeroMatrixUsed.inputGroup.d);

V_vect = linspace(0.001,max(options.speedVector),300);
Pdyn = V_vect.^2*options.rho*0.5;

model = modelForIntegration.model;
struData = modelForIntegration.struData;
nr = size(struData.suportTable,1);
ns = length(struData.surfNames);
nd = size(struData.Tgz,2) - nr - ns;
posd = 1:nd;

fixPoint = zeros(size(struData.Tgz,2),nSubspaces,nGapRatios,length(Pdyn));
fixpointAtNonlinearity = zeros(nNonlinearities,nSubspaces,nGapRatios,length(Pdyn));

DOF = obtainDOF(gapPoints,model);

for i = 1:nSubspaces
    
    offsetForce = zeros(size(struData.Tgz,1),1);
    offsetForce(DOF(subspaces(:,i)==1)) = stiffness(subspaces(:,i)==1);
    offsetForce(DOF(subspaces(:,i)==3)) = -stiffness(subspaces(:,i)==3);
    offsetForce = struData.Tgz'*offsetForce;
    
    tempStiffness = stiffness;
    tempStiffness(subspaces(:,i)==2)=0;
    [model_stiff, ~] = addNonlinearityStiffness(model, gapPoints, tempStiffness);
    struData_stiff = structuralPreprocessor(options.fidScreen, model_stiff, []);
    struStiffness = struData_stiff.Kzz(posd,posd);
    
    for j = 1:nGapRatios
        for k = 1:length(Pdyn)
            fixPoint(:,i,j,k) = (struStiffness-Pdyn(k)*aeroStiffness)\(Pdyn(k)*aeroForce*gapRatios(j) + offsetForce);
            fixpointAtNonlinearity(:,i,j,k) = struData.Tgz(DOF,:)*squeeze(fixPoint(:,i,j,k));
        end
    end
    
    clear tempStiffness
    
end

for y = 1:nNonlinearities
    figure
    hold on
    for i = 1:nSubspaces
        for j = 1:nGapRatios
            if j == 1
                plot(V_vect,squeeze(fixpointAtNonlinearity(y,i,j,:)),'bo')
            elseif j == 2
                plot(V_vect,squeeze(fixpointAtNonlinearity(y,i,j,:)),'ro')
            else
                plot(V_vect,squeeze(fixpointAtNonlinearity(y,i,j,:)),'go')
            end
        end
    end
    xlabel("Flow speed [m/s]")
    ylabel("Fixed points per unit gap size")
    title(strcat("Nonlinearity ",gapPoints{y,4}))
    h = gcf;
    saveas(h,strcat("Nonlinearity ",gapPoints{y,4}))
    close
end

return
