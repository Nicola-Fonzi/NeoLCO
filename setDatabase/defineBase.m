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
function reducedBasis = defineBase(model, struData, rigidModes, globalOptions, options, nModes, kFactor, FMass)
% This function always assume that the model is provided to NeoLCO so that
% the stiffness at the nonlinearity points is NULL. It will be later added.
% Also, no other options should be specified for the nonlinearity points,
% like SUPORT or AESURFD cards.
%
% This function merely complete the computation of the obtimal base, based
% on the results of the obtainOptimalBase function

struOpt = [];
eigOpt = globalOptions.eig;
fid = options.fidScreen;

maximumStiffnesses = cellfun(@(x) max(x), options.kNominal);

if nargin < 6

    SetValues = true;

    while SetValues
        
        FMass = input("Chosen value for the fictitious mass: ");
        kFactor = input("Chosen value for the stiffness factor value: ");
        nModes = input("Chosen number of modes: ");
    
        SetValues = input(strcat("You requested fictitious mass = ",num2str(FMass),...
            "; stiffness factor = ",num2str(kFactor),"; number modes = ",num2str(nModes),"... Ok? [0=no, 1=yes]")) == 0;
        
    end

end

[model_common, ~] = addNonlinearityStiffness(model, options.gapPoints, maximumStiffnesses*kFactor);
[model_common, ~] = addFictitiousMass(model_common, options.gapPoints, repmat(FMass,size(options.gapPoints,1),1));

struData_common = structuralPreprocessor(fid, model_common, struOpt);

eigOpt.UseFictmass = true;
eigOpt.NROOTS = nModes;
resultsEig_common = solve_eig_fun(fid, model_common, struData_common, eigOpt);

reducedBasis = defineReducedBasis(struData, rigidModes, 'all', resultsEig_common, 'all');

if ~isempty(model.Tabdmp1.ID(model.Tabdmp1.ID==model.Tabdmp1.caseControl))
    KDAMP = model.Tabdmp1.KDAMP;
    if KDAMP < 0
        error("Only real damping is supported in this context, please use KDAMP>0")
    else
        fprintf(fid,' - Building real damping matrix for viscous damping...');
        reducedBasis.Bmm = modalDamp(model, reducedBasis.V'*struData.Mzz*reducedBasis.V, reducedBasis.V'*struData.Kzz*reducedBasis.V,'Real');
    end
    fprintf(fid, 'done.\n');
end

resultsFolder = strcat("FM",num2str(FMass),"KF",num2str(kFactor),"NMODES",num2str(nModes));

mkdir(resultsFolder)
chdir(resultsFolder)

save('structuralModel.mat', 'model', 'struData', 'reducedBasis')

close all

return