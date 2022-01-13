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
function [modelForIntegration, options] = preprocessTimeMarchingLCO(model, struData, options, reducedBasis, aeroData, globalOptions)
%
% This function assumes that a modal base has already
% been defined, and the aerodynamic database has also been built with
% previous flutter and trim analysis
%
% The points where a nonlinearity must be added, should be in the free
% condition. That is, no stiffness should be associated with them. Also, if
% a scalar point is chosen as a nonlinearity point, there should not be
% associated an AESURFD. That is, a simple rigid motion must be allowed to
% the grid or scalar point where nonlinearity is present.
%
% The function will create the state-space model for the aerodynamic, it
% will create the simulink model for integration, and it will double check
% the aerodynamic approximation.
%
% The resulting model used for integration misses only some parameters,
% like the nominal stiffness at the nonlinearity points, and the other
% parameters defining the nonlinearity. Plus, the constant forcing terms,
% like the flight loads and constant loads. This is because, once the model
% is prepared, it can be reused changing the above values without
% rebuilding it. This implies that monitor points and nonlinearity points
% must be defined here and not changed later. This also include the type of
% nonlinearity and its behaviour. It is only excluded the behaviour of the
% wind speed (i.e. if it is constant, increasing, decreasing, ecc...) as
% that can be defined in the time marching directly.

% Set default options
iOpt = 0;
iOpt = iOpt+1; baseOpt.fidScreen = 1;                 descr{iOpt} = 'fid for screen printing. [1].';
% Options for the state-space realisation
iOpt = iOpt+1; baseOpt.DynVLM = false;                descr{iOpt} = 'Request the use of the VLM corrected matrices for unsteady analysis. [false].';
iOpt = iOpt+1; baseOpt.DynVLMtype = 'unsteady';       descr{iOpt} = 'Type of unsteady VLM correction to be applied';
iOpt = iOpt+1; baseOpt.selectionTrim = [];            descr{iOpt} = 'Selected trim case to be used when correcting the matrices with VLM results. [].';
iOpt = iOpt+1; baseOpt.machNumber = [];               descr{iOpt} = 'Mach number to be used for the state-space approximation of the aerodynamics. [].';
iOpt = iOpt+1; baseOpt.eigsopt.threshold = 0;         descr{iOpt} = 'Options for improvedMFDfun. See that function for details. [0].';
iOpt = iOpt+1; baseOpt.eigsopt.bound = -1e-4;         descr{iOpt} = 'Options for improvedMFDfun. See that function for details. [-1e-4].';
iOpt = iOpt+1; baseOpt.eigsopt.type = 'bound';        descr{iOpt} = 'Options for improvedMFDfun. See that function for details. ["bound"].';
iOpt = iOpt+1; baseOpt.eigsopt.method = 'eigshift';   descr{iOpt} = 'Options for improvedMFDfun. See that function for details. ["eigshift"].';
iOpt = iOpt+1; baseOpt.algROM = 'balance';            descr{iOpt} = 'Options for improvedMFDfun. See that function for details. ["balance"].';
iOpt = iOpt+1; baseOpt.optsLM = [0.1, 1.0e-6, 1.0e-6, 100];         descr{iOpt} = 'Options for improvedMFDfun. [tau tolg tolx maxIter] See that function for details. [0.1, 1.0e-6, 1.0e-6, 100].';
iOpt = iOpt+1; baseOpt.opt = {1, 1, 'lmfd', 2, 100};  descr{iOpt} = 'Options for improvedMFDfun. {mfd order, mfd algorithm, r3, numerator order, weight} See that function for details. [{1, 1, "lmfd", 2, 100}].';
% Options to produce graphical outputs
iOpt = iOpt+1; baseOpt.checkStateSpaceApproximation = false;        descr{iOpt} = 'This flag generates several plots that can be used to check the aerodynamic approximation. [false].';
% Options to build the time marching model
iOpt = iOpt+1; baseOpt.simulinkModel = 'model';       descr{iOpt} = 'Name of the simulink model to be generated. [model.slx]';
% Localise points of interest
iOpt = iOpt+1; baseOpt.gapPoints = {};                descr{iOpt} = 'Points where nonlinearities are present. The format is {point1,"s" or "g",dof (1,2,3,4,5 or 6),label;point2,...}. If in the same point we have more nonlinearity, the point must be repeated per each dof. {}.';
iOpt = iOpt+1; baseOpt.monitorPoints = {};            descr{iOpt} = 'Points to be used as monitor. For these points, the time history will be saved. The format is the same as above. {}.';
iOpt = iOpt+1; baseOpt.gapBehaviour = {};             descr{iOpt} = 'Each row will control the corresponding gap in gapPoints. The format is {"freeplay","static","";"freeplay","dynamic","increasing";"freeplay","dynamic","both"}. The previous instruction says: the first gap, with freeplay nonlinearity, does not change (static), the second changes (dynamic) increasing in size, the third increasing and then decreasing. Values provided in gap are used as interpolation values in case of "d". {}';

if nargin==0
    printOptionDescription(baseOpt, descr);
    return
end

% Process input options
options = setOptions(baseOpt, 'error', options);

if isempty(options.gapPoints) || isempty(options.gapBehaviour)
    error("You probably forgot to set the gapPoints option, identifying where nonlinearities are present, or their behaviour")
end

if length(options.selectionTrim)>1
    error("Only one trim condition can be specified to be used for the correction via VLM of the DLM matrices")
end

if any(cellfun(@(x) strcmp(x,"dynamic"),options.gapBehaviour(:,2))) && any(cellfun(@(x) strcmp(x,"static"),options.gapBehaviour(:,2)))
    error("A combination of static and dynamic gap must be set by requesting all dynamic gaps, but setting only one gap size for the ones we want to remain fixed")
end

nonlinearityBlock = strcat(options.simulinkModel,"/nonlinearity/");
nNonlinearities = size(options.gapPoints,1);

mkdir("TimeMarching")
chdir("TimeMarching")

%% Generate the DLM matrices if not present already

if ~isfield(aeroData, 'aeroMatrix_dlm')
    aeroOptions.gustShape = [];
    aeroData.aeroMatrix_dlm = getDLMmatrices(options.fidScreen, aeroData.dlmData, model, aeroData.lattice_dlm, ...
	                                aeroData.interpData_flat, ...
	                                aeroData.interpMat_dlm.IcList, struData, ...
	                                reducedBasis, aeroOptions);
end

%% State-space approximation of aerodynamics

[solution, chosenModes, machUsed, k, Ha, selectedTrim] = convertAeroInSS(aeroData, model, struData, reducedBasis, options);

%% Create the system for the simulink model

% Structural matrices
modelForIntegration.Mhh = reducedBasis.V(:,chosenModes)'*struData.Mzz*reducedBasis.V(:,chosenModes);
modelForIntegration.Khh = reducedBasis.V(:,chosenModes)'*struData.Kzz*reducedBasis.V(:,chosenModes);
modelForIntegration.Chh = reducedBasis.Bmm(chosenModes,chosenModes);
modelForIntegration.invMC = modelForIntegration.Mhh\modelForIntegration.Chh;
modelForIntegration.invMK = modelForIntegration.Mhh\modelForIntegration.Khh;
modelForIntegration.nStru = size(modelForIntegration.invMC,1);

% Obtain the mapping from modes to nonlinearity points
nonlinearityDOF = obtainDOF(options.gapPoints,model);
modelForIntegration.Ugap = struData.Tgz(nonlinearityDOF,:)*reducedBasis.V(:,chosenModes);


% Obtain mapping from modes to monitor points
if ~isempty(options.monitorPoints)
    fprintf(options.fidScreen, "WARNING: As of today, the monitor points are output in the global reference system\n")
    monitorDOF = obtainDOF(options.monitorPoints,model);
    modelForIntegration.Umonitor = struData.Tgz(monitorDOF,:)*reducedBasis.V(:,chosenModes);
end

% Aerodynamic matrices
modelForIntegration.D0 = solution.inoutresid.DD{1};
modelForIntegration.D1 = solution.inoutresid.DD{2};
modelForIntegration.D2 = solution.inoutresid.DD{3};
modelForIntegration.B0 = solution.inoutresid.BB{1};
modelForIntegration.B1 = solution.inoutresid.BB{2};
modelForIntegration.B2 = solution.inoutresid.BB{3};
modelForIntegration.Ca = solution.inoutresid.CC;
modelForIntegration.A = solution.inoutresid.AA;
modelForIntegration.nAero = size(solution.inoutresid.AA,1);
% lref must be the one used for the state-space representiation (i.e. the DLM one)
modelForIntegration.lref = aeroData.aeroMatrix_dlm.aero.lref;

modelForIntegration.reducedBasis.V = reducedBasis.V(:,chosenModes);
modelForIntegration.aeroData = aeroData;
modelForIntegration.machUsed = machUsed;
modelForIntegration.selectedTrim = selectedTrim;
modelForIntegration.struData = struData;
modelForIntegration.model = model;
modelForIntegration.globalOptions = globalOptions;

%% Now create the simulink model

templatePath = which("templateLCO.slx");
copyfile(templatePath,strcat(options.simulinkModel,".slx"))
load_system(options.simulinkModel)

centerPosition = [725 545];
dimensions = [30 30];

if any(cellfun(@(x) strcmp(x,"dynamic"),options.gapBehaviour(:,2)))
    dynamicNonlinearities = true;
else
    dynamicNonlinearities = false;
end

set_param(strcat(nonlinearityBlock,"divideNonlinearities"),"Outputs",num2str(nNonlinearities));
set_param(strcat(nonlinearityBlock,"mergeNonlinearities"),"Inputs",num2str(nNonlinearities));

for iGap = 1:nNonlinearities
    if dynamicNonlinearities
        switch(options.gapBehaviour{iGap,1})
            case "freeplay"
                add_block("dynamicFreeplay/dynamicGap",strcat(options.simulinkModel,"/nonlinearity/nonlinearity",num2str(iGap)),"position",[centerPosition(1), centerPosition(2)+2*dimensions(2)*(iGap-1), centerPosition(1)+dimensions(1), centerPosition(2)+dimensions(2)*(2*iGap-1)]);
                set_param(strcat(nonlinearityBlock,"nonlinearity",num2str(iGap),"/stiffness"),"Gain",strcat("-stiffness(",num2str(iGap),")"));
                set_param(strcat(nonlinearityBlock,"nonlinearity",num2str(iGap),"/freeplayVector"),"Value",strcat("freeplayVector{",num2str(iGap),"}"));
                set_param(strcat(nonlinearityBlock,"nonlinearity",num2str(iGap),"/interpolationType"),"Value",strcat("options.gapInterpolationType(",num2str(iGap),")"));
            otherwise
                error(strcat(options.gapBehaviour{iGap,1}," nonlinearity not yet supported"))
        end
        add_line(nonlinearityBlock,strcat("divideNonlinearities/",num2str(iGap)),strcat("nonlinearity",num2str(iGap),"/1"),"autorouting","smart");
        add_line(nonlinearityBlock,strcat("nonlinearity",num2str(iGap),"/1"),strcat("mergeNonlinearities/",num2str(iGap)),"autorouting","smart");
    else
        switch(options.gapBehaviour{iGap,1})
            case "freeplay"
                add_block("staticFreeplay/staticGap",strcat(options.simulinkModel,"/nonlinearity/nonlinearity",num2str(iGap)),"position",[centerPosition(1), centerPosition(2)+2*dimensions(2)*(iGap-1), centerPosition(1)+dimensions(1), centerPosition(2)+dimensions(2)*(2*iGap-1)]);
                set_param(strcat(nonlinearityBlock,"nonlinearity",num2str(iGap),"/gap"),"LowerValue",strcat("-freeplay(",num2str(iGap),")"));
                set_param(strcat(nonlinearityBlock,"nonlinearity",num2str(iGap),"/gap"),"UpperValue",strcat("freeplay(",num2str(iGap),")"));
                set_param(strcat(nonlinearityBlock,"nonlinearity",num2str(iGap),"/stiffness"),"Gain",strcat("-stiffness(",num2str(iGap),")"));
            otherwise
                error(strcat(options.gapBehaviour{iGap,1}," nonlinearity not yet supported"))
        end
        add_line(nonlinearityBlock,strcat("divideNonlinearities/",num2str(iGap)),strcat("nonlinearity",num2str(iGap),"/1"),"autorouting","smart");
        add_line(nonlinearityBlock,strcat("nonlinearity",num2str(iGap),"/1"),strcat("mergeNonlinearities/",num2str(iGap)),"autorouting","smart");
    end
end

save_system(options.simulinkModel)
close_system(options.simulinkModel)

%% Plot the approximation to double check

if options.checkStateSpaceApproximation

    compareTransferFunctionWithAeroSS(modelForIntegration,Ha,k)

end

return
