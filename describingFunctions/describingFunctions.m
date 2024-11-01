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
function [results, options] = describingFunctions(model, struData, aeroData, globalOptions, reducedBasis, aeroDatabaseOptions, options)

% Set default options
iOpt = 0;
iOpt = iOpt+1; baseOpt.fidScreen = 1;                 descr{iOpt} = 'Fid for screen printing. [1].';
iOpt = iOpt+1; baseOpt.useInApp = 0;                  descr{iOpt} = 'Utility flag to be used when called from App';
% General options
iOpt = iOpt+1; baseOpt.DynVLM = false;                descr{iOpt} = 'Request the use of the dynamic VLM. [false].';
iOpt = iOpt+1; baseOpt.selectionTrimDynVLM = [];      descr{iOpt} = 'Selected trim case to be used when correcting the matrices with VLM results. [].';
iOpt = iOpt+1; baseOpt.DynVLMtype = 'unsteady';       descr{iOpt} = 'Type of dyn VLM. ["unsteady"].';
iOpt = iOpt+1; baseOpt.gapPoints = {};                descr{iOpt} = 'Points where nonlinearities are present. The format is {point1,"s" or "g",dof (1,2,3,4,5 or 6),label;point2,...}. If in the same point we have more nonlinearity, the point must be repeated per each dof. {}.';
iOpt = iOpt+1; baseOpt.gap = {};                      descr{iOpt} = 'Cells array containing the peak-to-peak possible gaps at the nonlinearity points. The cells are in the same order as gapPoint IDs. The columns contain possible different values for the same point. The format is the same that we use for the kNominal. {}.';
iOpt = iOpt+1; baseOpt.kNominal = {};                 descr{iOpt} = 'Cells array containing the nominal stiffnesses at the nonlinearity points. The rows are in the same order as gapPoint IDs. The columns contain possible different values for the same point. The format is {[k1_point1,k2_point1];[k1_point2,k2_point2,k3_point2]}. {}.';
iOpt = iOpt+1; baseOpt.rho = 1.225;                   descr{iOpt} = 'Density used for the flutter analysis. [1.225].';
iOpt = iOpt+1; baseOpt.machNumber = [];               descr{iOpt} = 'Mach number to be used for the state-space approximation of the aerodynamics. [].';
iOpt = iOpt+1; baseOpt.amplitudeDefinition = 'rms';   descr{iOpt} = 'Algorithm to compute the LCO amplitude. "std", "rms" or "maxPeak". ["rms"].';
% Options used in case of PK flutter analysis
iOpt = iOpt+1; baseOpt.recomputeBase = false;         descr{iOpt} = 'Request the recomputation of the base per each equivalent sstiffness. [false]';
iOpt = iOpt+1; baseOpt.searchQuenchPoint = true;      descr{iOpt} = 'Request to dynamically adapt the equivalent stiffness step to search for the quanching point. [true].';
iOpt = iOpt+1; baseOpt.maxKeq = [];                   descr{iOpt} = 'Maximum equivalent stiffness to be considered. [].';
iOpt = iOpt+1; baseOpt.nKeq = [];                     descr{iOpt} = 'Total number of equaivalent stiffnesses to be considered. [].';
iOpt = iOpt+1; baseOpt.Vmax = 100;                    descr{iOpt} = 'Maximum speed for the flutter analysis. [100].';
iOpt = iOpt+1; baseOpt.Vmin = 10;                     descr{iOpt} = 'Minimum speed for the flutter analysis. [10].';
iOpt = iOpt+1; baseOpt.Vstep = 5;                     descr{iOpt} = 'Delta speed for the flutter analysis. [5].';
iOpt = iOpt+1; baseOpt.method = 'PK0';                descr{iOpt} = 'Method used for the flutter analysis. [PK0].';
iOpt = iOpt+1; baseOpt.modesPlot = 1:3;               descr{iOpt} = 'Modes to be plotted in the Vg diagrams. [1:3].';
iOpt = iOpt+1; baseOpt.axesUsed = 'body';             descr{iOpt} = 'Axes used in the flutter analysis. ["body"].';
% Options used in case of eigenvalue analysis
% For the introduction of steady loads
iOpt = iOpt+1; baseOpt.introduceFlightLoads = false;  descr{iOpt} = 'Flag to introduce flight steady loads in the analysis. [false].';
iOpt = iOpt+1; baseOpt.trimType = 'meanAxes';         descr{iOpt} = 'Type of trim output to be requested. ["meanAxes"],';
iOpt = iOpt+1; baseOpt.selectionTrim = 1;             descr{iOpt} = 'Selected trim ID to be used for the steady load calculation. [1].';
iOpt = iOpt+1; baseOpt.introduceStruLoads = false;    descr{iOpt} = 'Flag to introduce constant mechanical preloads';
iOpt = iOpt+1; baseOpt.struLoads = {};                descr{iOpt} = 'Option to specify the locations of the loads. The format is {point1,"s" or "g",dof (1,2,3,4,5 or 6),load;point2,...}. {}.';
iOpt = iOpt+1; baseOpt.introduceGravityLoads = false; descr{iOpt} = 'Flag to introduce gravity loads in the analysis. [false].';
iOpt = iOpt+1; baseOpt.gravDirection = [0,0,-1];      descr{iOpt} = 'Direction of the gravity loads, when applied. [0,0,-1].';

if nargin==0
    printOptionDescription(baseOpt, descr);
    return
end

% Process input options
options = setOptions(baseOpt, 'error', options);

if isempty(options.machNumber)
    % If not specified we use the first Mach number
    options.machUsed = 1;
else
    options.machUsed = aeroData.dlmData.aero.M==options.machNumber;
end

% Derived options
options.Mlist_dlm = globalOptions.aero.Mlist_dlm;
options.klist = globalOptions.aero.klist;
options.eigOpt = globalOptions.eig;
options.struOpt = [];

mkdir("DescribingFunctions")
chdir("DescribingFunctions")

if isempty(options.gapPoints)
    error("You probably forgot to set the gapPoints option, identifying where nonlinearities are present")
end

if size(options.gapPoints,1)==1
    results = describingFunctionPK(model, struData, aeroData, reducedBasis, globalOptions, aeroDatabaseOptions, options);
else
    error("Multiple nonlinearities not yet implemented");
end

if length(options.selectionTrim)>1 || length(options.selectionTrimDynVLM)>1
    error(strcat("Only one trim condition can be specified to be used for the correction via VLM of the DLM "...
        ,"matrices or for the introduction of flight loads"))
end

cd("..")


return