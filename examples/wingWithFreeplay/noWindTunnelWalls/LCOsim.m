clear all
close all
clc

%% General inputs

% Input file
filename_sma = 'flutter.dat';

% Model specific values
hingeScalarPoint = 54000;
stiffFlutterSpeed = 23.35;
kNominal = 2.025;
airDensity = 2.425;

% Inputs for the definition of the common base
optimalBaseOptions.nModes = 3;
optimalBaseOptions.fictitiousMass = [1e-4 1e-2 1e0 1e2 1e4];
optimalBaseOptions.stiffnessFactor = [0.00001 0.001 0.1 0.5 1];
optimalBaseOptions.optimumKnown = false;
optimalBaseOptions.gapPoints = {hingeScalarPoint,"s",0,"Aileron hinge"};
optimalBaseOptions.kNominal = {kNominal};

% Inputs for the aero database
aeroDatabaseOptions.selectionTrim = 'all';
aeroDatabaseOptions.DynVLM = false;
aeroDatabaseOptions.DynVLMtype = 'unsteady';

% Input for the time marching
preprocessTimeMarchingOptions.DynVLM = false;
preprocessTimeMarchingOptions.DynVLMtype = aeroDatabaseOptions.DynVLMtype;
preprocessTimeMarchingOptions.selectionTrim = 1; % This selection is used for the dyn VLM
preprocessTimeMarchingOptions.simulinkModel = 'SISO';
preprocessTimeMarchingOptions.gapPoints = optimalBaseOptions.gapPoints;
preprocessTimeMarchingOptions.gapBehaviour = {'freeplay','static',''}; %{'freeplay','dynamic','increasing'}; % the simulink model has to be built around the definition of static or dynamic gap, it cannot be decided later on.
preprocessTimeMarchingOptions.monitorPoints = {100,'g',5,'Pitch';...
                                               100,'g',3,'Plunge'};
timeMarchingOptions.gap = {2.3/180*pi}; %{[1.15, 2.3]/180*pi};
timeMarchingOptions.kNominal = optimalBaseOptions.kNominal;
timeMarchingOptions.gapInterpolationType = 0;
timeMarchingOptions.rho = airDensity;
timeMarchingOptions.speedVector = (1:18);
timeMarchingOptions.speedBehaviour = 'both'; % The behaviour of speed can be easily changed without rebuilding the model
timeMarchingOptions.speedInterpolationType = 0;
timeMarchingOptions.nFFTwindows = length(timeMarchingOptions.speedVector);
timeMarchingOptions.FFTwindowLength = 20;
timeMarchingOptions.overlapWindows = -10;
timeMarchingOptions.initialisationTime = 10;
timeMarchingOptions.introduceFlightLoads = false;
timeMarchingOptions.trimType = 'grounded';
timeMarchingOptions.selectionTrim = 1; % This selection is used for the flight loads

% Input for the describing function
describingFunctionOptions.gapPoints = optimalBaseOptions.gapPoints;
describingFunctionOptions.gap = timeMarchingOptions.gap;
describingFunctionOptions.kNominal = timeMarchingOptions.kNominal;
describingFunctionOptions.DynVLM = false;
describingFunctionOptions.DynVLMtype = aeroDatabaseOptions.DynVLMtype;
describingFunctionOptions.selectionTrim = 1; % This selection is used for the dyn VLM
describingFunctionOptions.recomputeBase = false;
describingFunctionOptions.searchQuenchPoint = true;
describingFunctionOptions.maxKeq = kNominal;
describingFunctionOptions.nKeq = 25;
describingFunctionOptions.Vmax = 18;
describingFunctionOptions.Vmin = 5;
describingFunctionOptions.Vstep = 0.25;
describingFunctionOptions.method = 'PK0';
describingFunctionOptions.rho = airDensity;
describingFunctionOptions.modesPlot = 1:3;
describingFunctionOptions.axesUsed = 'body';
describingFunctionOptions.introduceFlightLoads = false;
describingFunctionOptions.trimType = 'grounded';
describingFunctionOptions.selectionTrim = 1; % This selection is used for the flight loads

plotOptions.halfGapNormalisation = 0;
plotOptions.normalisationSpeed = stiffFlutterSpeed;
plotOptions.monitorNormalisation = 1;
plotOptions.torqueNormalisation = 0;

%% Loading of the model

inputData = readSmartcadFile(filename_sma);

%% Generate common basis

[rigidModes, struData, model, struData_stiff, model_stiff, options, optimalBaseOptions] = ...
    obtainOptimalBase(inputData, optimalBaseOptions);

reducedBasis = defineBase(model, struData, rigidModes, options, optimalBaseOptions);

%% Generate aerodynamic database

[aeroData, aeroDatabaseOptions] = buildFullAeroData(model, model_stiff, struData_stiff, options, aeroDatabaseOptions);

%% Different equivalent stifnesses

[describingFunctionResults, describingFunctionOptions] = describingFunctions(model, struData, aeroData, options, reducedBasis, aeroDatabaseOptions, describingFunctionOptions);

%% Time marching simulation

[modelForIntegration, preprocessTimeMarchingOptions] = preprocessTimeMarchingLCO(model, struData, preprocessTimeMarchingOptions, reducedBasis, aeroData, options);

[timeMarchingResults, timeMarchingOptions]=timeMarchingLCO(modelForIntegration, timeMarchingOptions, preprocessTimeMarchingOptions);

%% Now, we plot

mainPlotRoutine(describingFunctionResults, describingFunctionOptions, ...
    timeMarchingResults, timeMarchingOptions, preprocessTimeMarchingOptions, plotOptions)
