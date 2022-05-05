clear
close all
clc

%% General inputs

% Input file
%filename_sma = 'D:\Desktop\PhD\X-DIA\ElevatorOnlyLCO\numericalSimulations\Final_DuringGVT_plate\correlated\steelWeightSteelConnection\NeoCASS\flutter.dat';
filename_sma = 'flutter.dat';

% Solver info
ver = get_neolco_version(true);

% Model specific values
%tipLE = 40002;
%tipTE = 40001;
%midspanLE = 40006;
%midspanElevator = 40015;
hingeScalarPoint = 1;
kNominal = 48;
airDensity = 1.13;

% Inputs for the definition of the common base
optimalBaseOptions.nModes = [8 12 16 20];
optimalBaseOptions.fictitiousMass = [1e-4 1e-2 1e0 1e2 1e4];
optimalBaseOptions.stiffnessFactor = [0.00001 0.001 0.1 0.5 1];
%optimalBaseOptions.optimumKnown = false;
optimalBaseOptions.optimumKnown = false;
%optimalBaseOptions.gapPoints = {hingeScalarPoint,'s',0,'Elevator hinge'};
optimalBaseOptions.gapPoints = {hingeScalarPoint,'s',0,'Rudder hinge'};
optimalBaseOptions.kNominal = {kNominal};

% Inputs for the aero database
aeroDatabaseOptions.selectionTrim = 'all';
aeroDatabaseOptions.DynVLM = false;
aeroDatabaseOptions.DynVLMtype = 'unsteady';

% Inputs for the time marching approach
preprocessTimeMarchingOptions.DynVLM = false;
preprocessTimeMarchingOptions.DynVLMtype = aeroDatabaseOptions.DynVLMtype;
preprocessTimeMarchingOptions.selectionTrim = 1;
preprocessTimeMarchingOptions.simulinkModel = 'SISO';
preprocessTimeMarchingOptions.gapPoints = optimalBaseOptions.gapPoints;
preprocessTimeMarchingOptions.gapBehaviour = {'freeplay','static',''};
preprocessTimeMarchingOptions.checkStateSpaceApproximation = false; %%%%%%%%%%%
%preprocessTimeMarchingOptions.monitorPoints = {tipLE,'g',3,'Accelerometer at the tip LE';...
 %                                              tipTE,'g',3,'Accelerometer at the tip TE';...
  %                                             midspanLE,'g',3,'Accelerometer at midpan LE';...
   %                                            midspanElevator,'g',3,'Accelerometer at midpan elevator'};                                           
timeMarchingOptions.gap = {[0.1, 0.2, 0.4, 0.6, 1, 1.5, 2, 2.5, 3]/180*pi};
timeMarchingOptions.kNominal = optimalBaseOptions.kNominal;
timeMarchingOptions.rho = airDensity;
timeMarchingOptions.speedVector = (14:2:50);
timeMarchingOptions.speedBehaviour = 'both'; 
timeMarchingOptions.speedInterpolationType = 0;
timeMarchingOptions.nFFTwindows = length(timeMarchingOptions.speedVector);
timeMarchingOptions.FFTwindowLength = 20;
timeMarchingOptions.overlapWindows = -10;
timeMarchingOptions.initialisationTime = 10;
%timeMarchingOptions.amplitudeDefinition = 'std';
timeMarchingOptions.amplitudeDefinition = 'maxPeak';
timeMarchingOptions.introduceFlightLoads = false;
timeMarchingOptions.trimType = 'grounded';
timeMarchingOptions.selectionTrim = 1;

% Inputs for the describing function approach
describingFunctionOptions.gapPoints = optimalBaseOptions.gapPoints;
describingFunctionOptions.gap = timeMarchingOptions.gap;
describingFunctionOptions.kNominal = timeMarchingOptions.kNominal;
describingFunctionOptions.amplitudeDefinition = timeMarchingOptions.amplitudeDefinition;
describingFunctionOptions.DynVLM = false;
describingFunctionOptions.DynVLMtype = aeroDatabaseOptions.DynVLMtype;
describingFunctionOptions.selectionTrim = 1;
describingFunctionOptions.recomputeBase = false;
describingFunctionOptions.searchQuenchPoint = true;
%describingFunctionOptions.maxKeq = kNominal;
describingFunctionOptions.maxKeq = kNominal/20;
describingFunctionOptions.nKeq = 25;
describingFunctionOptions.Vmax = 70;
describingFunctionOptions.Vmin = 5;
describingFunctionOptions.Vstep = 0.5;
describingFunctionOptions.method = 'PK0';
describingFunctionOptions.rho = airDensity;
describingFunctionOptions.modesPlot = 1:8;
describingFunctionOptions.axesUsed = 'body';

plotOptions.monitorNormalisation = 1;
plotOptions.torqueNormalisation = 0;

%% Loading of the model

inputData = readSmartcadFile(filename_sma);

%% Generate common basis

[reducedBasis, struData, model, struData_stiff, model_stiff, options] = obtainOptimalBase(inputData, optimalBaseOptions);

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
