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
function [results, options] = timeMarchingLCO(modelForIntegration, options, preprocessTimeMarchingOptions)
%
% Important note: the computation of trim condition for the loads is velocity dependent,
% thus it should be recomputed per each flight condition. In practise, it
% is assumed that the trim condition specified is close enough to the
% conditions at which the system will operate and only one calculation is
% applied.

% Set default options
iOpt = 0;
% For screen output
iOpt = iOpt+1; baseOpt.fidScreen = 1;                 descr{iOpt} = 'fid for screen printing. [1]';
% For parameters of the problem
iOpt = iOpt+1; baseOpt.rho = 1.225;                   descr{iOpt} = 'Air density to be used for the different simulations. [1.225].';
iOpt = iOpt+1; baseOpt.speedVector = (10:2:56);       descr{iOpt} = 'Vector of speeds to analyse. [(10:2:56) m/s].';
iOpt = iOpt+1; baseOpt.speedBehaviour = 'both';       descr{iOpt} = 'Specify wether the speed will be increased "increasing", decreased "decreasing", or both "both". ["b"]';
iOpt = iOpt+1; baseOpt.speedInterpolationType = 0;    descr{iOpt} = 'Interpolation type to be used when changing the speed dynamically. Availables 0 for "previous", 1 for "linear", or 2 for "spline". [0]';
iOpt = iOpt+1; baseOpt.kNominal = {};                 descr{iOpt} = 'Cells array containing the nominal stiffnesses at the nonlinearity points. The rows are in the same order as gapPoint IDs. The columns contain possible different values for the same point. The format is {[k1_point1,k2_point1];[k1_point2,k2_point2,k3_point2]}. {}.';
iOpt = iOpt+1; baseOpt.gap = {};                      descr{iOpt} = 'Cells array containing the peak-to-peak possible gaps at the nonlinearity points. The cells are in the same order as gapPoint IDs. The columns contain possible different values for the same point. The format is the same that we use for the kNominal. {}.';
iOpt = iOpt+1; baseOpt.gapInterpolationType = [];     descr{iOpt} = 'Integer array containing the interpolation type to be used when changing the gap dynamically. One row per each gap point. Availables 0 for "previous", 1 for "linear", or 2 for "spline". [].';
iOpt = iOpt+1; baseOpt.initialAmplitude = 2;          descr{iOpt} = 'Factor that multiplies the freeplay to obtain the initial condition. If set to 2 it means that the initial amplitude will be twice the freeplay size. [2].';
% For output computation
iOpt = iOpt+1; baseOpt.initialisationTime = 10;       descr{iOpt} = 'Time over which the initialisation is assumed completed, the speed is kept constant equal to the first element of speedVector. [10 seconds].';
iOpt = iOpt+1; baseOpt.nFFTwindows = 1;               descr{iOpt} = 'Number of windows in which to compute the FFT. [1]';
iOpt = iOpt+1; baseOpt.FFTwindowLength = 20;          descr{iOpt} = 'Lenght in seconds of each FFT window';
iOpt = iOpt+1; baseOpt.timeStepForFFT = 1e-3;         descr{iOpt} = 'Simulink may use a variable time step, so we reconstruct the signal using a fixed time step with size specified in this variable. [1e-3].';
iOpt = iOpt+1; baseOpt.overlapWindows = 0;            descr{iOpt} = 'Time in seconds of overlapping of the FFT windows. Can be used negative to create a spacing. [0 s]';
iOpt = iOpt+1; baseOpt.maxFrequencyInterest = 150;    descr{iOpt} = 'Maximum frequency to retain in the FFT';
iOpt = iOpt+1; baseOpt.amplitudeDefinition = 'rms';   descr{iOpt} = 'Algorithm to compute the LCO amplitude. "std", "rms" or "maxPeak". ["rms"].';
iOpt = iOpt+1; baseOpt.plotMaxRealEigSS = true;       descr{iOpt} = 'Use this option to plot the linear flutter speed of the free and stiff system, as obtained via the eigenvalues of the problem ["true"].';
iOpt = iOpt+1; baseOpt.plotFixedPoints = true;        descr{iOpt} = 'Use this option to plot the fixed points of the nonlinear system. In a system with freeplay, for example, there can be three fixed points. ["true"].';
% For the introduction of steady loads
iOpt = iOpt+1; baseOpt.introduceFlightLoads = false;  descr{iOpt} = 'Flag to introduce flight steady loads in the analysis. [false].';
iOpt = iOpt+1; baseOpt.trimType = 'meanAxes';         descr{iOpt} = 'Type of trim output to be requested. ["meanAxes"],';
iOpt = iOpt+1; baseOpt.selectionTrim = 1;             descr{iOpt} = 'Selected trim ID to be used for the steady load calculation. [1].';
iOpt = iOpt+1; baseOpt.introduceStruLoads = false;    descr{iOpt} = 'Flag to introduce constant mechanical preloads';
iOpt = iOpt+1; baseOpt.struLoads = {};                descr{iOpt} = 'Option to specify the locations of the loads. The format is {point1,"s" or "g",dof (1,2,3,4,5 or 6),load;point2,...}. {}.';

if nargin==0
    printOptionDescription(baseOpt, descr);
    return
end

% Process input options
options = setOptions(baseOpt, 'error', options);

% Derived options
options.gapPoints = preprocessTimeMarchingOptions.gapPoints;
options.gapBehaviour = preprocessTimeMarchingOptions.gapBehaviour;
options.monitorPoints = preprocessTimeMarchingOptions.monitorPoints;
options.simulinkModel = preprocessTimeMarchingOptions.simulinkModel;

% Check if the required information is present
if size(options.kNominal,1)~=size(options.gapPoints,1)
    error("Nominal stifnesses are compulsory per each point")
elseif size(options.gap,1)~=size(options.gapPoints,1)
    error("Gap dimensions are compulsory per each point")
end

dynamicNonlinearities = any(cellfun(@(x) strcmp(x,"dynamic"),options.gapBehaviour(:,2)));
if dynamicNonlinearities
    if length(options.gapInterpolationType)~=size(options.gapPoints,1)
        error("Gap interpolation type must be defined per each nonlinearity if they are dynamic")
    end
end

% In case VLM correction of DLM matrices has been used, check that the same
% air density has been used. We already checked that the same Mach number
% has been used when merging the two matrices.
if ~isempty(modelForIntegration.selectedTrim)
    if options.rho~=modelForIntegration.aeroData.vlmData.DynVLM.trimRes.rho(modelForIntegration.selectedTrim)
        fprintf(2,"\nRequested air density is different from the one used to compute the trim solution for DLM correction\n")
        fprintf(2,"In time marching LCO function, while checking the model set up by the preprocessor\n")
    end
end

if length(options.selectionTrim)>1
    error("Only one trim condition can be used at a time for the introduction of flight loads")
end

%% Define simulations to be performed

modelForIntegration.rho = options.rho;

% Construct all the possible combinations of nominal stiffnesses
stiffnessCombinations = combvec(options.kNominal{:});
nstiffnessCombinations = size(stiffnessCombinations,2);

% Construct all the possible combinations of gap sizes, only done in case
% of all static gaps. It makes no sense otherwise
if dynamicNonlinearities
    ngapCombinations = 1;
    for i = 1:length(options.gap)
        sortedGapValues = sort(options.gap{i});
        if strcmp(options.gapBehaviour{i,3},"both")
            options.gapVector{i} = [sortedGapValues, sortedGapValues(end-1:-1:1)];
        elseif strcmp(options.gapBehaviour{i,3},"increasing")
            options.gapVector{i} = sortedGapValues;
        else
            options.gapVector{i} = sortedGapValues(end:-1:1);
        end
    end
else
    gapCombinations = combvec(options.gap{:});
    ngapCombinations = size(gapCombinations,2);
end

% Construct the speed vector
if strcmp(options.speedBehaviour,"both")
    speedFolder = 'increasingAndDecreasingSpeed';
    options.speedVector = [options.speedVector, options.speedVector(end-1:-1:1)];
    options.nFFTwindows = options.nFFTwindows*2 - 1;
elseif strcmp(options.speedBehaviour,"increasing")
    speedFolder = 'increasingSpeed';
else
    speedFolder = 'decreasingSpeed';
    options.speedVector = options.speedVector(end:-1:1);
end

% Total time of a single run. Obtained considering the initialisation time,
% the number of windows, their length and the overlapping time between windows.
options.totalWindowLength = options.FFTwindowLength-options.overlapWindows;
options.simulationTime = options.initialisationTime + (options.nFFTwindows*options.totalWindowLength) + (options.overlapWindows*double(options.overlapWindows>0));

% Initialise vectors
LCOfrequency.fVect = [];
LCOfrequency.pVect = [];
LCOfrequency = repmat(LCOfrequency, nstiffnessCombinations, ngapCombinations, options.nFFTwindows);
LCOamplitude = repmat({zeros(size(options.gapPoints,1),3)}, nstiffnessCombinations, ngapCombinations, options.nFFTwindows);
% These are the peak values
LCOtorque = repmat({zeros(size(options.gapPoints,1),3)}, nstiffnessCombinations, ngapCombinations, options.nFFTwindows);
LCOmonitor = repmat({zeros(size(options.monitorPoints,3),1)}, nstiffnessCombinations, ngapCombinations, options.nFFTwindows);

home = pwd;
mkdir(speedFolder)
chdir(speedFolder)
homeSpeed = pwd;

nKs = options.nFFTwindows;

%% Introduce a possible mechanical preload (like a constant torque on a hinge)

modelForIntegration.constantStruForce = zeros(size(modelForIntegration.reducedBasis.V,2),1);
if options.introduceStruLoads
    loadDOF = obtainDOF(options.struLoads, modelForIntegration.model);
    modelForIntegration.constantStruForce = modelForIntegration.reducedBasis.V'*modelForIntegration.struData.Tgz(loadDOF,:)'*cellfun(@(x) x,options.struLoads(:,4));
end

%% Perform a TRIM analysis to compute the trim loads

% We perform the trim analysis using the stiff system, as we have free
% dofs otherwise. This is ok as the displacements of the dof are usually
% very small. We remain in the region of linearity of trim
% Otherwise we would need to recompute the AICs in the deformed
% configuration which is really expensive form the computational point of
% view (it would be an iterative procedure)
%
% We must perform one analysis per each stiffness combination

modelForIntegration.constantAeroForce = zeros(size(modelForIntegration.reducedBasis.V,2),nstiffnessCombinations);

if options.introduceFlightLoads
    for i = 1:nstiffnessCombinations
        if isempty(modelForIntegration.globalOptions.trim.ID)
            error('Requested the introduction of flight loads, but no info about trim provided')
        end

        trimOptions.selectionList = options.selectionTrim;
        trimOptions.outputType = options.trimType;
        trimOptions.hmomSet = 'full';

        [model_stiff, ~] = addNonlinearityStiffness(modelForIntegration.model, options.gapPoints, stiffnessCombinations(:,i));
        struData_stiff = structuralPreprocessor(options.fidScreen, model_stiff, []);

        [resultsTrim(i), modelForIntegration.aeroData] = solve_lin_trim(options.fidScreen, model_stiff, struData_stiff, modelForIntegration.aeroData, modelForIntegration.globalOptions.trim, trimOptions);
        modelForIntegration.constantAeroForce(:,i) = modelForIntegration.reducedBasis.V'*[resultsTrim(i).FdTotal; resultsTrim(i).hingeMomStru; zeros(size(struData_stiff.Tgz,2)-length(resultsTrim(i).FdTotal)-length(resultsTrim(i).hingeMomStru),1)];
    end

    % Check if the Mach used to compute flight loads is the same used for the
    % DLM matrices. Check also if the density is the same
    trimData = selectTrimCondition(modelForIntegration.globalOptions.trim, options.selectionTrim);
    if isnan(trimData.Q)
        [rho, ~, ~, a, ~] = ISA_h(trimData.z);
        Vinf = trimData.Mach*a;
		qinf = 0.5*rho*Vinf^2;
    else
        qinf = trimData.Q;
        rho = 2*qinf;
    end
    if modelForIntegration.aeroData.dlmData.aero.M(modelForIntegration.machUsed)~=trimData.Mach
        error("The requested mach number for the time marching integration is different from the Mach number used to compute the constant forces via trim solution.")
    end
    if options.rho~=rho
        fprintf(2,"\nRequested air density is different from the one used to compute the trim solution for the introduction of loads\n")
        fprintf(2,"In the time marching LCO while introducing flight loads\n")
    end

    % Forces are to be used per unit dynamic pressure
    modelForIntegration.constantAeroForce = modelForIntegration.constantAeroForce/qinf;

    clear trimData rho Vinf a qinf
end

%% The actual integration

%The par loop can be exchanged with one of the inner ones, with no modification to the code, depending on the case at hand.
for i = 1:nstiffnessCombinations

    mkdir(strcat('StiffnessCombination',num2str(i)))
    chdir(strcat('StiffnessCombination',num2str(i)))
    homeStiffness = pwd;

    fileid = fopen('README.txt','w');
    fprintf(fileid,"The values of nominal stiffnesses at the gap points are: ");
    fprintf(fileid,num2str(stiffnessCombinations(:,i)));
    fclose(fileid);

    % Before to perform the integration in time, we can produce some useful
    % output, if requested
    if options.plotMaxRealEigSS
        plotMaxRealEigSS(modelForIntegration,stiffnessCombinations(:,i),options.gapPoints,options)
    end

    if options.plotFixedPoints && exist("resultsTrim","var")
        plotFixedPoints(modelForIntegration,stiffnessCombinations(:,i),options.gapPoints,resultsTrim(i),options)
    end

    for j = 1:ngapCombinations

        qInitial = zeros(2*modelForIntegration.nStru,1);

        if ~dynamicNonlinearities
            mkdir(strcat('GapCombination',num2str(j)))
            chdir(strcat('GapCombination',num2str(j)))

            fileid = fopen('README.txt','w');
            fprintf(fileid,"The values of gaps at the gap points are: ");
            fprintf(fileid,num2str(gapCombinations(:,j)));
            fclose(fileid);
            % Note that this is a pseudo inverse
            qInitial(modelForIntegration.nStru+1:end) = modelForIntegration.Ugap\gapCombinations(:,j)/2*options.initialAmplitude;
        else
            fileid = fopen('README.txt','a');
            fprintf(fileid,"\nA dynamic gap was requested, thus only a single time marching simulation is performed");
            fprintf(fileid,"\nIf more combinations of static gap are of interest, set all of them to static in the inputs");
            fprintf(fileid,"\nIn this way the code will compute all the possible combinations and perform a simulation per each of them");
            fclose(fileid);
            % Note that this is a pseudo inverse
            initialGapVector = cellfun(@(x) x(1),options.gapVector);
            qInitial(modelForIntegration.nStru+1:end) = modelForIntegration.Ugap\initialGapVector/2*options.initialAmplitude;
        end

        homeGap=pwd;

        for k = 1:nKs

            mkdir(strcat("outputWindow",num2str(k)))
            chdir(strcat("outputWindow",num2str(k)))

            % Time windows
            tmin = options.initialisationTime + options.totalWindowLength*(k-1) - options.overlapWindows*double(options.overlapWindows<0);
            tmax = tmin + options.FFTwindowLength;

            if k ~= 1
                % Indeces to use of the old simulation
                indexes = out.tout>= tmin;
                outOld.modes = out.modes.signals.values(indexes,:);
                outOld.torque = out.torque.signals.values(indexes,:);
                outOld.rotationInTime = out.rotationInTime.signals.values(indexes,:);
                outOld.dynamicPressure = out.dynamicPressure.signals.values(indexes);
                outOld.tout = out.tout(indexes);
            else
                outOld.modes = [];
                outOld.torque = [];
                outOld.rotationInTime = [];
                outOld.dynamicPressure = [];
                outOld.tout = [];
            end

            load_system(strcat(home,filesep,options.simulinkModel))
            workSpace = get_param(options.simulinkModel, 'ModelWorkspace');

            if ~dynamicNonlinearities
                % The gap is peak to peak, thus we divide by 2
                assignin(workSpace, "freeplay", gapCombinations(:,j)/2);
            else
                assignin(workSpace, "freeplayVector", options.gapVector);
            end

            if k == 1
                assignin(workSpace, "qInitial", qInitial);
                assignin(workSpace, "qInitialAero", zeros(modelForIntegration.nAero,1));
            else
                assignin(workSpace, "qInitial", [out.modesVelocity.signals.values,out.modes.signals.values(end,:)].');
                assignin(workSpace, "qInitialAero", out.modesAero.signals.values.');
            end

            % We need to assign the variables to the workspace because of the par loop
            assignin(workSpace, "stiffness", stiffnessCombinations(:,i));

            assignin(workSpace, "constantAeroForce", modelForIntegration.constantAeroForce(:,i));
            assignin(workSpace, "constantStruForce", modelForIntegration.constantStruForce);

            assignin(workSpace, "A", modelForIntegration.A);
            assignin(workSpace, "Ca", modelForIntegration.Ca);
            assignin(workSpace, "B0", modelForIntegration.B0);
            assignin(workSpace, "B1", modelForIntegration.B1);
            assignin(workSpace, "B2", modelForIntegration.B2);
            assignin(workSpace, "D0", modelForIntegration.D0);
            assignin(workSpace, "D1", modelForIntegration.D1);
            assignin(workSpace, "D2", modelForIntegration.D2);
            assignin(workSpace, "nAero", modelForIntegration.nAero);
            assignin(workSpace, "lref", modelForIntegration.lref);
            assignin(workSpace, "rho", modelForIntegration.rho);

            assignin(workSpace, "Ugap", modelForIntegration.Ugap);
            assignin(workSpace, "nStru", modelForIntegration.nStru);
            assignin(workSpace, "invMC", modelForIntegration.invMC);
            assignin(workSpace, "invMK", modelForIntegration.invMK);
            assignin(workSpace, "Mhh", modelForIntegration.Mhh);

            assignin(workSpace, "options", options);

            % Run the actual simulation
            if k == 1
                out = sim(options.simulinkModel,tmax);
            else
                out = sim(options.simulinkModel,[out.tout(end),tmax]);
                % Remove initial condition as this is already present in the
                % previous output
                out.modes.signals.values = out.modes.signals.values(2:end,:);
                out.torque.signals.values = out.torque.signals.values(2:end,:);
                out.rotationInTime.signals.values = out.rotationInTime.signals.values(2:end,:);
                out.dynamicPressure.signals.values = out.dynamicPressure.signals.values(2:end,:);
                out.unsteadyAerodynamicForces.signals.values = out.unsteadyAerodynamicForces.signals.values(2:end,:);
                out.quasisteadyAerodynamicForces.signals.values = out.quasisteadyAerodynamicForces.signals.values(2:end,:);
                out.tout = out.tout(2:end);
            end

            % Restrict window in case of negative overlap or initial time
            indexes = arrayfun(@(x) x>= tmin && x<tmax, out.tout);
            % The second condition basically removes the last point

            % Merge results
            processTorque = [outOld.torque;out.torque.signals.values(indexes,:)];
            processModes = [outOld.modes;out.modes.signals.values(indexes,:)];
            processRotation = [outOld.rotationInTime;out.rotationInTime.signals.values(indexes,:)];
            processTime = [outOld.tout;out.tout(indexes)];
            processDynamicPressure = [outOld.dynamicPressure;out.dynamicPressure.signals.values(indexes)];

            % Compute main outputs
            [LCOamplitude{i,j,k},LCOfrequency(i,j,k),LCOtorque{i,j,k},LCOmonitor{i,j,k}]=computeTimeMarchingOutputs(processRotation, processTime, processTorque, processModes, modelForIntegration.Umonitor, options);

            % Main graphical outputs
            plotTimeMarchingOutputs(processTime, processRotation, processTorque, processModes, processDynamicPressure, LCOfrequency(i,j,k), modelForIntegration, options)

            % Save raw data
            parsave('rawData.mat',out.modes.signals.values.',...
                out.unsteadyAerodynamicForces.signals.values.',out.quasisteadyAerodynamicForces.signals.values.',...
                out.torque.signals.values.',out.rotationInTime.signals.values.',out.dynamicPressure.signals.values.',out.tout.'),

            close_system(options.simulinkModel,0)

            chdir(homeGap)
        end
        out=[];
        outOld=[];
        chdir(homeStiffness)
    end
    chdir(homeSpeed)
end
chdir(home)

%% Export the output

results.LCOamplitude = LCOamplitude;
results.LCOfrequency = LCOfrequency;
results.LCOtorque = LCOtorque;
results.LCOmonitor = LCOmonitor;
results.stiffnessCombinations = stiffnessCombinations;
results.gapCombinations = gapCombinations;

chdir("..")

return
