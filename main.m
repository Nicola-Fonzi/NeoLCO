clear
close all
clc

%% General inputs

% INPUT FILE
filename_sma = 'D:\Desktop\PhD\X-DIA\ElevatorOnlyLCO\numericalSimulations\Final_DuringGVT_plate\NeoCASS\flutter.dat';

% SOLVER INFO
ver = get_neolco_version(true);

% ADDITIONAL VALUES FOR THE MODEL ANALYSIS
tipLE = 40002;
tipTE = 40001;
midspanLE = 40006;
midspanElevator = 40015;
hingeScalarPoint = 1;

% INPUTS FOR THE DESCRIBING FUNCTION
kNominal = 48;
maxKeq = 1.8;
maxNumberKeq = 25;
flutterOptions.Vmax = 70;
flutterOptions.Vmin = 5;
flutterOptions.Vstep = 0.5;
flutterOptions.method = 'PK0';
flutterOptions.rho = 1.13;
modesToPlotDF = 1:8;

% INPUT FOR THE TIME MARCHING
preprocessTimeMarchingOptions.simulinkModel = 'SISO';
preprocessTimeMarchingOptions.gapPoints = {hingeScalarPoint,'s',0,'Elevator hinge'};
preprocessTimeMarchingOptions.gapBehaviour = {'freeplay','static',''}; % Unfortunately, the simulink model has to be built around the definition of static or dynamic gap, it cannot be decided later on.
preprocessTimeMarchingOptions.monitorPoints = {tipLE,'g',3,'Accelerometer at the tip LE';...
                                               tipTE,'g',3,'Accelerometer at the tip TE';...
                                               midspanLE,'g',3,'Accelerometer at midpan LE';...
                                               midspanElevator,'g',3,'Accelerometer at midpan elevator'};
                                           
timeMarchingOptions.gap = {[0.1, 0.2, 0.4, 0.6, 1, 1.5, 2, 2.5, 3]/180*pi};
timeMarchingOptions.kNominal = {kNominal};   % In the MIMO case this will be a vector
timeMarchingOptions.rho = flutterOptions.rho;
timeMarchingOptions.speedVector = (14:2:50);
timeMarchingOptions.speedBehaviour = 'both'; % The behaviour of speed can be easily changed without rebuilding the model
timeMarchingOptions.speedInterpolationType = 0;
timeMarchingOptions.nFFTwindows = length(timeMarchingOptions.speedVector);
timeMarchingOptions.FFTwindowLength = 20;
timeMarchingOptions.overlapWindows = -10;
timeMarchingOptions.initialisationTime = 10;
timeMarchingoptions.amplitudeDefinition = 'std';
timeMarchingoptions.introduceFlightLoads = false;
timeMarchingoptions.trimType = 'grounded';
timeMarchingoptions.selectionTrim = 1;

% FLAGS DEFINING THE BEHAVIOUR OF THE CODE
recomputeBase = false;
searchFlutterQuenchingPoint = true;
optimalBaseOptions.NMODES = [8 12 16 20];
optimalBaseOptions.FM = [1e-4 1e-2 1e0 1e2 1e4];
optimalBaseOptions.KN = [0.00001 0.4 0.8 1.2 5 kNominal];
optimalBaseOptions.OptimumKnown = true;
optimalBaseOptions.FICTMASS_CMASS_ID = 34006;
HINGE_CELAS_ID = 34106;
optimalBaseOptions.HINGE_CELAS_ID = HINGE_CELAS_ID;


%% Loading of the models
inputData = readSmartcadFile(filename_sma);

inputData_stiff = inputData;
inputData_stiff.CELAS.ID = [inputData_stiff.CELAS.ID HINGE_CELAS_ID];
inputData_stiff.CELAS.value = [inputData_stiff.CELAS.value kNominal];
inputData_stiff.CELAS.Node = [inputData_stiff.CELAS.Node [hingeScalarPoint;0]];
inputData_stiff.CELAS.DOF = [inputData_stiff.CELAS.DOF [0;0]];
inputData_stiff.CELAS.PID = [inputData_stiff.CELAS.PID 0];

[model, options] = processInputData(inputData);

[model_stiff, options_stiff] = processInputData(inputData_stiff);

%% Gather missing options from the model
struOpt = [];

eigOpt = options.eig;

flutterOptions.axesUsed = 'body';
flutterOptions.Mlist_dlm = options.aero.Mlist_dlm;
flutterOptions.klist = options.aero.klist;

fid = options.FID;

%% Preprocessing

struData = structuralPreprocessor(fid, model, struOpt);

struData_stiff = structuralPreprocessor(fid, model_stiff, struOpt);

%% Recover of the DOF indeces required for postprocessing

[IDtable, gridDOF, spointDOF] = getIDtable(model.Node.ID, model.Spoint.ID);

HingePos = find(model.Spoint.ID == hingeScalarPoint,1);
HingeDOF = spointDOF(HingePos);


%% Generate common basis

reducedBasis = [];

if recomputeBase==0

    [reducedBasis, resultsFolder] = obtainOptimalBase(inputData, model, struData, options, model_stiff, struData_stiff, optimalBaseOptions, HingeDOF);

    reducedBasis.Bmm = modalDamp(model, reducedBasis.V'*struData.Mzz*reducedBasis.V, reducedBasis.V'*struData.Kzz*reducedBasis.V,'Real');

    mkdir("TimeMarching")
    chdir("TimeMarching")

    mkdir(resultsFolder)
    chdir(resultsFolder)

end


%% Different equivalent stifnesses

[flutterSpeed, flutterFrequency, KeqVect, resultsFlutter, aeroData] = describingFunctionsPK(maxKeq, maxNumberKeq, inputData, HINGE_CELAS_ID, hingeScalarPoint, struOpt, fid, eigOpt, ...
    recomputeBase, flutterOptions, reducedBasis, searchFlutterQuenchingPoint, modesToPlotDF);

if recomputeBase
    return
end


%% Here we try to reproduce the same using a time marching simulation

[modelForIntegration, preprocessTimeMarchingOptions] = preprocessTimeMarchingLCO(model, struData, preprocessTimeMarchingOptions, reducedBasis, aeroData, options);

[timeMarchingResults, timeMarchingOptions]=timeMarchingLCO(modelForIntegration, timeMarchingOptions, preprocessTimeMarchingOptions);

%% Now, we plot

% First, at each flutter speed we extract the amplitude from the Keq
% Create database
amplitudeRatioDB = linspace(1/5,1,1000);
index = 1;
for i = amplitudeRatioDB
    kRatioDB(index) = 1/pi*(pi - 2*asin(i) + sin(2*asin(i))) - 4/pi*i*cos(asin(i));
    index = index + 1;
end

for i = 1:length(KeqVect)
    for j = 1:length(kNominal)
        gapSizes = timeMarchingOptions.gap{1};
        for k = 1:length(gapSizes)
            FFF = KeqVect(i)/kNominal(j);
            amplitudeRatio = 1./interp1(kRatioDB,amplitudeRatioDB,FFF,'linear','extrap');
            if strcmp(timeMarchingOptions.amplitudeDefinition,'maxPeak')
                amplitude(i,j,k) = amplitudeRatio*gapSizes(k)/2;
            else
                amplitude(i,j,k) = amplitudeRatio/sqrt(2)*gapSizes(k)/2;
            end
        end
    end
end

figure
hold on
index=1;
for j = 1:length(kNominal)
    for k = 1:length(gapSizes)
        plotHysteresis(flutterSpeed/stiffFlutterSpeed,amplitude(:,j,k).'./gapSizes(k),'--','LineWidth',1.5)
        legendTitle{index} = ['Gap ',num2str(gapSizes(k)),' Kn ',num2str(kNominal(j)),' DF'] ;
        index=index+1;
    end
end

for i = 1:size(timeMarchingResults.stiffnessCombinations,2)
    for j = 1:size(timeMarchingResults.gapCombinations,2)
        for k = 1:timeMarchingOptions.nFFTwindows
            temp = timeMarchingResults.LCOamplitude{i,j,k};
            toPlot(k) = temp(3);
        end
        plotHysteresis(timeMarchingOptions.speedVector/stiffFlutterSpeed,toPlot./...
            timeMarchingResults.gapCombinations(:,j),'LineWidth',1.5)
        legendTitle{index} = ['Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM'] ;
        index=index+1;
    end
end
toPlot = [];
legend(legendTitle)
ylabel('\delta/\delta_{Gap}')
xlabel('Speed [m/s]')
saveas(gcf,"amplitude.fig")

figure
hold on
clear legendTitle
index=1;
plotHysteresis(flutterSpeed/stiffFlutterSpeed,flutterFrequency,'--','LineWidth',1.5)
legendTitle{index} = 'DF';
index=index+1;

for i = 1:size(timeMarchingResults.stiffnessCombinations,2)
    for j = 1:size(timeMarchingResults.gapCombinations,2)
        for k = 1:timeMarchingOptions.nFFTwindows
            [dummy , freq_index] = max(timeMarchingResults.LCOfrequency(i,j,k).pVect);
            freq(k) = timeMarchingResults.LCOfrequency(i,j,k).fVect(freq_index);
        end
        plotHysteresis(timeMarchingOptions.speedVector/stiffFlutterSpeed,freq,'LineWidth',1.5)
        legendTitle{index} = ['Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM'] ;
        index=index+1;
    end
end
legend(legendTitle)
ylabel('Frequency [Hz]')
xlabel('Speed [m/s]')
saveas(gcf,"frequency.fig")

figure
hold on
clear legendTitle
index=1;
for i = 1:size(timeMarchingResults.stiffnessCombinations,2)
    for j = 1:size(timeMarchingResults.gapCombinations,2)
        for k = 1:timeMarchingOptions.nFFTwindows
            temp = timeMarchingResults.LCOmonitor{i,j,k};
            toPlot(:,k) = temp(:,3);
        end
        plotHysteresis(timeMarchingOptions.speedVector/stiffFlutterSpeed,toPlot./...
            timeMarchingResults.gapCombinations(:,j),'LineWidth',1.5)
        legendTitle{index} = ['Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM'] ;
        index=index+1;
    end
end
legend(legendTitle)
ylabel('Monitor')
xlabel('Speed [m/s]')
saveas(gcf,"monitor.fig")

figure
hold on
clear legendTitle
index=1;
for i = 1:size(timeMarchingResults.stiffnessCombinations,2)
    for j = 1:size(timeMarchingResults.gapCombinations,2)
        for k = 1:timeMarchingOptions.nFFTwindows
            temp = timeMarchingResults.LCOmonitor{i,j,k};
            toPlot(:,k) = temp(:,1);
        end
        plotHysteresis(timeMarchingOptions.speedVector/stiffFlutterSpeed,toPlot./...
            timeMarchingResults.gapCombinations(:,j),'LineWidth',1.5)
        legendTitle{index} = ['Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM UpperLimit'] ;
        index=index+1;
        for k = 1:timeMarchingOptions.nFFTwindows
            temp = timeMarchingResults.LCOmonitor{i,j,k};
            toPlot(:,k) = temp(:,2);
        end
        plotHysteresis(timeMarchingOptions.speedVector/stiffFlutterSpeed,toPlot./...
            timeMarchingResults.gapCombinations(:,j),'LineWidth',1.5)
        legendTitle{index} = ['Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM LowerLimit'] ;
        index=index+1;
    end
end
legend(legendTitle)
ylabel('Monitor')
xlabel('Speed [m/s]')
saveas(gcf,"monitor.fig")
