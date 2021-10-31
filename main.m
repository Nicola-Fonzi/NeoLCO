clear
close all
clc

%% General inputs

% Input file
filename_sma = 'D:\Desktop\PhD\X-DIA\ElevatorOnlyLCO\numericalSimulations\Final_DuringGVT_plate\correlated\steelWeightSteelConnection\NeoCASS\flutter.dat';

% Solver info
ver = get_neolco_version(true);

% Model specific values
tipLE = 40002;
tipTE = 40001;
midspanLE = 40006;
midspanElevator = 40015;
hingeScalarPoint = 1;
kNominal = 48;
airDensity = 1.13;

% Inputs for the definition of the common base
optimalBaseOptions.nModes = [8 12 16 20];
optimalBaseOptions.fictitiousMass = [1e-4 1e-2 1e0 1e2 1e4];
optimalBaseOptions.stiffnessFactor = [0.00001 0.001 0.1 0.5 1];
optimalBaseOptions.optimumKnown = false;
optimalBaseOptions.gapPoints = {hingeScalarPoint,'s',0,'Elevator hinge'};
optimalBaseOptions.kNominal = {kNominal};

% Inputs for the aero database
aeroDatabaseOptions.selectionTrim = 'all';
aeroDatabaseOptions.DynVLM = false;
aeroDatabaseOptions.DynVLMtype = 'unsteady';

% Inputs for the time marching approach
preprocessTimeMarchingOptions.DynVLM = false;
preprocessTimeMarchingOptions.DynVLMtype = aeroDatabaseOptions;
preprocessTimeMarchingOptions.selectionTrim = 1;
preprocessTimeMarchingOptions.simulinkModel = 'SISO';
preprocessTimeMarchingOptions.gapPoints = optimalBaseOptions.gapPoints;
preprocessTimeMarchingOptions.gapBehaviour = {'freeplay','static',''};
preprocessTimeMarchingOptions.monitorPoints = {tipLE,'g',3,'Accelerometer at the tip LE';...
                                               tipTE,'g',3,'Accelerometer at the tip TE';...
                                               midspanLE,'g',3,'Accelerometer at midpan LE';...
                                               midspanElevator,'g',3,'Accelerometer at midpan elevator'};                                           
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
timeMarchingOptions.amplitudeDefinition = 'std';
timeMarchingOptions.introduceFlightLoads = false;
timeMarchingOptions.trimType = 'grounded';
timeMarchingOptions.selectionTrim = 1;

% Inputs for the describing function approach
recomputeBase = false;
searchFlutterQuenchingPoint = true;
maxKeq = 1.8;
maxNumberKeq = 25;
flutterOptions.Vmax = 70;
flutterOptions.Vmin = 5;
flutterOptions.Vstep = 0.5;
flutterOptions.method = 'PK0';
flutterOptions.rho = airDensity;
modesToPlotDF = 1:8;
flutterOptions.axesUsed = 'body';

%% Loading of the model

inputData = readSmartcadFile(filename_sma);

%% Generate common basis

[reducedBasis, struData, model, struData_stiff, model_stiff, resultsFolder, options] = obtainOptimalBase(inputData, optimalBaseOptions);

reducedBasis.Bmm = modalDamp(model, reducedBasis.V'*struData.Mzz*reducedBasis.V, reducedBasis.V'*struData.Kzz*reducedBasis.V,'Real');

mkdir(resultsFolder)
chdir(resultsFolder)

%% Generate aerodynamic database

[aeroData, options] = buildFullAeroData(model, model_stiff, struData_stiff, options, aeroDatabaseOptions);

%% Different equivalent stifnesses

flutterOptions.Mlist_dlm = options.aero.Mlist_dlm;
flutterOptions.klist = options.aero.klist;

if recomputeBase
    reducedBasis = [];
end

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
        plotHysteresis(flutterSpeed,amplitude(:,j,k).'./gapSizes(k)*2,'--','LineWidth',1.5)
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
        plotHysteresis(timeMarchingOptions.speedVector,toPlot./...
            timeMarchingResults.gapCombinations(:,j)*2,'LineWidth',1.5)
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
plotHysteresis(flutterSpeed,flutterFrequency,'--','LineWidth',1.5)
legendTitle{index} = 'DF';
index=index+1;

for i = 1:size(timeMarchingResults.stiffnessCombinations,2)
    for j = 1:size(timeMarchingResults.gapCombinations,2)
        for k = 1:timeMarchingOptions.nFFTwindows
            [dummy , freq_index] = max(timeMarchingResults.LCOfrequency(i,j,k).pVect);
            freq(k) = timeMarchingResults.LCOfrequency(i,j,k).fVect(freq_index);
        end
        plotHysteresis(timeMarchingOptions.speedVector,freq,'LineWidth',1.5)
        legendTitle{index} = ['Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM'] ;
        index=index+1;
    end
end
legend(legendTitle)
ylabel('Frequency [Hz]')
xlabel('Speed [m/s]')
saveas(gcf,"frequency.fig")

colorTable = [ ...
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.0000    0.5000    0.0000
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0.0000    0.0000    1.0000
    1.0000    0.0000    0.0000
    0.0000    1.0000    0.0000
    0.6000    0.6000    0.6000
    0.0000    0.0000    0.0000
    0.5000    0.0000    0.8000
    0.0000    0.4000    0.4000
    ];

clear legendTitle temp
index=1;
iColor = 1;
for i = 1:size(timeMarchingResults.stiffnessCombinations,2)
    for j = 1:size(timeMarchingResults.gapCombinations,2)
        for k = 1:timeMarchingOptions.nFFTwindows
            temp = timeMarchingResults.LCOmonitor{i,j,k};
            toPlotUp(:,k) = temp(:,1);
            toPlotDown(:,k) = temp(:,2);
        end
        for y = 1:size(toPlotUp,1)
            figure(y+100)
            hold on
            plotHysteresis(timeMarchingOptions.speedVector,toPlotUp(y,:),'LineWidth',1.5,'color',colorTable(iColor,:))
            plotHysteresis(timeMarchingOptions.speedVector,toPlotDown(y,:),'LineWidth',1.5,'color',colorTable(iColor,:))
        end
        iColor = iColor+1;
        legendTitle{index} = ['Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM UpperLimit'] ;
        index=index+1;
        legendTitle{index} = ['Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM LowerLimit'] ;
        index=index+1;
    end
end
for y = 1:size(toPlotUp,1)
    legend(legendTitle)
    ylabel(preprocessTimeMarchingOptions.monitorPoints(y,4))
    xlabel('Speed [m/s]')
    saveas(figure(y+100),strcat(preprocessTimeMarchingOptions.monitorPoints(y,4),".fig"))
end
