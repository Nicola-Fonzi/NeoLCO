clear all
close all
clc

%% General inputs

% Input file
filename_sma = 'flutter.dat';

% Solver info
ver = get_neolco_version(true);

% Model specific values
hingeScalarPoint = 54000;
stiffFlutterSpeed = 23.35;
kNominal = 2.025;
airDensity = 2.425;

% Inputs for the definition of the common base
optimalBaseOptions.nModes = 3;
optimalBaseOptions.fictitiousMass = [1e-4 1e-2 1e0 1e2 1e4];
optimalBaseOptions.stiffnessFactor = [0.00001 0.001 0.1 0.5 1];
optimalBaseOptions.optimumKnown = true;
optimalBaseOptions.gapPoints = {hingeScalarPoint,"s",0,"Aileron hinge"};
optimalBaseOptions.kNominal = {kNominal};

% Inputs for the aero database
aeroDatabaseOptions.selectionTrim = 'all';
aeroDatabaseOptions.DynVLM = false;
aeroDatabaseOptions.DynVLMtype = 'unsteady';

% Input for the time marching
preprocessTimeMarchingOptions.DynVLM = false;
preprocessTimeMarchingOptions.DynVLMtype = aeroDatabaseOptions;
preprocessTimeMarchingOptions.selectionTrim = 1; % This selection is used for the dyn VLM
preprocessTimeMarchingOptions.simulinkModel = 'SISO';
preprocessTimeMarchingOptions.gapPoints = optimalBaseOptions.gapPoints;
preprocessTimeMarchingOptions.gapBehaviour = {'freeplay','static',''}; %{'freeplay','dynamic','increasing'}; % the simulink model has to be built around the definition of static or dynamic gap, it cannot be decided later on.
preprocessTimeMarchingOptions.monitorPoints = {100,'g',5,'Pitch';...
                                               100,'g',3,'Plunge'};
timeMarchingOptions.gap = {2.3/180*pi}; %{[1.15, 2.3]/180*pi};
timeMarchingOptions.kNominal = optimalBaseOptions.kNominal;
%timeMarchingOptions.gapInterpolationType = 0;
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
describingFunctionsOptions.gapPoints = optimalBaseOptions.gapPoints;
describingFunctionsOptions.DynVLM = false;
describingFunctionsOptions.selectionTrim = 1; % This selection is used for the dyn VLM
describingFunctionsOptions.recomputeBase = false;
describingFunctionsOptions.searchQuenchPoint = true;
describingFunctionsOptions.maxKeq = kNominal;
describingFunctionsOptions.maxNKeq = 25;
describingFunctionsOptions.Vmax = 18;
describingFunctionsOptions.Vmin = 5;
describingFunctionsOptions.Vstep = 0.25;
describingFunctionsOptions.method = 'PK0';
describingFunctionsOptions.rho = airDensity;
describingFunctionsOptions.modesPlot = 1:3;
describingFunctionsOptions.axesUsed = 'body';

%% Loading of the model

inputData = readSmartcadFile(filename_sma);

%% Generate common basis

[reducedBasis, struData, model, struData_stiff, model_stiff, options] = obtainOptimalBase(inputData, optimalBaseOptions);

%% Generate aerodynamic database

[aeroData, aeroDatabaseOptions] = buildFullAeroData(model, model_stiff, struData_stiff, options, aeroDatabaseOptions);

%% Different equivalent stifnesses

[describingFunctionResults, describingFunctionsOptions] = describingFunctions(model, struData, aeroData, options, reducedBasis, aeroDatabaseOptions, describingFunctionsOptions);

%% Time marching simulation

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

for i = 1:length(describingFunctionResults.KeqVect)
    for j = 1:length(kNominal)
        gapSizes = timeMarchingOptions.gap{1};
        for k = 1:length(gapSizes)
            FFF = describingFunctionResults.KeqVect(i)/kNominal(j);
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
        plotHysteresis(describingFunctionResults.speedVector/stiffFlutterSpeed,amplitude(:,j,k).'./gapSizes(k),'--','LineWidth',1.5)
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
plotHysteresis(describingFunctionResults.speedVector/stiffFlutterSpeed,describingFunctionResults.frequencyVector,'--','LineWidth',1.5)
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
clear legendTitle temp
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
        clear temp toPlot
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
saveas(gcf,"monitor2.fig")
