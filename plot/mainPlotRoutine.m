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
function mainPlotRoutine(describingFunctionResults, describingFunctionOptions, ...
    timeMarchingResults, timeMarchingOptions, options)

% Set default options
iOpt = 0;
iOpt = iOpt+1; baseOpt.fidScreen = 1;                 descr{iOpt} = 'Fid for screen printing. [1].';
iOpt = iOpt+1; baseOpt.halfGapNormalisation = 1;      descr{iOpt} = 'Normalise the LCO amplitude using half gap size, 1=yes, 0=no. [1].';
iOpt = iOpt+1; baseOpt.singleGapDF = 1;               descr{iOpt} = 'Plot one gap only for the describing functions, 1=yes, 0=no. [1].';
iOpt = iOpt+1; baseOpt.normalisationSpeed = 1;        descr{iOpt} = 'Normalisation speed to be used in the plots. [1].';
iOpt = iOpt+1; baseOpt.monitorNormalisation = 1;      descr{iOpt} = 'Normalise also the amplitude at the monitor points, 1=yes, 0=no. [1],';
iOpt = iOpt+1; baseOpt.plotType = 'std';              descr{iOpt} = 'Quantity to plot, standard deviation ("std") or bifurcation plot ("bfc"). ["std"].';

if nargin==0
    printOptionDescription(baseOpt, descr);
    return
end

% Process input options
options = setOptions(baseOpt, 'error', options);

if isempty(describingFunctionResults) || isempty(describingFunctionOptions)
    options.useDF = 0;
else
    options.useDF = 1;
end
if isempty(timeMarchingResults) || isempty(timeMarchingOptions)
    options.useTM = 0;
else
    options.useTM = 1;
end

%
% Amplitude
%
index(1:100)=1;

if options.useDF
    for m = 1:size(describingFunctionOptions.gapPoints,1)  % Per each nonlinearity point
        figure(1000+m)
        hold on
        for i = 1:size(describingFunctionResults.stiffnessCombinations,2)  % Per each stiffness combination
            if options.singleGapDF
                maxj = 1;
            else
                maxj = size(describingFunctionResults.gapCombinations,2);
            end
            for j = 1:maxj  % Per each gap combination, or for one gap combination only
                gapToPlot = describingFunctionResults.gapCombinations(m,j)*(1+options.halfGapNormalisation);
                if options.plotType == 'std'
                    for k = 1:describingFunctionOptions.nKeq
                        ytoPlot(k) = describingFunctionResults.LCOamplitude{i,j,k}(m,3);
                    end
                    xtoPlot = describingFunctionResults.speedVector/options.normalisationSpeed;
                else
                    for k = 1:describingFunctionOptions.nKeq
                        ytoPlot(k) = describingFunctionResults.LCOamplitude{i,j,k}(m,1);
                    end
                    ytoPloy(k+1) = nan;
                    for k = describingFunctionOptions.nKeq+2:describingFunctionOptions.nKeq*2+1
                        ytoPlot(k) = describingFunctionResults.LCOamplitude{i,j,k-describingFunctionOptions.nKeq-1}(m,2);
                    end
                    xtoPlot = [describingFunctionResults.speedVector, nan, describingFunctionResults.speedVector]/options.normalisationSpeed;
                end
                plot(xtoPlot,ytoPlot/gapToPlot,'LineWidth',1.5)
                if options.singleGapDF
                    string = [describingFunctionOptions.gapPoints{m,4},...
                        ' Stiffness ', num2str(describingFunctionResults.stiffnessCombinations(:,i).'),' DF'];
                else
                    string = [describingFunctionOptions.gapPoints{m,4},...
                        ' Gap ',num2str(describingFunctionResults.gapCombinations(:,j).'),' Stiffness ',...
                        num2str(describingFunctionResults.stiffnessCombinations(:,i).'),' DF'];
                end
                legendTitle{m,index(m)} =  string;
                index(m)=index(m)+1;
            end
        end
    end
end

clear toPlot gapToPlot string

if options.useTM
    for m = 1:size(preprocessTimeMarchingOptions.gapPoints,1)  % Per each nonlinearity point
        figure(1000+m)
        hold on
        for i = 1:size(timeMarchingResults.stiffnessCombinations,2)  % Per each stiffness combination
            for j = 1:size(timeMarchingResults.gapCombinations,2)  % Per each gap combination
                gapToPlot = timeMarchingResults.gapCombinations(m,j)*(1+options.halfGapNormalisation);
                if options.plotType == 'std'
                    for k = 1:timeMarchingOptions.nFFTwindows
                        ytoPlot(k) = timeMarchingResults.LCOamplitude{i,j,k}(m,3);
                    end
                    xtoPlot = timeMarchingResults.speedVector/options.normalisationSpeed;
                else
                    for k = 1:timeMarchingOptions.nFFTwindows
                        ytoPlot(k) = timeMarchingResults.LCOamplitude{i,j,k}(m,1);
                    end
                    ytoPloy(k+1) = nan;
                    for k = timeMarchingOptions.nFFTwindows+2:timeMarchingOptions.nFFTwindows*2+1
                        ytoPlot(k) = timeMarchingResults.LCOamplitude{i,j,k-timeMarchingOptions.nFFTwindows-1}(m,2);
                    end
                    xtoPlot = [timeMarchingResults.speedVector, nan, timeMarchingResults.speedVector]/options.normalisationSpeed;
                end
                plotHysteresis(xtoPlot,ytoPlot/gapToPlot,index(m))
                string = [timeMarchingOptions.gapPoints{m,4},...
                    ' Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',...
                    num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM'];
                legendTitle{m,index(m)} =  string;
                index(m)=index(m)+1;
            end
        end
    end
end

clear toPlot gapToPlot string

for m = 1:size(legendTitle,1)
    figure(1000+m)
    legend(legendTitle(m,:))
    ylabel('\delta/\delta_{Gap}')
    xlabel('Speed [m/s]')
    try
        string = strcat(describingFunctionOptions.gapPoints{m,4},".fig");
    catch
        string = strcat(timeMarchingOptions.gapPoints{m,4},".fig");
    end
    saveas(figure(1000+m),string)
end

clear legendTitle

%
% Frequency
%
index=1;

if options.useDF
    figure(2000)
    hold on
    for i = 1:size(describingFunctionResults.stiffnessCombinations,2)  % Per each stiffness combination
        if options.singleGapDF
            maxj = 1;
        else
            maxj = size(describingFunctionResults.gapCombinations,2);
        end
        for j = 1:maxj  % Per each gap combination, or for one gap combination only
            ytoPlot = describingFunctionResults.LCOfrequency(i,j,:);
            plot(describingFunctionResults.speedVector/options.normalisationSpeed,ytoPlot,'LineWidth',1.5)
            if options.singleGapDF
                string = ['Stiffness ', num2str(describingFunctionResults.stiffnessCombinations(:,i).'),' DF'];
            else
                string = ['Gap ',num2str(describingFunctionResults.gapCombinations(:,j).'),' Stiffness ',...
                    num2str(describingFunctionResults.stiffnessCombinations(:,i).'),' DF'];
            end
            legendTitle{index} =  string;
            index=index+1;
        end
    end
end

clear toPlot string

if options.useTM
    figure(2000)
    hold on
    for i = 1:size(timeMarchingResults.stiffnessCombinations,2)  % Per each stiffness combination
        for j = 1:size(timeMarchingResults.gapCombinations,2)  % Per each gap combination
            for k = 1:timeMarchingOptions.nFFTwindows
                [~ , freq_index] = max(timeMarchingResults.LCOfrequency(i,j,k).pVect);
                ytoPlot(k) = timeMarchingResults.LCOfrequency(i,j,k).fVect(freq_index);
            end
            plotHysteresis(timeMarchingResults.speedVector/options.normalisationSpeed,ytoPlot,index)
            string = ['Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',...
                num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM'];
            legendTitle{index} =  string;
            index=index+1;
        end
    end
end

clear toPlot string

figure(2000)
legend(legendTitle)
ylabel('Frequency [Hz]')
xlabel('Speed [m/s]')
string = "frequency.fig";
saveas(figure(2000),string)

clear legendTitle

%
% Monitor
%
index(1:100)=1;

if options.useTM
    for m = 1:size(preprocessTimeMarchingOptions.monitorPoints,1)  % Per each monitor point
        figure(3000+m)
        hold on
        for i = 1:size(timeMarchingResults.stiffnessCombinations,2)  % Per each stiffness combination
            for j = 1:size(timeMarchingResults.gapCombinations,2)  % Per each gap combination
                if options.monitorNormalisation
                    gapToPlot = timeMarchingResults.gapCombinations(m,j)*(1+options.halfGapNormalisation);
                else
                    gapToPlot = 1;
                end
                if options.plotType == 'std'
                    for k = 1:timeMarchingOptions.nFFTwindows
                        ytoPlot(k) = timeMarchingResults.LCOmonitor{i,j,k}(m,3);
                    end
                    xtoPlot = timeMarchingResults.speedVector/options.normalisationSpeed;
                else
                    for k = 1:timeMarchingOptions.nFFTwindows
                        ytoPlot(k) = timeMarchingResults.LCOmonitor{i,j,k}(m,1);
                    end
                    ytoPloy(k+1) = nan;
                    for k = timeMarchingOptions.nFFTwindows+2:timeMarchingOptions.nFFTwindows*2+1
                        ytoPlot(k) = timeMarchingResults.LCOmonitor{i,j,k-timeMarchingOptions.nFFTwindows-1}(m,2);
                    end
                    xtoPlot = [timeMarchingResults.speedVector, nan, timeMarchingResults.speedVector]/options.normalisationSpeed;
                end
                plotHysteresis(timeMarchingResults.speedVector/options.normalisationSpeed,ytoPlot/gapToPlot,index(m))
                string = [timeMarchingOptions.monitorPoints{m,4},...
                    ' Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',...
                    num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM'];
                legendTitle{m,index(m)} =  string;
                index(m)=index(m)+1;
            end
        end
    end
end

clear toPlot gapToPlot string

for m = 1:size(legendTitle,1)
    figure(3000+m)
    legend(legendTitle(m,:))
    ylabel(timeMarchingOptions.monitorPoints{m,4})
    xlabel('Speed [m/s]')
    string = strcat(timeMarchingOptions.monitorPoints{m,4},".fig");
    saveas(figure(3000+m),string)
end

clear legendTitle

return