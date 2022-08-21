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
    timeMarchingResults, timeMarchingOptions, preprocessTimeMarchingOptions, options)

% Set default options
iOpt = 0;
iOpt = iOpt+1; baseOpt.fidScreen = 1;                 descr{iOpt} = 'Fid for screen printing. [1].';
iOpt = iOpt+1; baseOpt.useInApp = 0;                  descr{iOpt} = 'Utility flag to be used when called from App';
iOpt = iOpt+1; baseOpt.plotForApp = 1;                descr{iOpt} = 'Index of the nonlinearity to be plotted also on App';
iOpt = iOpt+1; baseOpt.halfGapNormalisation = 1;      descr{iOpt} = 'Normalise the LCO amplitude using half gap size, 1=yes, 0=no. [1].';
iOpt = iOpt+1; baseOpt.singleGapDF = 1;               descr{iOpt} = 'Plot one gap only for the describing functions, 1=yes, 0=no. [1].';
iOpt = iOpt+1; baseOpt.normalisationSpeed = 1;        descr{iOpt} = 'Normalisation speed to be used in the plots. [1].';
iOpt = iOpt+1; baseOpt.monitorNormalisation = 0;      descr{iOpt} = 'Normalise also the amplitude at the monitor points, 1=yes, 0=no. [0],';
iOpt = iOpt+1; baseOpt.torqueNormalisation = 0;       descr{iOpt} = 'Normalise also the amplitude of the force at the nonlinearity points, 1=yes, 0=no. [0],';
iOpt = iOpt+1; baseOpt.plotType = 'single';           descr{iOpt} = 'Quantity to plot, single value ("single") or bifurcation plot ("bfc"). ["single"].';

if nargin==0
    printOptionDescription(baseOpt, descr);
    return
end

% Process input options
options = setOptions(baseOpt, 'error', options);

if options.useInApp
    allfigs = findall(0,'Type', 'figure');
    appHandle = findall(allfigs, 'Name', 'NeoLCO');
end

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
                gapToPlot = describingFunctionResults.gapCombinations(m,j)/(1+options.halfGapNormalisation);
                if strcmp(options.plotType,'single')
                    for k = 1:describingFunctionOptions.nKeq
                        ytoPlot(k) = subsref(describingFunctionResults.LCOamplitude{i,j,k}, struct('type', '()', 'subs', {{m, 3}}));
                    end
                    xtoPlot = describingFunctionResults.speedVector/options.normalisationSpeed;
                else
                    for k = 1:describingFunctionOptions.nKeq
                        ytoPlot(k) = subsref(describingFunctionResults.LCOamplitude{i,j,k}, struct('type', '()', 'subs', {{m, 1}}));
                    end
                    ytoPlot(k+1) = nan;
                    for k = describingFunctionOptions.nKeq+2:describingFunctionOptions.nKeq*2+1
                        ytoPlot(k) = subsref(describingFunctionResults.LCOamplitude{i,j,k-describingFunctionOptions.nKeq-1}, ...
                            struct('type', '()', 'subs', {{m, 2}}));
                    end
                    xtoPlot = [describingFunctionResults.speedVector, nan, describingFunctionResults.speedVector]/options.normalisationSpeed;
                end
                plot(xtoPlot(:),ytoPlot(:)/gapToPlot,'LineWidth',1.5)
                if options.singleGapDF
                    string = strcat(describingFunctionOptions.gapPoints{m,4},...
                        ' Stiffness ', num2str(describingFunctionResults.stiffnessCombinations(:,i).'),' DF');
                else
                    string = strcat(describingFunctionOptions.gapPoints{m,4},...
                        ' Gap ',num2str(describingFunctionResults.gapCombinations(:,j).'),' Stiffness ',...
                        num2str(describingFunctionResults.stiffnessCombinations(:,i).'),' DF');
                end
                legendTitle{m,index(m)} =  string;
                index(m)=index(m)+1;
            end
        end
    end
end

clear xtoPlot ytoPlot gapToPlot string

if options.useTM
    for m = 1:size(preprocessTimeMarchingOptions.gapPoints,1)  % Per each nonlinearity point
        figure(1000+m)
        hold on
        for i = 1:size(timeMarchingResults.stiffnessCombinations,2)  % Per each stiffness combination
            for j = 1:size(timeMarchingResults.gapCombinations,2)  % Per each gap combination
                gapToPlot = timeMarchingResults.gapCombinations(m,j)/(1+options.halfGapNormalisation);
                if strcmp(options.plotType,'single')
                    for k = 1:timeMarchingOptions.nFFTwindows
                        ytoPlot(k) = subsref(timeMarchingResults.LCOamplitude{i,j,k}, struct('type', '()', 'subs', {{m, 3}}));
                    end
                    xtoPlot = timeMarchingOptions.speedVector/options.normalisationSpeed;
                else
                    for k = 1:timeMarchingOptions.nFFTwindows
                        ytoPlot(k) = subsref(timeMarchingResults.LCOamplitude{i,j,k}, struct('type', '()', 'subs', {{m, 1}}));
                    end
                    ytoPlot(k+1) = nan;
                    for k = timeMarchingOptions.nFFTwindows+2:timeMarchingOptions.nFFTwindows*2+1
                        ytoPlot(k) = subsref(timeMarchingResults.LCOamplitude{i,j,k-timeMarchingOptions.nFFTwindows-1}, ...
                            struct('type', '()', 'subs', {{m, 2}}));
                    end
                    xtoPlot = [timeMarchingOptions.speedVector, nan, timeMarchingOptions.speedVector]/options.normalisationSpeed;
                end
                plotHysteresis(xtoPlot(:),ytoPlot(:)/gapToPlot,index(m))
                string = strcat(preprocessTimeMarchingOptions.gapPoints{m,4},...
                    ' Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',...
                    num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM');
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
        string = strcat(preprocessTimeMarchingOptions.gapPoints{m,4},".fig");
    end
    saveas(figure(1000+m),string)
    if options.useInApp
        if m==options.plotInApp
            fig = gcf;
            copyobj(fig.Children.Children,appHandle.Children(4).Children(6).Children(19))
        end
    end
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
            xtoPlot = describingFunctionResults.speedVector/options.normalisationSpeed;
            plot(xtoPlot(:),ytoPlot(:),'LineWidth',1.5)
            if options.singleGapDF
                string = strcat('Stiffness ', num2str(describingFunctionResults.stiffnessCombinations(:,i).'),' DF');
            else
                string = strcat('Gap ',num2str(describingFunctionResults.gapCombinations(:,j).'),' Stiffness ',...
                    num2str(describingFunctionResults.stiffnessCombinations(:,i).'),' DF');
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
            xtoPlot = timeMarchingOptions.speedVector/options.normalisationSpeed;
            plotHysteresis(xtoPlot(:),ytoPlot(:),index)
            string = strcat('Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',...
                num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM');
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
if options.useInApp
    fig = gcf;
    copyobj(fig.Children.Children,appHandle.Children(4).Children(6).Children(20))
end
saveas(figure(2000),string)

clear legendTitle

if options.useTM && isempty(timeMarchingResults.LCOmonitor{1,1,1})==0
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
                    if options.monitorNormalisation && size(preprocessTimeMarchingOptions.gapPoints,1)==1
                        gapToPlot = timeMarchingResults.gapCombinations(1,j)/(1+options.halfGapNormalisation);
                    elseif optiona.monitorNormalisation && size(preprocessTimeMarchingOptions.gapPoints,1)>1
                        error("The normalisation of the monitor points can only be done if a single nonlinearity is present")
                    else
                        gapToPlot = 1;
                    end
                    if strcmp(options.plotType,'single')
                        for k = 1:timeMarchingOptions.nFFTwindows
                            ytoPlot(k) = timeMarchingResults.LCOmonitor{i,j,k}(m,3);
                        end
                        xtoPlot = timeMarchingOptions.speedVector/options.normalisationSpeed;
                    else
                        for k = 1:timeMarchingOptions.nFFTwindows
                            ytoPlot(k) = timeMarchingResults.LCOmonitor{i,j,k}(m,1);
                        end
                        ytoPlot(k+1) = nan;
                        for k = timeMarchingOptions.nFFTwindows+2:timeMarchingOptions.nFFTwindows*2+1
                            ytoPlot(k) = timeMarchingResults.LCOmonitor{i,j,k-timeMarchingOptions.nFFTwindows-1}(m,2);
                        end
                        xtoPlot = [timeMarchingOptions.speedVector, nan, timeMarchingOptions.speedVector]/options.normalisationSpeed;
                    end
                    plotHysteresis(xtoPlot(:),ytoPlot(:)/gapToPlot,index(m))
                    string = strcat(preprocessTimeMarchingOptions.monitorPoints{m,4},...
                        ' Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',...
                        num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM');
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
        ylabel(preprocessTimeMarchingOptions.monitorPoints{m,4})
        xlabel('Speed [m/s]')
        string = strcat(preprocessTimeMarchingOptions.monitorPoints{m,4},".fig");
        if options.useInApp
            if m==options.plotInApp
                fig = gcf;
                copyobj(fig.Children.Children,appHandle.Children(4).Children(6).Children(21))
            end
        end
        saveas(figure(3000+m),string)
    end

    clear legendTitle

    %
    % Torque
    %
    index(1:100)=1;

    if options.useTM
        for m = 1:size(preprocessTimeMarchingOptions.gapPoints,1)  % Per each monitor point
            figure(4000+m)
            hold on
            for i = 1:size(timeMarchingResults.stiffnessCombinations,2)  % Per each stiffness combination
                for j = 1:size(timeMarchingResults.gapCombinations,2)  % Per each gap combination
                    if options.torqueNormalisation
                        gapToPlot = timeMarchingResults.gapCombinations(m,j)/(1+options.halfGapNormalisation);
                    else
                        gapToPlot = 1;
                    end
                    if strcmp(options.plotType,'single')
                        for k = 1:timeMarchingOptions.nFFTwindows
                            ytoPlot(k) = squeeze(timeMarchingResults.LCOtorque{i,j,k}(m,3));
                        end
                        xtoPlot = timeMarchingOptions.speedVector/options.normalisationSpeed;
                    else
                        for k = 1:timeMarchingOptions.nFFTwindows
                            ytoPlot(k) = timeMarchingResults.LCOtorque{i,j,k}(m,1);
                        end
                        ytoPlot(k+1) = nan;
                        for k = timeMarchingOptions.nFFTwindows+2:timeMarchingOptions.nFFTwindows*2+1
                            ytoPlot(k) = timeMarchingResults.LCOtorque{i,j,k-timeMarchingOptions.nFFTwindows-1}(m,2);
                        end
                        xtoPlot = [timeMarchingOptions.speedVector, nan, timeMarchingOptions.speedVector]/options.normalisationSpeed;
                    end
                    plotHysteresis(xtoPlot(:),ytoPlot(:)/gapToPlot,index(m))
                    string = strcat(preprocessTimeMarchingOptions.gapPoints{m,4},...
                        ' Gap ',num2str(timeMarchingResults.gapCombinations(:,j).'),' Stiffness ',...
                        num2str(timeMarchingResults.stiffnessCombinations(:,i).'),' TM');
                    legendTitle{m,index(m)} =  string;
                    index(m)=index(m)+1;
                end
            end
        end
    end

    clear toPlot gapToPlot string

    for m = 1:size(legendTitle,1)
        figure(4000+m)
        legend(legendTitle(m,:))
        ylabel(strcat("Force at ",preprocessTimeMarchingOptions.gapPoints{m,4}))
        xlabel('Speed [m/s]')
        string = strcat("Force at ",preprocessTimeMarchingOptions.gapPoints{m,4},".fig");
        if options.useInApp
            if m==options.plotInApp
                fig = gcf;
                copyobj(fig.Children.Children,appHandle.Children(4).Children(6).Children(22))
            end
        end
        saveas(figure(4000+m),string)
    end

    clear legendTitle

end

return