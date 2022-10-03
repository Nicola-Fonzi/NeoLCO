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
function results = describingFunctionPK(model, struData, aeroData, reducedBasis, globalOptions, aeroDatabaseOptions, options)

%% Initialise variables

KeqVect = [];
indexK = 1;
deltaK = options.maxKeq/options.nKeq;
Keq_prev = -deltaK;
reachedMinKeq = false;

if length(options.selectionTrim)>1
    error("Only one trim condition can be specified to be used for the introduction of flight loads")
end
if length(options.selectionTrimDynVLM)>1
    error("Only one trim condition can be specified to be used for the correction via VLM of the DLM matrices")
end

%% Loop to solve the flutter equations for the harmonic part

while true

    Keq = Keq_prev+deltaK;

    [model_descrFun, ~] = addNonlinearityStiffness(model, options.gapPoints, Keq);

    % We need to rebuild the matrices as we changed the stiffness
    struData_descrFun = structuralPreprocessor(options.fidScreen, model_descrFun, options.struOpt);

    % recompute the basis
    if options.recomputeBase
        options.eigOpt.UseFictmass = false;
        resultsEig_descrFun = solve_eig_fun(options.fidScreen, model_descrFun, struData_descrFun, options.eigOpt);
        reducedBasis = defineReducedBasis(struData_descrFun, resultsEig_descrFun, 'all');
        if ~isempty(model.Tabdmp1.ID(model.Tabdmp1.ID==model.Tabdmp1.caseControl))
            KDAMP = model.Tabdmp1.KDAMP;
            if KDAMP < 0
                error("Only real damping is supported in this context, please use KDAMP>0")
            else
                fprintf(options.fidScreen,' - Building real damping matrix for viscous damping...');
                reducedBasis.Bmm = modalDamp(model, reducedBasis.V'*struData.Mzz*reducedBasis.V, ...
                    reducedBasis.V'*struData.Kzz*reducedBasis.V,'Real');
            end
            fprintf(options.fidScreen, 'done.\n');
        end
        if isfield(aeroData,"aeroMatrix_dlm")
            aeroData = rmfield(aeroData,"aeroMatrix_dlm");
        end
    end

    % Flutter computation
    if options.DynVLM
        [resultsFlutter_descrFun, aeroData] = solve_linflutt_fun(model_descrFun, struData_descrFun, reducedBasis, aeroData, options, globalOptions.trim, aeroDatabaseOptions);
    else
        [resultsFlutter_descrFun, aeroData] = solve_linflutt_fun(model_descrFun, struData_descrFun, reducedBasis, aeroData, options);
    end

    if options.searchQuenchPoint
        if any(any(real(resultsFlutter_descrFun.modes.eig)>=0)) || reachedMinKeq || Keq==0
            resultsFlutter(indexK).modes = resultsFlutter_descrFun.modes;
            indexK = indexK+1;
            KeqVect = [KeqVect, Keq];
            Keq_prev = Keq;
            deltaK = options.maxKeq/options.nKeq;
            reachedMinKeq = false;
        else
            deltaK = deltaK/2;
            if deltaK<(options.maxKeq/options.nKeq/1000)
                reachedMinKeq = true;
                deltaK = options.maxKeq/options.nKeq;
            end
        end
    else
        resultsFlutter(indexK).modes = resultsFlutter_descrFun.modes;
        indexK = indexK+1;
        KeqVect = [KeqVect, Keq];
        Keq_prev = Keq;
    end
    if indexK>options.nKeq
        break
    end
end

%% Plot the quasi-linear harmonic results

plotOption.modeRequested = options.modesPlot;
plotOption.plotEnvelope = true;
plotOption.envelopeLabel = "Equivalent stiffness";
plotOption.envelopeMeasureUnit = "Nm";
plotOption.envelopeList = KeqVect;

[speedVectorLCO, frequencyLCO] = plotVgDiagrams(resultsFlutter, plotOption);

handles=findall(0,'type','figure');

h = figure(handles(1).Number);
saveas(h,"EnvelopeHarmonicLCO.fig")
close(h)

indexK = 1;
for i = 1:handles(2).Number
    h = figure(i);
    figIndex = mod(i,3);
    if figIndex==0
        figIndex=3;
    end
    saveas(h,strcat("Keq",num2str(KeqVect(indexK)),"_",num2str(figIndex),"_harmonic.fig"))
    close(h)
    if figIndex==3
        indexK = indexK+1;
    end
end

clear h handles

%% Loop to solve the constant part, if required

if options.introduceFlightLoads || options.introduceStruLoads || options.introduceGravityLoads

    % First, we introduce GRAV or LOAD cards
    if options.introduceGravityLoads
        model.grav.SID = 1;
        model.grav.direction = reshape(options.gravDirection,[3,1]);
        model.grav.amplitude = 9.8;
        trimOptions.gravSet = 1;
    end
    if options.introduceStruLoads
        model.load.SID = 1;
        nPointsLoaded = size(options.struLoads,1);
        model.load.amplitude = cellfun(@(x) x, options.struLoads(:,4));
        model.load.directionType = zeros(1, nPointsLoaded);
        model.load.directionType(cellfun(@(x) x=='s', options.struLoads(:,2))) = -1;
        model.load.directionData = zeros(4, nPointsLoaded);
        for loadedPoint = 1:nPointsLoaded
            if any(options.struLoads{loadedPoint,3} == [1,4])
                model.load.directionData(1,loadedPoint) = 1;
            elseif any(options.struLoads{loadedPoint,3} == [2,5])
                model.load.directionData(2,loadedPoint) = 1;
            else
                model.load.directionData(3,loadedPoint) = 1;
            end
        end
        model.load.isMoment = zeros(1, nPointsLoaded);
        model.load.isMoment(cellfun(@(x) x==4 || x==5 || x==6, options.struLoads(:,3))) = 1;
        model.load.node = cellfun(@(x) x, options.struLoads(:,1));
        model.load.followerType = zeros(1, nPointsLoaded);
        trimOptions.loadSet = 1;
    end

    % Check if the Mach used to compute flight loads is the same used for the
    % DLM matrices. Check also if the density is the same
    if options.introduceFlightLoads
        if isempty(globalOptions.trim.ID)
            error('Requested the introduction of flight loads, but no info about trim provided')
        end
    else
        fprintf(options.fidScreen, "In order to introduce the structural loads and/or gravity loads, a trim card\n");
        fprintf(options.fidScreen, "is added with no angle of attack, no side-slip angle, no rotations.");
        globalOptions.trim.ID = 1;
        globalOptions.trim.name = 'noAero';
        globalOptions.trim.type = 0;
        globalOptions.trim.Mach = aeroData.dlmData.aero.M(options.machUsed);
        globalOptions.trim.z = nan;
        globalOptions.trim.Q = 1; % It is overwritten next
        globalOptions.trim.aelinkIDlist = 'all'; % All the surface links are used
        % If the model is free with gravity or loads it may have acceleration, it is thus taken into account
        globalOptions.trim.data(1).labelList = {'SIDES','ANGLEA','ROLL','PITCH','YAW'};
        globalOptions.trim.data(1).valueList = [0,0,0,0,0];
        options.selectionTrim = 1;
    end
    trimDataDF = selectTrimCondition(globalOptions.trim, options.selectionTrim);
    if aeroData.dlmData.aero.M(options.machUsed)~=trimDataDF.Mach
        error("The requested mach number for the descrbing function analysis is different from the Mach number used to compute the aerodynamic matrices.")
    end
    trimDataDF.z = nan;

    speedVectorBias = options.Vmin:options.Vstep:options.Vmax;
    KeqVectBias = 0:(options.maxKeq/options.nKeq):options.maxKeq;
    for i = 1:length(speedVectorBias)
        % Loop on the speed vector and on the values of equivalent stiffness
        trimDataDF.Q = 0.5*options.rho*speedVectorBias(i)^2;

        for j = 1:length(KeqVectBias)

            trimOptions.selectionList = options.selectionTrim;
            trimOptions.outputType = options.trimType;
            trimOptions.hmomSet = 'full';
            trimOptions.fidScreen = options.fidScreen;

            [model_stiff, ~] = addNonlinearityStiffness(model, options.gapPoints, KeqVectBias(j));
            struData_stiff = structuralPreprocessor(options.fidScreen, model_stiff, []);

            [resultsTrim, ~] = solve_lin_trim(options.fidScreen, model_stiff, struData_stiff, aeroData, trimDataDF, trimOptions);
            [~, Pos] = obtainDOF(options.gapPoints, model);
            biasVect(j,i) = resultsTrim.SDispl(Pos, 1);
        end
    end
    clear trimDataDF

%% Plot the quasi-linear bias results
    for j = 1:length(KeqVectBias)
        plot(speedVectorBias,biasVect(j,:),'LineWidth',2)
        xlabel("V [m/s]","FontSize",12)
        ylabel("Bias","FontSize",12)
        grid minor
        h = gcf;
        saveas(h,strcat("Keq",num2str(KeqVectBias(j)),"_bias.fig"))
        close(h)
    end

    if options.introduceGravityLoads
        model.grav.SID = [];
        model.grav.amplitude = [];
        model.grav.direction = zeros(3,0);
    end
    if options.introduceStruLoads
        model.load.SID = [];
        model.load.node = [];
        model.load.ismoment = [];
        model.load.followerType = [];
        model.load.amplitude = [];
        model.load.directionData = zeros(4,0);
        model.load.directionType = [];
    end

end

%% Reconstruct LCO amplitude
kNominal = options.kNominal{1};
gap = options.gap{1};

for i = 1:length(kNominal)
    for j = 1:length(gap)
        if options.introduceFlightLoads || options.introduceStruLoads || options.introduceGravityLoads
            % Range of values for A and B
            AmplitudeDB = gap(j)/2*linspace(1e-12,5,1000);
            BiasDB = gap(j)/2*linspace(-5,5,1000);
            % Static stiffness as a function of A and B. Matrices that have
            % a different A per each column, and B changing row-wise.
            beta = (gap(j)/2 - BiasDB.')./AmplitudeDB;
            gamma = (-gap(j)/2 - BiasDB.')./AmplitudeDB;
            BoverA = BiasDB.'./AmplitudeDB;

            Ks = kNominal(i)./BoverA.*( -(gamma+beta)/2 - g(gamma) + g(beta));
            % Dynamic stiffness as a function of A and B. Matrices that have
            % a different A per each column, and B changing row-wise.
            Kd = kNominal(i)*( 1 + (f(gamma) - f(beta)));
            for k = 1:length(KeqVect)
                if ~isnan(speedVectorLCO(k))

                    % First, we find the combinations A1,B1 that provides an
                    % equivalent stiffness equal to the one that provokes
                    % flutter at our speed
                    figure
                    plot(AmplitudeDB, ones(1,length(AmplitudeDB))*KeqVect(k))
                    hold on
                    plot(AmplitudeDB, Kd(1:100:end,:))
                    xlabel("LCO amplitude [rad]","FontSize",12)
                    ylabel("Dynamic stiffness","FontSize",12)
                    legend([{"Equivalent stiffness for flutter"},...
                        strcat("B=",string(num2cell(BiasDB(1:100:end))))],"FontSize",12)
                    h = gcf;
                    saveas(h,strcat("DynamicSolutionsKnominal",num2str(kNominal(i)),...
                        "Gap",num2str(gap(j)),"KeqDyn",num2str(KeqVect(k)),".fig"))
                    close(h)
                    A1 = [];
                    B1 = [];
                    K1 = [];
                    for m = 1:length(BiasDB)
                        [a1, k1] = polyxpoly(AmplitudeDB, ones(1,length(AmplitudeDB))*KeqVect(k), AmplitudeDB, Kd(m,:));
                        K1 = [K1; k1];
                        A1 = [A1; a1];
                        B1 = [B1; BiasDB(m)*ones(length(a1),1)];
                    end

                    % Then, we can plot the SR at that speed
                    slicedBias = interp2(speedVectorBias, KeqVectBias, biasVect, speedVectorLCO(k), KeqVectBias);
                    figure
                    plot(KeqVectBias, slicedBias, "LineWidth", 2)
                    hold on
                    % The intersection of this curve with all the possible Ks
                    % provides A2 and B2
                    plot(Ks(:,1:100:end), BiasDB, "LineWidth", 2)
                    ylabel("LCO bias [rad]","FontSize",12)
                    xlabel("Static stiffness","FontSize",12)
                    legend([{"Static Response"}, strcat("A=",string(num2cell(AmplitudeDB(1:100:end))))],"FontSize",12)
                    h = gcf;
                    saveas(h,strcat("StaticSolutionsKnominal",num2str(kNominal(i)),...
                        "Gap",num2str(gap(j)),"KeqDyn",num2str(KeqVect(k)),".fig"))
                    close(h)
                    A2 = [];
                    B2 = [];
                    K2 = [];
                    for m = 1:length(AmplitudeDB)
                        [b2, k2] = polyxpoly(slicedBias, KeqVectBias , BiasDB, Ks(:,m));
                        K2 = [K2; k2];
                        A2 = [A2; AmplitudeDB(m)*ones(length(b2),1)];
                        B2 = [B2; b2];
                    end

                    % We can now find the intersection
                    [A3,B3]=polyxpoly(A1,B1,A2,B2);
                    figure
                    plot(A1,B1,'LineWidth',2)
                    hold on
                    plot(A2,B2,'LineWidth',2)
                    plot(A3,B3,'ro','LineWidth',2)
                    xlabel("A","FontSize",12)
                    ylabel("B","FontSize",12)
                    legend("Dynamic solutions","Static solutions")
                    h = gcf;
                    saveas(h,strcat("IntersectSolutionsKnominal",num2str(kNominal(i)),...
                        "Gap",num2str(gap(j)),"KeqDyn",num2str(KeqVect(k)),".fig"))
                    close(h)
                    if isempty(A3)
                        LCOamplitude{i,j,k} = nan(1,3);
                        LCObias{i,j,k} = nan;
                        LCOfrequency(i,j,k) = nan;
                    else
                        LCOamplitude{i,j,k} = zeros(1,3, length(A3));
                        if strcmp(options.amplitudeDefinition,'maxPeak')
                            LCOamplitude{i,j,k}(1,1,:) = B3 + A3;
                            LCOamplitude{i,j,k}(1,2,:) = B3 - A3;
                            LCOamplitude{i,j,k}(1,3,:) = A3;
                        elseif strcmp(options.amplitudeDefinition,'rms')
                            LCOamplitude{i,j,k}(1,1,:) = B3 + A3/sqrt(2);
                            LCOamplitude{i,j,k}(1,2,:) = B3 - A3/sqrt(2);
                            LCOamplitude{i,j,k}(1,3,:) = sqrt(B3.^2 + (A3/sqrt(2)).^2);
                        else
                            LCOamplitude{i,j,k}(1,1,:) = A3/sqrt(2);
                            LCOamplitude{i,j,k}(1,2,:) = -A3/sqrt(2);
                            LCOamplitude{i,j,k}(1,3,:) = A3/sqrt(2);
                        end
                        LCObias{i,j,k} = B3;
                        LCOfrequency(i,j,k) = frequencyLCO(k);
                    end
                else
                    LCOamplitude{i,j,k} = nan(1,3);
                    LCObias{i,j,k} = nan;
                    LCOfrequency(i,j,k) = nan;
                end
            end
        else
            amplitudeRatioDB = linspace(1/5,1,1000);
            for m = 1:length(amplitudeRatioDB)
                kRatioDB(m) = 1/pi*(pi - 2*asin(amplitudeRatioDB(m)) + sin(2*asin(amplitudeRatioDB(m)))) - ...
                    4/pi*amplitudeRatioDB(m)*cos(asin(amplitudeRatioDB(m)));
            end
            for k = 1:length(KeqVect)
                if ~isnan(speedVectorLCO(k))
                    Kd = KeqVect(k)/kNominal(i);
                    % We intersect the describing function with the value of Kd
                    % for which we have flutter
                    amplitudeRatio = 1./interp1(kRatioDB,amplitudeRatioDB,Kd,'linear','extrap');
                    if strcmp(options.amplitudeDefinition,'maxPeak')
                        LCOamplitude{i,j,k} = repmat(amplitudeRatio*gap(j)/2,1,3);
                        LCOamplitude{i,j,k}(2) = -LCOamplitude{i,j,k}(2);
                    else
                        LCOamplitude{i,j,k} = repmat(amplitudeRatio/sqrt(2)*gap(j)/2,1,3);
                        LCOamplitude{i,j,k}(2) = -LCOamplitude{i,j,k}(2);
                   end
                    LCOfrequency(i,j,k)=frequencyLCO(k);
                    LCObias{i,j,k} = nan;
                else
                    LCOamplitude{i,j,k} = nan(1,3);
                    LCObias{i,j,k} = nan;
                    LCOfrequency(i,j,k) = nan;
                end
            end
        end
    end
end

results.KeqVect = KeqVect;
results.speedVector = speedVectorLCO;
results.LCOamplitude = LCOamplitude;
results.LCOfrequency = LCOfrequency;
results.LCObias = LCObias;
results.stiffnessCombinations = kNominal;
results.gapCombinations = gap;

end

function value = f(gamma)

i = gamma < -1;
j = gamma > 1;
k = ~(i | j);

value = nan(size(gamma));
value(i) = -1/2;
value(j) = 1/2;
value(k) = 1/pi*(real(asin(gamma(k))) + gamma(k).*sqrt(1-gamma(k).^2));

end

function value = g(gamma)

i = gamma < -1;
j = gamma > 1;
k = ~(i | j);

value = nan(size(gamma));
value(i) = abs(gamma(i))/2;
value(j) = abs(gamma(j))/2;
value(k) = 1/pi*(gamma(k).*real(asin(gamma(k))) + sqrt(1-gamma(k).^2));

end