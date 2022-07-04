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
    string = num2str(figIndex);
    saveas(h,strcat("Keq",num2str(KeqVect(indexK)),"_",string,"_harmonic.fig"))
    close(h)
    if figIndex==3
        indexK = indexK+1;
    end
end

clear h handles

%% Loop to solve the constant part, if required
% We exploit the fact that the KeqVect is already computed from the loop
% before. A mechanical preload is included here automatically

if options.introduceFlightLoads
    % Check if the Mach used to compute flight loads is the same used for the
    % DLM matrices. Check also if the density is the same
    if isempty(globalOptions.trim.ID)
        error('Requested the introduction of flight loads, but no info about trim provided')
    end
    trimDataDF = selectTrimCondition(globalOptions.trim, options.selectionTrim);
    if aeroData.dlmData.aero.M(options.machUsed)~=trimDataDF.Mach
        error("The requested mach number for the time marching integration is different from the Mach number used to compute the constant forces via trim solution.")
    end
    trimDataDF.z = nan;

    speedVectorBias = options.Vmin:options.Vstep:options.Vmax;
    for i = 1:length(speedVectorBias)
        % Loop on the speed vector and on the values of equivalent stiffness
        trimDataDF.Q = 0.5*options.rho*speedVectorBias(i)^2;

        for j = 1:length(KeqVect)

            trimOptions.selectionList = options.selectionTrim;
            trimOptions.outputType = options.trimType;
            trimOptions.hmomSet = 'full';

            [model_stiff, ~] = addNonlinearityStiffness(model, options.gapPoints, KeqVect(j));
            struData_stiff = structuralPreprocessor(options.fidScreen, model_stiff, []);

            [resultsTrim, ~] = solve_lin_trim(options.fidScreen, model_stiff, struData_stiff, aeroData, trimDataDF, trimOptions);
            [~, Pos] = obtainDOF(options.gapPoints, model);
            biasVect(j,i) = resultsTrim.SDispl(Pos, 1);
        end
    end
    clear trimDataDF
end

%% Plot the quasi-linear bias results
for j = 1:length(KeqVect)
    plot(speedVectorBias,biasVect(:,j),'LineWidth',2)
    h = gcf;
    saveas(h,strcat("Keq",num2str(KeqVect(j)),"_bias.fig"))
    close(h)
end

%% Reconstruct LCO amplitude
kNominal = options.kNominal{1};
gap = options.gap{1};

for i = 1:length(kNominal)
    for j = 1:length(gap)
        if options.introduceFlightLoads
            % Range of values for A and B
            AmplitudeDB = gap(j)*linspace(1,5,1000);
            BiasDB = AmplitudeDB;
            % Static stiffness as a function of A and B. Matrices that have
            % a different A per each column, and B changing row-wise.
            beta = (gap(j)/2 - BiasDB.')./AmplitudeDB;
            gamma = (-gap(j)/2 - BiasDB.')./AmplitudeDB;
            BoverA = BiasDB.'./AmplitudeDB;
            Ks = kNominal(i)*( 1 + 1/pi./BoverA*( beta.*asin(beta) + sqrt(1-beta.^2) -gamma.*asin(gamma) - sqrt(1-gamma.^2) ) );
            % Dynamic stiffness as a function of A and B. Matrices that have
            % a different A per each column, and B changing row-wise.
            Kd = kNominal(i)*( 1 - 1/pi*( asin(beta) + beta.*sqrt(1-beta.^2) - asin(gamma) - gamma.*sqrt(1-gamma.^2) ) );
            for k = 1:length(speedVectorLCO)
                % First, we find the combinations A1,B1 that provides an
                % equivalent stiffness equal to the one that provokes
                % flutter at our speed
                A1 = [];
                B1 = [];
                K1 = [];
                for m = 1:length(BiasDB)
                    [a1, k1] = polyxpoly(AmplitudeDB, ones(1,length(KeqVect))*KeqVect(k), AmplitudeDB, Kd(m,:));
                    K1 = [K1; k1];
                    A1 = [A1; a1];
                    B1 = [B1; BiasDB(m)*ones(length(a1),1)];
                end
                % Then, we can plot the SR at that speed
                slicedBias = interp2(speedVectorBias, KeqVect, biasVect, speedVectorLCO(k), KeqVect);
                figure
                plot(KeqVect, slicedBias, "LineWidth", 2)
                hold on
                % The intersection of this curve with all the possible Ks
                % provides A2 and B2
                plot(Ks(:,1:100:end), BiasDB, "LineWidth", 2)
                legend([{"SR"}, string(num2cell(AmplitudeDB(1:100:end)))])
                A2 = [];
                B2 = [];
                K2 = [];
                for m = 1:length(AmplitudeDB)
                    [b2, k2] = polyxpoly(slicedBias, KeqVect, BiasDB, Ks(:,m));
                    K2 = [K2; k2];
                    A2 = [A2; AmplitudeDB(m)*ones(length(b2),1)];
                    B2 = [B2: b2];
                end
                % We can now plot A2,B2 to double check
                figure
                plot(A1,B1,'LineWidth',2)
                hold on
                plot(A2,B2,'LineWidth',2)
                % The intersection of this curve with the curve A1,B1
                % gives the possible solutions
                [A3,B3]=polyxpoly(A1,B1,A2,B2);
            end
        else
            amplitudeRatioDB = linspace(1/5,1,1000);
            for m = 1:length(amplitudeRatioDB)
                kRatioDB(m) = 1/pi*(pi - 2*asin(amplitudeRatioDB(m)) + sin(2*asin(amplitudeRatioDB(m)))) - ...
                    4/pi*i*cos(asin(amplitudeRatioDB(m)));
            end
            for k = 1:length(speedVectorLCO)
                Kd = KeqVect(k)/kNominal(i);
                % We intersect the describing function with the value of Kd
                % for which we have flutter
                amplitudeRatio = 1./interp1(kRatioDB,amplitudeRatioDB,Kd,'linear','extrap');
                if strcmp(options.amplitudeDefinition,'maxPeak')
                    LCOamplitude{i,j,k} = repmat(amplitudeRatio*gap(j)/2,1,3);
                else
                    LCOamplitude{i,j,k} = repmat(amplitudeRatio/sqrt(2)*gap(j)/2,1,3);
                end
                LCOfrequency(i,j,k)=frequencyLCO(k);
            end
        end
    end
end

results.KeqVect = KeqVect;
results.speedVector = speedVectorLCO;
results.LCOamplitude = LCOamplitude;
results.LCOfrequency = LCOfrequency;
results.bias = nan;
results.stiffnessCombinations = kNominal;
results.gapCombinations = gap;

return