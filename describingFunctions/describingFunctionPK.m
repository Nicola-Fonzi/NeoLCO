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

%% Loop to solve the flutter equations

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
        selectedTrim = aeroData.vlmData.DynVLM.trimRes.data.ID == options.selectionTrim;
        aeroData.vlmData.DynVLM.DynMatrices.Qa = aeroData.vlmData.DynVLM.DynMatrices.Qa(:,:,:,:,selectedTrim);
        trimData = selectTrimCondition(globalOptions.trim, options.selectionTrim);
        [resultsFlutter_descrFun, aeroData] = solve_linflutt_fun(model_descrFun, struData_descrFun, reducedBasis, aeroData, options, trimData, aeroDatabaseOptions);
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

%% Plot the quasi-linear flutter results

plotOption.modeRequested = options.modesPlot;
plotOption.plotEnvelope = true;
plotOption.envelopeLabel = "Equivalent stiffness";
plotOption.envelopeMeasureUnit = "Nm";
plotOption.envelopeList = KeqVect;

[speedVector, frequency] = plotVgDiagrams(resultsFlutter, plotOption);

handles=findall(0,'type','figure');

h = figure(handles(1).Number);
saveas(h,"Envelope.fig")
close(h)

indexK = 1;
for i = 1:handles(2).Number
    h = figure(i);
    figIndex = mod(i,3);
    if figIndex==0
        figIndex=3;
    end
    string = num2str(figIndex);
    saveas(h,strcat("Keq",num2str(KeqVect(indexK)),"_",string,".fig"))
    close(h)
    if figIndex==3
        indexK = indexK+1;
    end
end

clear h handles

%% Reconstruct LCO amplitude
amplitudeRatioDB = linspace(1/5,1,1000);
index = 1;
for i = amplitudeRatioDB
    kRatioDB(index) = 1/pi*(pi - 2*asin(i) + sin(2*asin(i))) - 4/pi*i*cos(asin(i));
    index = index + 1;
end

kNominal = options.kNominal{1};
gap = options.gap{1};

for i = 1:length(kNominal)
    for j = 1:length(gap)
        for k = 1:length(KeqVect)
            FFF = KeqVect(k)/kNominal(i);
            amplitudeRatio = 1./interp1(kRatioDB,amplitudeRatioDB,FFF,'linear','extrap');
            if strcmp(options.amplitudeDefinition,'maxPeak')
                LCOamplitude{i,j,k} = repmat(amplitudeRatio*gap(j)/2,1,3);
            else
                LCOamplitude{i,j,k} = repmat(amplitudeRatio/sqrt(2)*gap(j)/2,1,3);
            end
            LCOfrequency(i,j,k)=frequency(k);
        end
    end
end

results.KeqVect = KeqVect;
results.LCOamplitude = LCOamplitude;
results.speedVector = speedVector;
results.LCOfrequency = LCOfrequency;
results.stiffnessCombinations = kNominal;
results.gapCombinations = gap;

return