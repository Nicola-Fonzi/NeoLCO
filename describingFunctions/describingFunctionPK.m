function results = describingFunctionPK(model, struData, aeroData, reducedBasis, globalOptions, aeroDatabaseOptions, options)

KeqVect = [];
indexK = 1;
deltaK = options.maxKeq/options.maxNKeq;
Keq_prev = -deltaK;
reachedMinKeq = false;

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
        if any(any(real(resultsFlutter_descrFun.modes.eig)>=0)) || reachedMinKeq
            resultsFlutter(indexK).modes = resultsFlutter_descrFun.modes;
            indexK = indexK+1;
            KeqVect = [KeqVect, Keq];
            Keq_prev = Keq;
            deltaK = options.maxKeq/options.maxNKeq;
            reachedMinKeq = false;
        else
            deltaK = deltaK/2;
            if deltaK<(options.maxKeq/options.maxNKeq/1000)
                reachedMinKeq = true;
                deltaK = options.maxKeq/options.maxNKeq;
            end
        end
    else
        resultsFlutter(indexK).modes = resultsFlutter_descrFun.modes;
        indexK = indexK+1;
        KeqVect = [KeqVect, Keq];
        Keq_prev = Keq;
    end
    if indexK>options.maxNKeq
        break
    end
end

plotOption.modeRequested = options.modesPlot;
plotOption.plotEnvelope = true;
plotOption.envelopeLabel = "Equivalent stiffness";
plotOption.envelopeMeasureUnit = "Nm";
plotOption.envelopeList = KeqVect;

[flutterSpeed, flutterFrequency] = plotVgDiagrams(resultsFlutter, plotOption);

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

results.KeqVect = KeqVect;
results.speedVector = flutterSpeed;
results.frequencyVector = flutterFrequency;

return