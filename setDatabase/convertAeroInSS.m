function [solution, chosenModes, machUsed, k, Ha, selectedTrim] = convertAeroInSS(aeroData, model, struData, reducedBasis, options)

if options.DynVLM
    trimData = aeroData.vlmData.DynVLM.trimRes.data;
    selectedTrim = trimData.ID == options.selectionTrim;
    % We need to be sure that the same Mach number is used for both the VLM
    % matrices and the DLM matrices to be corrected
    if ~isempty(options.machNumber)
        if options.machNumber ~= trimData.Mach(selectedTrim)
            error("Requested mach number for DLM which is different from the Mach number used for the trim solution")
        end
    end
    machUsed = aeroData.dlmData.aero.M==trimData.Mach(selectedTrim);

    aeroSystem.Qhh = aeroData.aeroMatrix_dlm.aeroMatrix.Qhh(:,:,:,machUsed); % NOTE: This is inertial frame, a conversion is needed in theory
    aeroSystem = mergeVLMandDLMmatrices(model, aeroSystem, aeroData.vlmData, selectedTrim, struData, reducedBasis, options.DynVLMtype);
    
else
    selectedTrim = [];
    if isempty(options.machNumber)
        % If not specified we use the first Mach number
        machUsed = 1;
    else
        machUsed = aeroData.dlmData.aero.M==options.machNumber;
    end
    aeroSystem.Qhh = aeroData.aeroMatrix_dlm.aeroMatrix.Qhh(:,:,:,machUsed);
end
k = aeroData.aeroMatrix_dlm.aero.k;

chosenModes = [];
for modeIndex = 1:size(aeroSystem.Qhh,1)
    plotLinearDispl(model,[],struData,reducedBasis,modeIndex,0.1);
    use = input("Use this mode?");
    while isempty(use)
      disp("Wrong selection, please use [1] if you want to use the mode, or [0] if not")
      use = input("Use this mode?");
    end  
    while use~=1 && use~=0
      disp("Wrong selection, please use [1] if you want to use the mode, or [0] if not")
      use = input("Use this mode?");
    end
    if use
        chosenModes = [chosenModes modeIndex];
        h = gcf;
        saveas(h,strcat("Mode",num2str(modeIndex),"--used.fig"));
        close
    else
        h = gcf;
        saveas(h,strcat("Mode",num2str(modeIndex),".fig"));
        close
    end
end

Ha = aeroSystem.Qhh(chosenModes,chosenModes,:);

solution = improvedMFDfun(k,Ha,options.opt,options.optsLM,options.eigsopt,options.algROM,[],[]);

return