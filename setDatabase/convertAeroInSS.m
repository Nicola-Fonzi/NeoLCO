%*******************************************************************************
%                                                                              *
%                    _   _            _     ____ ___                           *
%                   | \ | | ___  ___ | |   / ___/ _ \                          *
%                   |  \| |/ _ \/ _ \| |  | |  | | | |                         *
%                   | |\  |  __/ (_) | |__| |__| |_| |                         *
%                   |_| \_|\___|\___/|_____\____\___/                          *
%                                                                              *
%                                                                              *
% Copyright (C) 2020 - 2024                                                    *
%                                                                              *
% Nicola Fonzi (nicola.fonzi@outlook.com)                                      *
%                                                                              *
%                                                                              *
% This file is part of NeoLCO Software (github.com/Nicola-Fonzi/NeoLCO).       *
% The use of this software is licensed based on the licence distributed        *
% together with the source code. If you have not received the license please   *
% contact the copywright owner before using the software.                      *
%                                                                              *
%*******************************************************************************
function [solution, machUsed, k, Ha, selectedTrim, options] = convertAeroInSS(aeroData, model, struData, reducedBasis, options)

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

plotInApp = options.useInApp && isempty(options.chosenModes);
askModes = isempty(options.chosenModes);

if plotInApp
    allfigs = findall(0,'Type', 'figure');
    appHandle = findall(allfigs, 'Name', 'NeoLCO');
    h = appHandle.Children(4).Children(5).Children.Children(1).Children(25);
    h.Visible = 'on';
    appHandle.Children(4).Children(5).Children.Children(1).Children(13).Visible = 'on';
    appHandle.Children(4).Children(5).Children.Children(1).Children(7).Visible = 'on';
else
    h = figure();
end

for modeIndex = 1:size(aeroSystem.Qhh,1)
    plotLinearDispl(model,[],struData,reducedBasis,modeIndex,0.1,h);
    if askModes
        if plotInApp
            while(appHandle.Children(4).Children(5).Children.Children(1).Children(13).Value ~=1) && ...
                    (appHandle.Children(4).Children(5).Children.Children(1).Children(7).Value ~=1)
                pause(0.1);
            end
            use = appHandle.Children(4).Children(5).Children.Children(1).Children(13).Value;
            appHandle.Children(4).Children(5).Children.Children(1).Children(13).Value = 0;
            appHandle.Children(4).Children(5).Children.Children(1).Children(7).Value = 0;
        else
            use = input("Use this mode?");
            while isempty(use) || (use~=1 && use~=0)
                disp("Wrong selection, please use [1] if you want to use the mode, or [0] if not")
                use = input("Use this mode?");
            end
        end
        if use
            options.chosenModes = [options.chosenModes modeIndex];
            if options.useInApp
                dummyfig = figure;
                dummyax = axes;
                copyobj(h.Children,dummyax)
                saveas(dummyfig,strcat("Mode",num2str(modeIndex),"--used.fig"));
                close(dummyfig)
            else
                saveas(h,strcat("Mode",num2str(modeIndex),"--used.fig"));
            end
        else
            if options.useInApp
                dummyfig = figure;
                dummyax = axes;
                copyobj(h.Children,dummyax)
                saveas(dummyfig,strcat("Mode",num2str(modeIndex),".fig"));
                close(dummyfig)
            else
                saveas(h,strcat("Mode",num2str(modeIndex),".fig"));
            end
        end
        if plotInApp
            cla(h)
        else
            cla(h.Children)
        end
    else
        if any(modeIndex==options.chosenModes)
            saveas(h,strcat("Mode",num2str(modeIndex),"--used.fig"));
        else
            saveas(h,strcat("Mode",num2str(modeIndex),".fig"));
        end
        cla(h.Children)
    end
end


Ha = aeroSystem.Qhh(options.chosenModes,options.chosenModes,:);

if options.useInApp
    appHandle.Children(4).Children(5).Children.Children(1).Children(25).Visible = 'off';
    appHandle.Children(4).Children(5).Children.Children(1).Children(13).Visible = 'off';
    appHandle.Children(4).Children(5).Children.Children(1).Children(7).Visible = 'off';
end

if ~plotInApp
    close(h)
end

solution = improvedMFDfun(k,Ha,options.opt,options.optsLM,options.eigsopt,options.algROM,options.orderROM,[]);

return