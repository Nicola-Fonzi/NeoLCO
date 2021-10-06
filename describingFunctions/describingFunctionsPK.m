%*******************************************************************************
% Copyright (C) 2020 - 2021                                                    *
%                                                                              *
% Nicola Fonzi (nicola.fonzi@polimi.it)                                        *
%                                                                              *
% Politecnico di Milano, Dipartimento di Ingegneria Aerospaziale               *
% Via La Masa 34, 20156 Milano - ITALY                                         *
%                                                                              *
% This file is part of NeoLCO Software (github.com/Nicola-Fonzi/NeoLCO)        *
%                                                                              *
%*******************************************************************************
%                                                                              *
%                                                                              *
%                                                                              *
% Version: 2.0.0                                                               *
%                                                                              *
%                                                                              *
%                                                                              *
%*******************************************************************************
function [flutterSpeed, flutterFrequency, KeqVect, resultsFlutter, aeroData] = describingFunctionsPK(maxKeq, maxNumberKeq, inputData, HINGE_CELAS_ID, hingeScalarPoint, struOpt, fid, eigOpt, ...
    recomputeBase, flutterOptions, reducedBasis, searchFlutterQuenchingPoint, modesToPlotDF)

mkdir("DescribingFunctions")
chdir("DescribingFunctions")

aeroData = [];
KeqVect = [];
indexK = 1;
deltaK = maxKeq/maxNumberKeq;
Keq_prev = -deltaK;
reachedMinKeq = false;

while true
    
    Keq = Keq_prev+deltaK;
    
    %Note that the node 1 can be used only in inputData, in model we need
    %to provide the NeoCASS index
    inputData_descrFun = inputData;
    inputData_descrFun.CELAS.ID = [inputData.CELAS.ID HINGE_CELAS_ID];
    inputData_descrFun.CELAS.value = [inputData.CELAS.value Keq];
    inputData_descrFun.CELAS.Node = [inputData.CELAS.Node [hingeScalarPoint;0]];
    inputData_descrFun.CELAS.DOF = [inputData.CELAS.DOF [0;0]];
    inputData_descrFun.CELAS.PID = [inputData.CELAS.PID 0];
    
    [model_descrFun, ~] = processInputData(inputData_descrFun);
    
    % We need to rebuild the matrices as we changed the stiffness
    struData_descrFun = structuralPreprocessor(fid, model_descrFun, struOpt);
    
    % recompute the basis
    if recomputeBase
        eigOpt.UseFictmass = false;
        resultsEig_descrFun = solve_eig_fun(fid, model_descrFun, struData_descrFun, eigOpt);
        reducedBasis = defineReducedBasis(struData_descrFun, resultsEig_descrFun, 'all');
        aeroData = [];
    end
    
    % Flutter computation
    [resultsFlutter_descrFun, aeroData] = solve_linflutt_fun(model_descrFun, struData_descrFun, reducedBasis, aeroData, flutterOptions);
    
    if searchFlutterQuenchingPoint
        if any(any(real(resultsFlutter_descrFun.modes.eig)>=0)) || reachedMinKeq
            resultsFlutter.modes(indexK) = resultsFlutter_descrFun.modes;
            indexK = indexK+1;
            KeqVect = [KeqVect, Keq];
            Keq_prev = Keq;
            deltaK = maxKeq/maxNumberKeq;
            reachedMinKeq = false;
        else
            deltaK = deltaK/2;
            if deltaK<(maxKeq/maxNumberKeq/1000)
                reachedMinKeq = true;
                deltaK = maxKeq/maxNumberKeq;
            end
        end
    else
        resultsFlutter.modes(indexK) = resultsFlutter_descrFun.modes;
        indexK = indexK+1;
        KeqVect = [KeqVect, Keq];
        Keq_prev = Keq;
    end
    if indexK>maxNumberKeq
        break
    end
end

[flutterSpeed, flutterFrequency] = plotVgDiagrams(resultsFlutter, modesToPlotDF, false);
figure
plot(KeqVect,flutterSpeed,'LineWidth',2)
xlabel("Equivalen stiffness [Nm]")
ylabel("Flutter speed [m/s]")
figure
plot(KeqVect,flutterFrequency,'LineWidth',2)
xlabel("Equivalen stiffness [Nm]")
ylabel("Flutter frequency [Hz]")

handles=findall(0,'type','figure');

for i = 1:handles(1).Number
    h = figure(i);
    saveas(h,strcat("figure",num2str(i),".fig"))
    close(h)
end

clear h handles

cd("..")


return