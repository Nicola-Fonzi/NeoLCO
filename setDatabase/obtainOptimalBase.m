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
% You are not authorized to use, distribute, or modify this file in any way,   *
% unless explicitly decleared otherwise by the copyright owner.                *
%                                                                              *
%*******************************************************************************
function [reducedBasis, struData, model, struData_stiff, model_stiff, globalOptions] = obtainOptimalBase(inputData, options)
% This function always assume that the model is provided to NeoLCO so that
% the stiffness at the nonlinearity points is NULL. It will be later added.
% Also, no other options should be specified for the nonlinearity points,
% like SUPORT or AESURFD cards.
%
% The way this function works is by taking the extrema in stiffness values
% at the nonlinearity points and compute the base which correctly reproduce
% the frequency response of the stiff and free system.
%

% Set default options
iOpt = 0;
iOpt = iOpt+1; baseOpt.fidScreen = 1;                 descr{iOpt} = 'fid for screen printing. [1].';
iOpt = iOpt+1; baseOpt.kNominal = {};                 descr{iOpt} = 'Cells array containing the nominal stiffnesses at the nonlinearity points. The rows are in the same order as gapPoint IDs. The columns contain possible different values for the same point. The format is {[k1_point1,k2_point1];[k1_point2,k2_point2,k3_point2]}. {}.';
iOpt = iOpt+1; baseOpt.gapPoints = {};                descr{iOpt} = 'Points where nonlinearities are present. The format is {point1,"s" or "g",dof (1,2,3,4,5 or 6),label;point2,...}. If in the same point we have more nonlinearity, the point must be repeated per each dof. {}.';
iOpt = iOpt+1; baseOpt.referenceNModes = 100;         descr{iOpt} = 'Number of modes to be retained and used as a reference for the "correct" system. [100].';
iOpt = iOpt+1; baseOpt.optimumKnown = 0;              descr{iOpt} = 'Option that can be used in case the optimum base is known. 0 for false, 1 for true. [0]';
iOpt = iOpt+1; baseOpt.nModes = [];                   descr{iOpt} = 'Vector of modes to be tested. []';
iOpt = iOpt+1; baseOpt.stiffnessFactor = [];          descr{iOpt} = 'Vector of factors that multiplies the nominal stiffness to be tested. [].';
iOpt = iOpt+1; baseOpt.fictitiousMass = [];           descr{iOpt} = 'Vector of fictitious masses to be tested. []';
iOpt = iOpt+1; baseOpt.maxFreqInterest = 100;         descr{iOpt} = 'Max frequency that is of interest. Inside this band the error will be computed. [100 Hz]';

if nargin==0
    printOptionDescription(baseOpt, descr);
    return
end

% Process input options
options = setOptions(baseOpt, 'error', options);

[model, globalOptions] = processInputData(inputData);

maximumStiffnesses = cellfun(@(x) max(x), options.kNominal);

[model_stiff, ~] = addNonlinearityStiffness(model, options.gapPoints, maximumStiffnesses);

%% Gather missing options from the model
struOpt = [];

eigOpt = globalOptions.eig;

fid = globalOptions.FID;

%% Preprocessing

struData = structuralPreprocessor(fid, model, struOpt);

struData_stiff = structuralPreprocessor(fid, model_stiff, struOpt);

%% Recover of the DOF indeces required for postprocessing

nonlinearityDOF = obtainDOF(options.gapPoints, model);


%% Create the nominal bases

eigOpt.NROOTS = options.referenceNModes;

eigOpt.UseFictmass = false;
resultsEig = solve_eig_fun(fid, model, struData, eigOpt);

% We extract the modes at zero frequency that are not supported, as these
% contain the nonlinearity rigid modes
rigidModes.V = resultsEig.V(:,resultsEig.Omega<1);
rigidModes.suportTable = resultsEig.suportTable;
rigidModes.surfNames = resultsEig.surfNames;

eigOpt.UseFictmass = false;
resultsEig_stiff = solve_eig_fun(fid, model_stiff, struData_stiff, eigOpt);


%% Loop over all the possibilities and test them one by one

home = pwd;

if ~options.optimumKnown
    
    mkdir("comparisonDifferentBases")
    chdir("comparisonDifferentBases")
    homeDir = pwd;
    
    frequencyResponseErrorStiff = nan(length(options.nModes),length(options.stiffnessFactor),length(options.fictitiousMass));
    frequencyResponseErrorFree = nan(length(options.nModes),length(options.stiffnessFactor),length(options.fictitiousMass));
    frequencyResponseErrorStiff_RB = nan(length(options.nModes),length(options.stiffnessFactor),length(options.fictitiousMass));
    frequencyResponseErrorFree_RB = nan(length(options.nModes),length(options.stiffnessFactor),length(options.fictitiousMass));
    
    for i = options.nModes
        mkdir(strcat("nModes",num2str(i)))
        chdir(strcat("nModes",num2str(i)))
        homeModes = pwd;
        
        eigOpt.NROOTS = i;
        
        for j = options.stiffnessFactor
            
            mkdir(strcat("stiffnessFactor",num2str(j)))
            chdir(strcat("stiffnessFactor",num2str(j)))
            homeStiffness = pwd;
            
            for k = options.fictitiousMass
                
                mkdir(strcat("fictitiousMass",num2str(k)))
                chdir(strcat("fictitiousMass",num2str(k)))
                
                [model_common, ~] = addNonlinearityStiffness(model, options.gapPoints, maximumStiffnesses*j);
                [model_common, ~] = addFictitiousMass(model_common, options.gapPoints, repmat(k,size(options.gapPoints,1),1));
                
                struData_common = structuralPreprocessor(fid, model_common, struOpt);
                
                eigOpt.UseFictmass = true;
                resultsEig_common = solve_eig_fun(fid, model_common, struData_common, eigOpt);
                
                [D_free, V_free, D_stiff, V_stiff] = recomputeEigResults(struData, struData_stiff, resultsEig_common, j);
                
                % We add the rigid body modes of the free system, as they
                % include the static modes.
                % The orthogonality is set with respect to the mass ---> same mass matrix for free and stiff
                reducedBasis_common = defineReducedBasis(struData, rigidModes, 'all', resultsEig_common, 'all');
                
                [frequencyResponseErrorFree, frequencyResponseErrorStiff, frequencyResponseErrorFree_RB, frequencyResponseErrorStiff_RB] ...
                    = plotErrors(resultsEig, resultsEig_stiff, resultsEig_common, reducedBasis_common, D_free, D_stiff, V_free, ...
                    V_stiff, struData, struData_stiff, k, i, j, nonlinearityDOF, options,...
                    frequencyResponseErrorFree, frequencyResponseErrorStiff, frequencyResponseErrorFree_RB, frequencyResponseErrorStiff_RB);
                
                chdir(homeStiffness)
            end
            chdir(homeModes)
        end
        chdir(homeDir)
    end
    
    plotSurfError(frequencyResponseErrorStiff,options,"Stiff")
    plotSurfError(frequencyResponseErrorFree,options,"Free")
    plotSurfError(frequencyResponseErrorFree_RB,options,"FreeRB")
    plotSurfError(frequencyResponseErrorStiff_RB,options,"StiffRB")
    
    chdir(home)
    
end

SetValues = true;

while SetValues
    
    fictitiousMassRequired = input("Chosen value for the fictitious mass: ");
    stiffnessFactorRequired = input("Chosen value for the stiffness factor value: ");
    nModesRequired = input("Chosen number of modes: ");
    
    SetValues = input(strcat("You requested fictitious mass = ",num2str(fictitiousMassRequired),"; stiffness factor = ",num2str(stiffnessFactorRequired),"; number modes = ",num2str(nModesRequired),"... Ok? [0=no, 1=yes]")) == 0;
    
end

eigOpt.NROOTS = nModesRequired;

[model_common, ~] = addNonlinearityStiffness(model, options.gapPoints, maximumStiffnesses*stiffnessFactorRequired);
[model_common, ~] = addFictitiousMass(model_common, options.gapPoints, repmat(fictitiousMassRequired,size(options.gapPoints,1),1));

struData_common = structuralPreprocessor(fid, model_common, struOpt);

eigOpt.UseFictmass = true;
resultsEig_common = solve_eig_fun(fid, model_common, struData_common, eigOpt);

reducedBasis = defineReducedBasis(struData, rigidModes, 'all', resultsEig_common, 'all');

if ~isempty(model.Tabdmp1.ID(model.Tabdmp1.ID==model.Tabdmp1.caseControl))
    KDAMP = model.Tabdmp1.KDAMP;
    if KDAMP < 0
        error("Only real damping is supported in this context, please use KDAMP>0")
    else
        fprintf(fid,' - Building real damping matrix for viscous damping...');
        reducedBasis.Bmm = modalDamp(model, reducedBasis.V'*struData.Mzz*reducedBasis.V, reducedBasis.V'*struData.Kzz*reducedBasis.V,'Real');
    end
    fprintf(fid, 'done.\n');
end

resultsFolder = strcat("FM",num2str(fictitiousMassRequired),"KF",num2str(stiffnessFactorRequired),"NMODES",num2str(nModesRequired));

mkdir(resultsFolder)
chdir(resultsFolder)

close all

return

function [D_free, V_free, D_stiff, V_stiff] = recomputeEigResults(struData, struData_stiff, resultsEig_common, stiffnessFactor)

warning off;
[V, D] = eig(resultsEig_common.V'*full(struData_stiff.Kzz)*resultsEig_common.V,...
    resultsEig_common.V'*full(struData.Mzz)*resultsEig_common.V);
warning on;
eigenvalue = diag(D);
[eigenvalue, indexROTT] = sort(eigenvalue);
V = V(:,indexROTT);
D_stiff = sqrt(abs(eigenvalue))/2/pi;
V_stiff = resultsEig_common.V*V;

if stiffnessFactor
    [V, D] = eig(resultsEig_common.V'*full(struData.Kzz)*resultsEig_common.V,...
        resultsEig_common.V'*full(struData.Mzz)*resultsEig_common.V);
    warning on;
    eigenvalue = diag(D);
    [eigenvalue, indexROTT] = sort(eigenvalue);
    V = V(:,indexROTT);
    D_free = sqrt(abs(eigenvalue))/2/pi;
    V_free = resultsEig_common.V*V;
else
    V_free = resultsEig_common.V;
    D_free = resultsEig_common.Omega/2/pi;
end

return


function [frequencyResponseErrorFree, frequencyResponseErrorStiff, frequencyResponseErrorFree_RB, frequencyResponseErrorStiff_RB]...
    = plotErrors(resultsEig, resultsEig_stiff, resultsEig_common, reducedBasis_common, D_free, D_stiff, V_free, ...
    V_stiff, struData, struData_stiff, fictitiousMass, nModes, stiffnessFactor, nonlinearityDOF, options, ...
    frequencyResponseErrorFree, frequencyResponseErrorStiff, frequencyResponseErrorFree_RB, frequencyResponseErrorStiff_RB)

figure
subplot(4,2,[1 2])
semilogy(abs(resultsEig.Omega(1:length(D_free))/2/pi - D_free)./(resultsEig.Omega(1:length(D_free))/2/pi)*100,'LineWidth',2)
hold on
semilogy(abs(resultsEig_stiff.Omega(1:length(D_stiff))/2/pi - D_stiff)./(resultsEig_stiff.Omega(1:length(D_stiff))/2/pi)*100,'LineWidth',2)
grid on
xlabel("Mode index")
ylabel("Error frequency [%]")
legend("Free","Stiff")

subplot(4,2,3)
mac_free = computeModalMAC(resultsEig.V, V_free);
semilogy(abs(diag(mac_free)))
xlabel("Mode index")
ylabel("MAC diagonal element")
title("Free")

subplot(4,2,4)
mac_stiff = computeModalMAC(resultsEig_stiff.V, V_stiff);
semilogy(abs(diag(mac_stiff)))
xlabel("Mode index")
ylabel("MAC diagonal")
title("Stiff")

wmax = 2*pi*options.maxFreqInterest;
w = linspace(5,wmax,100);

subplot(4,2,5)

frf_free_common = computeFRF(w,struData,resultsEig_common,nonlinearityDOF);
frf_free = computeFRF(w,struData,resultsEig,nonlinearityDOF);
frfError_free = computeFRFerror(frf_free_common, frf_free, "full");
semilogy(w/2/pi,frfError_free,'LineWidth',2)
xlabel("Frequency [Hz]")
ylabel("FRF error [%]")
title("Free")

frequencyResponseErrorFree(nModes==options.nModes,stiffnessFactor==options.stiffnessFactor,fictitiousMass==options.fictitiousMass) = computeFRFerror(frf_free_common,frf_free,"condensed");

subplot(4,2,6)

frf_stiff_common = computeFRF(w,struData_stiff,resultsEig_common,nonlinearityDOF);
frf_stiff = computeFRF(w,struData_stiff,resultsEig_stiff,nonlinearityDOF);
frfError_stiff = computeFRFerror(frf_stiff_common, frf_stiff, "full");
semilogy(w/2/pi,frfError_stiff,'LineWidth',2)
xlabel("Frequency [Hz]")
ylabel("FRF error [%]")
title("Stiff")

frequencyResponseErrorStiff(nModes==options.nModes,stiffnessFactor==options.stiffnessFactor,fictitiousMass==options.fictitiousMass) = computeFRFerror(frf_stiff_common,frf_stiff,"condensed");

subplot(4,2,7)
frf_free_common_RB = computeFRF(w,struData,reducedBasis_common,nonlinearityDOF);
frfError_free_RB = computeFRFerror(frf_free_common_RB, frf_free, "full");
semilogy(w/2/pi,frfError_free_RB,'LineWidth',2)
xlabel("Frequency [Hz]")
ylabel("FRF error [%]")
title("Free")

frequencyResponseErrorFree_RB(nModes==options.nModes,stiffnessFactor==options.stiffnessFactor,fictitiousMass==options.fictitiousMass) = computeFRFerror(frf_free_common_RB,frf_free,"condensed");

subplot(4,2,8)
frf_stiff_common_RB = computeFRF(w, struData_stiff, reducedBasis_common, nonlinearityDOF);
frfError_stiff_RB = computeFRFerror(frf_stiff_common_RB, frf_stiff, "full");
semilogy(w/2/pi,frfError_stiff_RB,'LineWidth',2)
xlabel("Frequency [Hz]")
ylabel("FRF error [%]")
title("Stiff")

frequencyResponseErrorStiff_RB(nModes==options.nModes,stiffnessFactor==options.stiffnessFactor,fictitiousMass==options.fictitiousMass) = computeFRFerror(frf_stiff_common_RB,frf_stiff,"condensed");

sgtitle(strcat("Fictitious mass = ",num2str(fictitiousMass),"; Number modes = ",num2str(nModes),"; stiffness Factor = ",num2str(stiffnessFactor)))

saveas(gcf,"MAC.fig")

return

function FRF = computeFRF(w,struData,reducedBasis,DOF)
FRF = zeros(length(DOF),length(DOF),length(w));
for i=1:length(w)
    FRF(:,:,i) = struData.Tgz(DOF,:)*reducedBasis.V*((-reducedBasis.V'*struData.Mzz*reducedBasis.V*w(i)^2 ...
        + reducedBasis.V'*struData.Kzz*reducedBasis.V)\(struData.Tgz(DOF,:)*reducedBasis.V)');
end
return

function err = computeFRFerror(frf,frfReference,type)

if size(frf) ~= size(frfReference)
    error("Wrong size for the two frequency response matrices")
end

err = nan(size(frf,3),1);
for w = 1:size(frf,3)
    err(w) = sqrt(sum(sum(((frf(:,:,w)-frfReference(:,:,w))./frfReference(:,:,w)*100).^2/(size(frf,1)*size(frf,2)))));
end

if strcmp(type,"condensed")
    err = sum(err)/sqrt(size(frf,3));
end

return

function plotSurfError(data,options,label)

for imodes = 1:length(options.nModes)
    
    figure
    surf(log10(options.fictitiousMass),log10(options.stiffnessFactor),log10(squeeze(data(imodes,:,:))),'LineWidth',2)
    xlabel("Fictitious mass")
    ylabel("Stiffness factor")
    title(strcat(label,"; Number modes = ",num2str(options.nModes(imodes))));
    theCurrentAxis = gca;
    for i = 1:length(theCurrentAxis.XTickLabel)
        Scaled(i)=str2double(theCurrentAxis.XTickLabel{i});
        Label{i} = num2str(10^Scaled(i));
    end
    set(gca,'XTickLabel',Label);
    clear Scaled Label
    for i = 1:length(theCurrentAxis.YTickLabel)
        Scaled(i)=str2double(theCurrentAxis.YTickLabel{i});
        Label{i} = num2str(10^Scaled(i));
    end
    set(gca,'YTickLabel',Label);
    clear Scaled Label
    for i = 1:length(theCurrentAxis.ZTickLabel)
        Scaled(i)=str2double(theCurrentAxis.ZTickLabel{i});
        Label{i} = num2str(10^Scaled(i));
    end
    set(gca,'ZTickLabel',Label);
    clear Scaled Label
    
    saveas(gcf,strcat(label,"nModes=",num2str(options.nModes(imodes)),".fig"))
    
end

return