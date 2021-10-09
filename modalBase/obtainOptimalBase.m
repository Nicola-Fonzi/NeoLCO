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
function [reducedBasis, resultsFolder] = obtainOptimalBase(inputData, model, struData, options, model_stiff, struData_stiff, parameters, HingeDOF)


%% Set the missing options

struOpt = [];

eigOpt = options.eig;

fid = options.FID;


%% Create the nominal bases

eigOpt.NROOTS = 80; % Large number of roots as the "correct" system

eigOpt.UseFictmass = false;
resultsEig = solve_eig_fun(fid, model, struData, eigOpt);

eigOpt.UseFictmass = false;
resultsEig_stiff = solve_eig_fun(fid, model_stiff, struData_stiff, eigOpt);


%% Loop over all the possibilities and test them one by one

if parameters.OptimumKnown==0
    
    mkdir("Comparison_different_bases")
    chdir("Comparison_different_bases")
    
    frequencyResponseErrorStiff = nan(length(parameters.NMODES),length(parameters.KN),length(parameters.FM));
    frequencyResponseErrorFree = nan(length(parameters.NMODES),length(parameters.KN),length(parameters.FM));
    
    wmax = inf;
    
    for i = parameters.NMODES
        mkdir(strcat("NMODES",num2str(i)))
        chdir(strcat("NMODES",num2str(i)))
        
        eigOpt.NROOTS = i;
        
        for j = parameters.KN
            
            mkdir(strcat("KN",num2str(j)))
            chdir(strcat("KN",num2str(j)))
            
            for k = parameters.FM
                
                wmax=inf;
                
                mkdir(strcat("FM",num2str(k)))
                chdir(strcat("FM",num2str(k)))
                
                
                inputData_common = inputData;
                inputData_common.CMASS.value(inputData_common.CMASS.ID==parameters.FICTMASS_CMASS_ID) = k;
                
                if j
                    inputData_common.CELAS.ID = [inputData_common.CELAS.ID parameters.HINGE_CELAS_ID];
                    inputData_common.CELAS.value = [inputData_common.CELAS.value j];
                    inputData_common.CELAS.Node = [inputData_common.CELAS.Node [parameters.HINGE_POINT_ID;0]];
                    inputData_common.CELAS.DOF = [inputData_common.CELAS.DOF [0;0]];
                    inputData_common.CELAS.PID = [inputData_common.CELAS.PID 0];
                end
                
                [model_common, ~] = processInputData(inputData_common);
                
                struData_common = structuralPreprocessor(fid, model_common, struOpt);
                
                eigOpt.UseFictmass = true;
                resultsEig_common = solve_eig_fun(fid, model_common, struData_common, eigOpt);
                
                warning off;
                [V, D] = eig(resultsEig_common.V'*full(struData_stiff.Kzz)*resultsEig_common.V,...
                    resultsEig_common.V'*full(struData.Mzz)*resultsEig_common.V);
                warning on;
                eigenvalue = diag(D);
                [eigenvalue, indexROTT] = sort(eigenvalue);
                V = V(:,indexROTT);
                D_stiff = sqrt(abs(eigenvalue))/2/pi;
                V_stiff = resultsEig_common.V*V;
                
                if j
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
                
                V_struct.V = [resultsEig.V(:,1) resultsEig_common.V];
                V_struct.surfNames = resultsEig.surfNames;
                V_struct.suportTable = resultsEig.suportTable;
                reducedBasis_common = defineReducedBasis(struData, V_struct, 'all'); % The orthogonality is set with respect to the mass ---> same for free and stiff
                
                figure
                subplot(4,2,1)
                semilogy(abs(resultsEig.Omega(1:length(D_free))/2/pi - D_free),'LineWidth',2)
                hold on
                semilogy(abs(resultsEig_stiff.Omega(1:length(D_free))/2/pi - D_stiff),'LineWidth',2)
                grid on
                xlabel("Mode index")
                ylabel("Delta natural frequency [Hz]")
                legend("Free","Stiff")
                
                
                subplot(4,2,2)
                title(strcat("FM = ",num2str(k),"; NMODES = ",num2str(i),"; KN = ",num2str(j)))
                
                subplot(4,2,3)
                for index = 1:length(D_free)
                    for jndex = 1:length(D_free)
                        mac_free(index,jndex) = (resultsEig.V(:,index).'*V_free(:,jndex))^2/norm(V_free(:,jndex))^2/norm(resultsEig.V(:,index))^2;
                    end
                end
                semilogy(abs(diag(mac_free)))
                xlabel("Mode index")
                ylabel("MAC")
                title("Free")
                
                subplot(4,2,4)
                for index = 1:length(D_free)
                    for jndex = 1:length(D_free)
                        mac_stiff(index,jndex) = (resultsEig_stiff.V(:,index).'*V_stiff(:,jndex))^2/norm(V_stiff(:,jndex))^2/norm(resultsEig_stiff.V(:,index))^2;
                    end
                end
                semilogy(abs(diag(mac_stiff)))
                xlabel("Mode index")
                ylabel("MAC")
                title("Stiff")
                
                subplot(4,2,5)
                w = linspace(5,(D_free(end)-10)*2*pi,1000);
                if max(w)<wmax
                    wmax = max(w);
                end
                
                for index=1:length(w)
                    frf_free_common(index) = struData.Tgz(HingeDOF,:)*resultsEig_common.V*((-resultsEig_common.V'*struData.Mzz*resultsEig_common.V*w(index)^2 ...
                        + resultsEig_common.V'*struData.Kzz*resultsEig_common.V)\(struData.Tgz(HingeDOF,:)*resultsEig_common.V)');
                end
                for index=1:length(w)
                    frf_free(index) = struData.Tgz(HingeDOF,:)*resultsEig.V*((-resultsEig.V'*struData.Mzz*resultsEig.V*w(index)^2 + ...
                        resultsEig.V'*struData.Kzz*resultsEig.V)\(struData.Tgz(HingeDOF,:)*resultsEig.V)');
                end
                semilogy(w/2/pi,abs(frf_free_common-frf_free),'LineWidth',2)
                title("FRF difference for free case")
                xlabel("Frequency [Hz]")
                ylabel("FRF")
                
                frequencyResponseErrorFree(i==parameters.NMODES,j==parameters.KN,k==parameters.FM) = norm(frf_free_common(1:find(w<=wmax,1,'last'))-frf_free(1:find(w<=wmax,1,'last')))/norm(frf_free(1:find(w<=wmax,1,'last')));
                
                subplot(4,2,6)
                for index=1:length(w)
                    frf_stiff_common(index) = struData_stiff.Tgz(HingeDOF,:)*resultsEig_common.V*((-resultsEig_common.V'*struData_stiff.Mzz*resultsEig_common.V*w(index)^2 ...
                        + resultsEig_common.V'*struData_stiff.Kzz*resultsEig_common.V)\(struData_stiff.Tgz(HingeDOF,:)*resultsEig_common.V)');
                end
                for index=1:length(w)
                    frf_stiff(index) = struData_stiff.Tgz(HingeDOF,:)*resultsEig_stiff.V*((-resultsEig_stiff.V'*struData_stiff.Mzz*resultsEig_stiff.V*w(index)^2 + ...
                        resultsEig_stiff.V'*struData_stiff.Kzz*resultsEig_stiff.V)\(struData_stiff.Tgz(HingeDOF,:)*resultsEig_stiff.V)');
                end
                semilogy(w/2/pi,abs(frf_stiff_common-frf_stiff),'LineWidth',2)
                title("FRF difference for stiff case")
                xlabel("Frequency [Hz]")
                ylabel("FRF")
                
                frequencyResponseErrorStiff(i==parameters.NMODES,j==parameters.KN,k==parameters.FM) = norm(frf_stiff_common(1:find(w<=wmax,1,'last'))-frf_stiff(1:find(w<=wmax,1,'last')))/norm(frf_stiff(1:find(w<=wmax,1,'last')));
                
                subplot(4,2,7)
                for index=1:length(w)
                    frf_free_common_RB(index) = struData.Tgz(HingeDOF,:)*reducedBasis_common.V*((-reducedBasis_common.V'*struData.Mzz*reducedBasis_common.V*w(index)^2 ...
                        + reducedBasis_common.V'*struData.Kzz*reducedBasis_common.V)\(struData.Tgz(HingeDOF,:)*reducedBasis_common.V)');
                end
                semilogy(w/2/pi,abs(frf_free_common_RB-frf_free),'LineWidth',2)
                title("FRF difference for free case")
                xlabel("Frequency [Hz]")
                ylabel("FRF")
                
                frequencyResponseErrorFree_RB(i==parameters.NMODES,j==parameters.KN,k==parameters.FM) = norm(frf_free_common_RB(1:find(w<=wmax,1,'last'))-frf_free(1:find(w<=wmax,1,'last')))/norm(frf_free(1:find(w<=wmax,1,'last')));
                
                subplot(4,2,8)
                for index=1:length(w)
                    frf_stiff_common_RB(index) = struData_stiff.Tgz(HingeDOF,:)*reducedBasis_common.V*((-reducedBasis_common.V'*struData_stiff.Mzz*reducedBasis_common.V*w(index)^2 ...
                        + reducedBasis_common.V'*struData_stiff.Kzz*reducedBasis_common.V)\(struData_stiff.Tgz(HingeDOF,:)*reducedBasis_common.V)');
                end
                semilogy(w/2/pi,abs(frf_stiff_common_RB-frf_stiff),'LineWidth',2)
                title("FRF difference for stiff case")
                xlabel("Frequency [Hz]")
                ylabel("FRF")
                
                frequencyResponseErrorStiff_RB(i==parameters.NMODES,j==parameters.KN,k==parameters.FM) = norm(frf_stiff_common_RB(1:find(w<=wmax,1,'last'))-frf_stiff(1:find(w<=wmax,1,'last')))/norm(frf_stiff(1:find(w<=wmax,1,'last')));
                
                saveas(gcf,"MAC.fig")
                
                chdir("..")
            end
            chdir("..")
        end
        chdir("..")
    end
    
    for index = 1:length(parameters.NMODES)
        figure
        surf(log10(parameters.FM),log10(parameters.KN),log10(squeeze(frequencyResponseErrorStiff(index,:,:))),'LineWidth',2)
        xlabel("FM")
        ylabel("KN")
        title(strcat("Stiff","; NMODES = ",num2str(parameters.NMODES(index))));
        theCurrentAxis = gca;
        for jndex = 1:length(theCurrentAxis.XTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.XTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'XTickLabel',Label);
        clear Scaled Label
        for jndex = 1:length(theCurrentAxis.YTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.YTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'YTickLabel',Label);
        clear Scaled Label
        for jndex = 1:length(theCurrentAxis.ZTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.ZTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'ZTickLabel',Label);
        clear Scaled Label
        
        saveas(gcf,strcat("Stiff","NMODES=",num2str(parameters.NMODES(index)),".fig"))
        
    end
    
    for index = 1:length(parameters.NMODES)
        figure
        surf(log10(parameters.FM),log10(parameters.KN),log(squeeze(frequencyResponseErrorFree(index,:,:))),'LineWidth',2)
        xlabel("FM")
        ylabel("KN")
        title(strcat("Free","; NMODES = ",num2str(parameters.NMODES(index))));
        theCurrentAxis = gca;
        for jndex = 1:length(theCurrentAxis.XTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.XTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'XTickLabel',Label);
        clear Scaled Label
        for jndex = 1:length(theCurrentAxis.YTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.YTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'YTickLabel',Label);
        clear Scaled Label
        for jndex = 1:length(theCurrentAxis.ZTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.ZTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'ZTickLabel',Label);
        clear Scaled Label
        
        saveas(gcf,strcat("Free","NMODES=",num2str(parameters.NMODES(index)),".fig"))
        
    end
    
    
    for index = 1:length(parameters.NMODES)
        figure
        surf(log10(parameters.FM),log10(parameters.KN),log(squeeze(frequencyResponseErrorFree_RB(index,:,:))),'LineWidth',2)
        xlabel("FM")
        ylabel("KN")
        title(strcat("Free RB","; NMODES = ",num2str(parameters.NMODES(index))));
        theCurrentAxis = gca;
        for jndex = 1:length(theCurrentAxis.XTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.XTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'XTickLabel',Label);
        clear Scaled Label
        for jndex = 1:length(theCurrentAxis.YTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.YTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'YTickLabel',Label);
        clear Scaled Label
        for jndex = 1:length(theCurrentAxis.ZTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.ZTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'ZTickLabel',Label);
        clear Scaled Label
        
        saveas(gcf,strcat("FreeRB","NMODES=",num2str(parameters.NMODES(index)),".fig"))
        
    end
    
    for index = 1:length(parameters.NMODES)
        figure
        surf(log10(parameters.FM),log10(parameters.KN),log(squeeze(frequencyResponseErrorStiff_RB(index,:,:))),'LineWidth',2)
        xlabel("FM")
        ylabel("KN")
        title(strcat("Stiff RB","; NMODES = ",num2str(parameters.NMODES(index))));
        theCurrentAxis = gca;
        for jndex = 1:length(theCurrentAxis.XTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.XTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'XTickLabel',Label);
        clear Scaled Label
        for jndex = 1:length(theCurrentAxis.YTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.YTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'YTickLabel',Label);
        clear Scaled Label
        for jndex = 1:length(theCurrentAxis.ZTickLabel)
            Scaled(jndex)=str2double(theCurrentAxis.ZTickLabel{jndex});
            Label{jndex} = num2str(10^Scaled(jndex));
        end
        set(gca,'ZTickLabel',Label);
        clear Scaled Label
        
        saveas(gcf,strcat("StiffRB","NMODES=",num2str(parameters.NMODES(index)),".fig"))
        
    end
    
    chdir("..")
    
end

SetValues = true;

while SetValues
    
    FM_Required = input("Chosen value for the fictitious mass: ");
    KN_Required = input("Chosen value for the stiffness value: ");
    NMODES_Required = input("Chosen number of modes: ");
    
    SetValues = input(strcat("You requested FM = ",num2str(FM_Required),"; KN = ",num2str(KN_Required),"; NMODES = ",num2str(NMODES_Required),"... Ok?")) == 0;
    
end

eigOpt.NROOTS = NMODES_Required;

inputData_common = inputData;
inputData_common.CMASS.value(inputData_common.CMASS.ID==parameters.FICTMASS_CMASS_ID) = FM_Required;

if KN_Required
    inputData_common.CELAS.ID = [inputData_common.CELAS.ID parameters.HINGE_CELAS_ID];
    inputData_common.CELAS.value = [inputData_common.CELAS.value KN_Required];
    inputData_common.CELAS.Node = [inputData_common.CELAS.Node [1;0]];
    inputData_common.CELAS.DOF = [inputData_common.CELAS.DOF [0;0]];
    inputData_common.CELAS.PID = [inputData_common.CELAS.PID 0];
end

[model_common, ~] = processInputData(inputData_common);

struData_common = structuralPreprocessor(fid, model_common, struOpt);

eigOpt.UseFictmass = true;
resultsEig_common = solve_eig_fun(fid, model_common, struData_common, eigOpt);

V_struct.V = [resultsEig.V(:,1) resultsEig_common.V];
V_struct.surfNames = resultsEig.surfNames;
V_struct.suportTable = resultsEig.suportTable;
V_struct.SOL = "Reduced basis";

reducedBasis = defineReducedBasis(struData, V_struct, 'all');

resultsFolder = strcat("FM",num2str(FM_Required),"KN",num2str(KN_Required),"NMODES",num2str(NMODES_Required));

close all

return