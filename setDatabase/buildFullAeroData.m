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
function [aeroData, options] = buildFullAeroData(model, model_stiff, struData_stiff, globalOptions, options)
%
% This function builds the aerodynamic database for all the subsequent
% analyses. The minimum required is the generation of the dlm matrices that
% are later used for both the flutter analyses and the generation of a
% state space model for aerodynamics. However, if a preload is to be
% included, also VLM matrices are built.
%
% Please note that in case a VLM correction of the DLM matrices is to be
% used, (i.e. dynVLM in NeoCASS), this function has to perform the trim
% analyses for all the conditions required, in order to save the necessary
% matrices. There is NO CHECK if the trim condition later used for the time
% marching simulation or the describing function simulation is the same
% used to create the database. IT IS ASSUMED that the user will create the
% database with the same trim condition later used for the solution. Also,
% in this case, only one stiffness combination should be used, as the trim
% would be different otherwise. This is corresponding to the configuration
% contained into model_stiff.
%
% Important note: the VLM correction to DLM matrices is velocity dependent,
% thus it should be recomputed per each flight condition. In practise, it
% is assumed that the trim condition specified is close enough to the
% conditions at which the system will operate and only one correction is
% applied.

% Set default options
iOpt = 0;
iOpt = iOpt+1; baseOpt.fidScreen = 1;               descr{iOpt} = 'fid for screen printing. [1].';
iOpt = iOpt+1; baseOpt.meshType = 'flat';           descr{iOpt}=  'Mesh type "fullflat" or "flat". ["flat"].';
iOpt = iOpt+1; baseOpt.noLinCurvature = false;      descr{iOpt} = 'Avoid the use of the linearized mesh curvature';
iOpt = iOpt+1; baseOpt.DynVLM = false;              descr{iOpt} = 'Use VLM correction to DLM matrices';
iOpt = iOpt+1; baseOpt.DynVLMtype = 'quasi-steady'; descr{iOpt} = 'Type of correction used if DynVLM is set to true: ''steady'' or ''quasi-steady'' or ''unsteady'' ';
iOpt = iOpt+1; baseOpt.selectionTrim = [];          descr{iOpt} = 'If more trim conditions are specified in the input file, this allows to select one. [].';

if nargin==0
    printOptionDescription(baseOpt, descr);
    return
end

% Process input options
options = setOptions(baseOpt, 'error', options);

options.Mlist_dlm = globalOptions.aero.Mlist_dlm;
options.klist = globalOptions.aero.klist;

aeroData = setAeroDatabase(globalOptions.FID, model, options);

if options.DynVLM
    % In this case we have to perform a trim solution so that we can
    % generate the required matrices later to be added to the DLM ones.
    
    if isempty(globalOptions.trim.ID)
        error('Requested the introduction of VLM correction to DLM, but no info about trim provided')
    end
    
    options.Mlist_vlm = globalOptions.trim.Mach;
    
    trimOptions.selectionList = options.selectionTrim;
    trimOptions.outputType = 'meanAxes';
    trimOptions.klist = globalOptions.aero.klist;
    trimOptions.DynVLM = true;
    trimOptions.DynVLMtype = options.DynVLMtype;
    
    [~, aeroData] = solve_lin_trim(globalOptions.FID, model_stiff, struData_stiff, aeroData, globalOptions.trim, trimOptions);

    aeroData = addVLMdynamicCorrection(globalOptions.FID, model_stiff, aeroData, globalOptions.trim, trimOptions);
    
end

return