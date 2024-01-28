clear
close all
clc
filename_sma = 'flutter.dat';
ver = get_neocass_version();
fprintf('\n\nNeoCASS version: %s\n\n', ver);
inputData = readSmartcadFile(filename_sma);
[model, options] = processInputData(inputData);
% Set options
fid = options.FID;
struOpt = [];
eigOpt = options.eig;
aeroOpt = options.aero;
flutterOptions.Mlist_dlm = options.aero.Mlist_dlm;
flutterOptions.klist = options.aero.klist;
flutterOptions.load = [];
flutterOptions.Vmax = 30;
flutterOptions.Vmin = 2;
flutterOptions.Vstep = 0.1;
flutterOptions.rho = 2.4; % V flutter linear 23.45
flutterOptions.method = 'PK0';
% Structural preprocessor
struData = structuralPreprocessor(fid, model, struOpt);
resultsEig = solve_eig_fun(fid, model, struData, eigOpt);
reducedBasis = defineReducedBasis(struData, resultsEig, 'all');
aeroData = [];
% Flutter computation
[resultsFlutter, aeroData] = solve_linflutt_fun(model, struData, reducedBasis, aeroData, flutterOptions);
plotVgDiagrams(resultsFlutter,[])
plotNeoModel(model)
