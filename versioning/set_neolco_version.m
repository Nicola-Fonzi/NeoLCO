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
function set_neolco_version(version)
%**************************************************************************
% This file creates a header to all the source code files, and it also se a
% tag in the git history for the new version number.
% Input:
% A string strictly of this format x.x.xxx
%**************************************************************************

home = pwd;

versionNum = split(version,'.');
versionNum = cellfun(@(x) str2double(x),versionNum);

version_old = get_neolco_version();
versionNum_old = split(version_old,'.');
versionNum_old = cellfun(@(x) str2double(x),versionNum_old);

if versionNum_old(1) > versionNum(1) || ...
        versionNum_old(1) == versionNum(1) && versionNum_old(2) > versionNum(2) || ...
        versionNum_old(1) == versionNum(1) && versionNum_old(2) == versionNum(2) && versionNum_old(3) > versionNum(3)
    error('The old version is newer that the version you are trying to set')
end

thisFile = which('set_neolco_version');
[versionDir,~,~] = fileparts(thisFile);
mainDir = extractBetween(versionDir,'',strcat(filesep,'versioning'));
walkOnFiles(mainDir{1},version);

cd(home);

return

function walkOnFiles(mainDir,version)
subfolders = (genpath(mainDir));
if ispc()
    separator = ';';
else
    separator = ':';
end
subfolders = split(subfolders,separator);

for iSubfolder = 1:length(subfolders)
    if ~isempty(subfolders{iSubfolder})
        chdir(subfolders{iSubfolder})
        files = dir(pwd);
        for iFile = 1:length(files)
            if ~isfolder(files(iFile).name) && my_isfunction(files(iFile).name)==1 && ~contains(pwd,".git")
                content = fileread(files(iFile).name);
                % Remove previous header
                contentStripped = extractAfter(content,'function');
                fid = fopen(files(iFile).name,'w');
                if ~isempty(contentStripped)
                    % Write new header
                    if strcmp(files(iFile).name,"get_neolco_version.m")
                        printHeader(fid,version,true);
                    else
                        printHeader(fid,version);
                    end
                    % Write code
                    fwrite(fid,contentStripped);
                else
                    % Write code
                    fwrite(fid,content);
                end
                fclose(fid);
            end
        end
    end
end

return

function printHeader(fid,version,longerHeader)
fprintf(fid,'%%*******************************************************************************\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%%                    _   _            _     ____ ___                           *\n');
fprintf(fid,'%%                   | \\ | | ___  ___ | |   / ___/ _ \\                          *\n');
fprintf(fid,'%%                   |  \\| |/ _ \\/ _ \\| |  | |  | | | |                         *\n');
fprintf(fid,'%%                   | |\\  |  __/ (_) | |__| |__| |_| |                         *\n');
fprintf(fid,'%%                   |_| \\_|\\___|\\___/|_____\\____\\___/                          *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% Copyright (C) 2020 - 2024                                                    *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% Nicola Fonzi (nicola.fonzi@outlook.com)                                      *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% This file is part of NeoLCO Software (github.com/Nicola-Fonzi/NeoLCO).       *\n');
fprintf(fid,'%% The use of this software is licensed based on the licence distributed        *\n');
fprintf(fid,'%% together with the source code. If you have not received the license please   *\n');
fprintf(fid,'%% contact the copywright owner before using the software.                      *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%%*******************************************************************************\n');
if nargin == 3 && longerHeader
    fprintf(fid,'%%                                                                              *\n');
    fprintf(fid,'%%                                                                              *\n');
    fprintf(fid,'%%                                                                              *\n');
    fprintf(fid,'%% Version: %-68s*\n',version);
    fprintf(fid,'%%                                                                              *\n');
    fprintf(fid,'%%                                                                              *\n');
    fprintf(fid,'%%                                                                              *\n');
    fprintf(fid,'%%*******************************************************************************\n');
end
fprintf(fid,'function');
return

function ID = my_isfunction(FUNNAME)
[filepath,name,ext] = fileparts(FUNNAME);
if strcmp(ext, '.mlapp')
    ID = 0;
    return
end
try    
    nargin(FUNNAME) ; % nargin errors when FUNNAME is not a function
    ID = 1  + isa(FUNNAME, 'function_handle') ; % 1 for m-file, 2 for handle
catch ME
    % catch the error of nargin
    switch (ME.identifier)        
        case 'MATLAB:nargin:isScript'
            ID = -1 ; % script
        case 'MATLAB:narginout:notValidMfile'
            ID = -2 ; % probably another type of file, or it does not exist
        case 'MATLAB:narginout:functionDoesnotExist'
            ID = -3 ; % probably a handle, but not to a function
        case 'MATLAB:narginout:BadInput'
            ID = -4 ; % probably a variable or an array
        otherwise
            ID = 0 ; % unknown cause for error
    end
end
return