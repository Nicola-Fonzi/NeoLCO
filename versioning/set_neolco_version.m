%*******************************************************************************
% Copyright (C) 2008 - 2021                                                    *
%                                                                              *
% Sergio Ricci (sergio.ricci@polimi.it)                                        *
%                                                                              *
% Politecnico di Milano, Dipartimento di Ingegneria Aerospaziale               *
% Via La Masa 34, 20156 Milano - ITALY                                         *
%                                                                              *
% This file is part of NeoCASS Software (www.neocass.org)                      *
%                                                                              *
% NeoCASS is free software; you can redistribute it and/or                     *
% modify it under the terms of the GNU General Public                          *
% License as published by the Free Software Foundation;                        *
% either version 3, or (at your option) any later version.                     *
%                                                                              *
% NeoCASS is distributed in the hope that it will be useful,                   *
% but WITHOUT ANY WARRANTY; without even the implied                           *
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR                      *
% PURPOSE.  See the GNU General Public License for more                        *
% details.                                                                     *
%                                                                              *
% You should have received a copy of the GNU General Public                    *
% License along with NeoCASS; see the file GNU GENERAL                         *
% PUBLIC LICENSE.TXT.  If not, see <http://www.gnu.org/licenses/>.             *
%*******************************************************************************
%                                                                              *
%                                                                              *
%                                                                              *
% Version: 3.0.0                                                               *
%                                                                              *
%                                                                              *
%                                                                              *
%*******************************************************************************
%                                                                              *
% Authors:                                                                     *
%                                                                              *
%                      Sergio Ricci         <sergio.ricci@polimi.it>           *
%                      Alessandro Degaspari <alessandro.degaspari@polimi.it>   *
%                      Luca Riccobene       <luca.riccobene@polimi.it>         *
%                      Federico Fonte       <federico.fonte@polimi.it>         *
%                      Francesco Toffol     <francesco.toffol@polimi.it>       *
%                      Nicola Fonzi         <nicola.fonzi@polimi.it>           *
%                                                                              *
%                                                                              *
% Politecnico di Milano, Dipartimento di Ingegneria Aerospaziale               *
% Via La Masa 34, 20156 Milano - ITALY                                         *
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

version_old = get_neocass_version();
versionNum_old = split(version_old,'.');
versionNum_old = cellfun(@(x) str2double(x),versionNum_old);

if versionNum_old(1) > versionNum(1) || ...
        versionNum_old(1) == versionNum(1) && versionNum_old(2) > versionNum(2) || ...
        versionNum_old(1) == versionNum(1) && versionNum_old(2) == versionNum(2) && versionNum_old(3) > versionNum(3)
    error('The old version is newer that the version you are trying to set')
end

thisFile = which('set_neocass_version');
[versionDir,~,~] = fileparts(thisFile);
mainDir = extractBetween(versionDir,'',strcat(filesep,'version'));
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
            if ~isfolder(files(iFile).name)
                content = fileread(files(iFile).name);
                % Remove previous header
                contentStripped = extractAfter(content,'function');
                fid = fopen(files(iFile).name,'w');
                if ~isempty(contentStripped)
                    % Write new header
                    printHeader(fid,version);
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

function printHeader(fid,version)
fprintf(fid,'%%*******************************************************************************\n');
fprintf(fid,'%% Copyright (C) 2020 - 2021                                                    *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% Nicola Fonzi (nicola.fonzi@polimi.it)                                        *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% Politecnico di Milano, Dipartimento di Ingegneria Aerospaziale               *\n');
fprintf(fid,'%% Via La Masa 34, 20156 Milano - ITALY                                         *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% This file is part of NeoLCO Software                                         *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% NeoLCO is free software; you can redistribute it and/or                      *\n');
fprintf(fid,'%% modify it under the terms of the GNU General Public                          *\n');
fprintf(fid,'%% License as published by the Free Software Foundation;                        *\n');
fprintf(fid,'%% either version 3, or (at your option) any later version.                     *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% NeoCASS is distributed in the hope that it will be useful,                   *\n');
fprintf(fid,'%% but WITHOUT ANY WARRANTY; without even the implied                           *\n');
fprintf(fid,'%% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR                      *\n');
fprintf(fid,'%% PURPOSE.  See the GNU General Public License for more                        *\n');
fprintf(fid,'%% details.                                                                     *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% You should have received a copy of the GNU General Public                    *\n');
fprintf(fid,'%% License along with NeoLCO; see the file GNU GENERAL                          *\n');
fprintf(fid,'%% PUBLIC LICENSE.TXT.  If not, see <http://www.gnu.org/licenses/>.             *\n');
fprintf(fid,'%%*******************************************************************************\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%% Version: %-68s*\n',version);
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%%                                                                              *\n');
fprintf(fid,'%%*******************************************************************************\n');
fprintf(fid,'function');
return
