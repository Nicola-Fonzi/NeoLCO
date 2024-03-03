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
%                                                                              *
%                                                                              *
%                                                                              *
% Version: 3.0.0                                                               *
%                                                                              *
%                                                                              *
%                                                                              *
%*******************************************************************************
function version = get_neolco_version(printHeader)

% print on screen the header as this is called at the beginning of a
% simulation
thisFile = which('get_neolco_version');
content = fileread(thisFile);

if nargin > 0
    if printHeader
        content = extractBefore(content,"function");
        disp(content);
    end
end

version = extractBetween(content,"Version: ","*");
version = strip(version," ");
version = version{1};

end
