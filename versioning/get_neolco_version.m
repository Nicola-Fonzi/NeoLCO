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
