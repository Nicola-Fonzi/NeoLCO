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
function version = get_neolco_version(printHeader)

% print on screen the header as this is called at the beginning of a
% simulation
thisFile = which('get_neocass_version');
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
