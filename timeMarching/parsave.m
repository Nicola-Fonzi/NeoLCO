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
function parsave(fname, modes, quasiSteadyAerodynamicForces, unsteadyAerodynamicForces, torque, rotationInTime, dynamicPressure, tout)
  save(fname, "modes", "quasiSteadyAerodynamicForces", "unsteadyAerodynamicForces","torque","rotationInTime","dynamicPressure","tout")
end