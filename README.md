# NeoLCO

## How to use the code

NeoLCO is a NeoCASS add-on for structurally nonlinear simulations. In order to use the code, you must have access to Matlab and Simulink, as the entire library is designed
to work in that environment. You also need to download NeoCASS from the official repository.
Only NeoCASS version 3.0.0 or higher is supported.

## Scope of the code

The framework is intended to be used as a general-purpose tool to study LCO and other nonlinear phenomena when concentrated nonlinearities are preset in the system.
The idea is to start from a linear structure, in Nastran/NeoCASS format, where some selected locations, either grid points or scalar points, are left completely free from
stiffness, support or definition of an aerodynamic surface. These points are then identified as nonlinearity locations and will later be used to introduce
concentrated nonlinearities.

At that point, different simulations can be performed semi-automatically to obtain the different nonlinear behaviours.
