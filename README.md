# HeatTest

This repository is for the heat test required for the computational physics internship. 

## Model
For simplicity I chose to model this system as heat transfer throughout a rod of 
lenght L with lateral insulation where the cross-sectional area of the rod is 
thicker at the storage tank. Thus the density, thermal conductivities, and 
specific heats vary as a function of x. For continuity I prescribe that the heat
at the endpoints be equal, i.e. U(xl) = U(xr). Additionally, the pump direction
and the solar panel give a heat flux rightward through the rod so
d/dx U(xl) = f(t) where f(t) > 0 during the day, and negative at night. It is
my hope that this model captures the behavior of the system.

## Simulation
For this system, I will use backward Euler to ensure consistency and stability 
considerations of the solution. 
## Dependencies
This project depends on an installation of `Python 3.7`. The standard `numpy` 
library is assumed installed.

## Running