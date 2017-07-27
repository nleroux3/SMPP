# SMPP

SMPP (Snow Model with Preferential flow Paths) is a 2D snowmelt model that estimates heat and mass fluxes through a snowpack. 
It includes the formation of preferential flow paths.
Further description of the model is available in Leroux and Pomeroy (2017)


Reference:
Leroux, N. R., and Pomeroy, J. W.: Modelling capillary hysteresis effects on preferential flow through melting and cold layered snowpacks. Advances in Water Resources, 107, 250-264, 2017

## Getting Started

SMPP is coded in Fortran 90. A linux executable `SMPP` is produced by running the script 'compil.sh'. 
Running the script 'compil.sh' will also create the output folders if they are not present.

### Prerequisites

The [gfortran] (https://gcc.gnu.org/wiki/GFortran) compiler is used to compile the model.
The Makefile in the folder src can be edited to use other compilers. 

## Running the tests

A test can be run by running the executable "SMPP". 
The input file "inputs.txt" contains the lab snow data from Waldner et al. (2014).
This file can be modified to model rain on snow or snowmelt of any structured snowpack.

## Model outputs

All model ouputs are written in the folder "outputs".
The folders "density", "grain_d", tempera" and "theta_w" contain 2D model outputs (density, optical grain diameter, temperature and water content, respectively) that can be visualized with Paraview.
Other outputs file are written in this folder, such as "Tss.dat" containing snow surface temperature results and "outflow.dat" containing model outflow results.

## Cleaning model ouputs

By running the script 'clean.sh', the files in the folder 'outputs' will be removed.

## Author

* **Nicolas Leroux** - [nleroux3](https://github.com/nleroux3)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


