# Makefile

# List of the files (put the modules first) 
OBJ =  declarations.o main.o condition_Darcy.o read_inputs.o mesh.o timestep.o SMPP.o boundaries.o VTSWriter2D.o calc_parameters.o thinning_melt.o
OBJ2 = check_nan.o refreezing.o write_outputs.o solver_mass_heat.o conversion_grain.o hysteresis.o

# Name of the compiler
#F90 =  g95 
F90 =  gfortran 

# Flags
FLAGS = -g -O3 -Wall 

# Name of the executable
EXE = SMPP

# Edition of the links 
prog_principal:	$(OBJ) $(OBJ2)
	        $(F90) $(FLAGS) $(OBJ) $(OBJ2) -o $(EXE)


# Compilation of all the functions
declarations.o: declarations.f90
	$(F90) $(FLAGS) -c declarations.f90

main.o: main.f90
	$(F90) $(FLAGS) -c main.f90
	
SMPP.o: SMPP.f90
	$(F90) $(FLAGS) -c SMPP.f90

condition_Darcy.o: condition_Darcy.f90
	$(F90) $(FLAGS) -c condition_Darcy.f90


read_inputs.o: read_inputs.f90
	$(F90) $(FLAGS) -c read_inputs.f90

	
mesh.o: mesh.f90
	$(F90) $(FLAGS) -c mesh.f90

timestep.o: timestep.f90
	$(F90) $(FLAGS) -c timestep.f90

boundaries.o: boundaries.f90
	$(F90) $(FLAGS) -c boundaries.f90
	
	
calc_parameters.o: calc_parameters.f90
	$(F90) $(FLAGS) -c calc_parameters.f90
	
VTSWriter2D.o: VTSWriter2D.f90
	$(F90) $(FLAGS) -c VTSWriter2D.f90

thinning_melt.o: thinning_melt.f90
	$(F90) $(FLAGS) -c thinning_melt.f90	

check_nan.o: check_nan.f90
	$(F90) $(FLAGS) -c check_nan.f90	

refreezing.o: refreezing.f90
	$(F90) $(FLAGS) -c refreezing.f90	

write_outputs.o: write_outputs.f90
	$(F90) $(FLAGS) -c write_outputs.f90	

solver_mass_heat.o: solver_mass_heat.f90
	$(F90) $(FLAGS) -c solver_mass_heat.f90	

conversion_grain.o: conversion_grain.f90
	$(F90) $(FLAGS) -c conversion_grain.f90	

hysteresis.o: hysteresis.f90
	$(F90) $(FLAGS) -c hysteresis.f90


# Cleaning
clean:
	rm -f $(OBJ) $(EXE) $(OBJ2) *.mod
	cd ../outputs/water_output && rm -f *.vts *.pvd
	cd ../outputs/tempe_output && rm -f *.vts *.pvd
	cd ../outputs/densi_output && rm -f *.vts *.pvd
	cd ../outputs/dopti_output && rm -f *.vts *.pvd
	cd ../outputs/spher_output && rm -f *.vts *.pvd
	cd ../outputs/ && rm -f *.dat

