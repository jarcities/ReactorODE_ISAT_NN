#!/bin/bash

#fermi and oneapi stuff
# ml unload shared DefaultModules gcc/13.1.0 slurm/slurm/23.02.8 # FOR SDSU FERMI ONLY
# source ../../intel/oneapi/setvars.sh # FOR SDSU FERMI ONLY
# source $HOME/Code/intel/oneapi/setvars.sh # FOR LOCAL LINUX ONLY
 
#compile cpp with linking cantera
# icpx -O3 -L../cantera/build/lib -I../cantera/include -c reactor.cpp -lcantera
icpx -O3 -L$CONDA_PREFIX/lib -I$CONDA_PREFIX/include -c reactor.cpp -lcantera
#debugging
# icpx -g -fsanitize=address -O3 -L$CONDA_PREFIX/lib -I$CONDA_PREFIX/include -c reactor.cpp -lcantera


#compile f90 linking cantera and cpp 
# ifx -O3 -L../cantera/build/lib -I../cantera/include -o main.exe driver_samples.f90 reactor.o \
#	-lstdc++ -lcantera -Bstatic -L../ISAT/lib -I../ISAT/isatab_ser -lisat7_ser -L -mkl -Bdynamic -mkl -lstdc++
ifx -O3 -L$CONDA_PREFIX/lib -I$CONDA_PREFIX/include -o main.exe driver_samples.f90 reactor.o \
	-lstdc++ -lcantera -Bstatic -L../isat/ISAT/lib -I../isat/ISAT/isatab_ser -lisat7_ser -qmkl \
	-Bdynamic -qmkl -Wl,-rpath,$CONDA_PREFIX/lib

#run program
# ./main.exe
