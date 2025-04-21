
#!/bin/bash
rm *.o
rm *.mod
rm main.exe

# COMMENT OUT IF NECESSARY
# ml unload shared DefaultModules gcc/13.1.0 slurm/slurm/23.02.8 # FOR SDSU FERMI ONLY

# SOURCE VARIABLES FROM INTEL ONEAPI
# source ../../intel/oneapi/setvars.sh # FOR SDSU FERMI ONLY
source $HOME/intel/oneapi/setvars.sh # FOR LOCAL LINUX ONLY
 
# USED A CONDA ENV WITH CANTERA INSTALLED IN ENV AND OTHER SELF ADDED PATHS
#icpx -O3 -L../cantera/build/lib -I../cantera/include -c reactor.cpp -lcantera
icpx -O3 -L$CONDA_PREFIX/lib -I$CONDA_PREFIX/include -c reactor.cpp -lcantera
#ifx -O3 -L../cantera/build/lib -I../cantera/include -o main.exe driver_samples.f90 reactor.o \
#	-lstdc++ -lcantera -Bstatic -L../ISAT/lib -I../ISAT/isatab_ser -lisat7_ser -L -mkl -Bdynamic -mkl -lstdc++
ifx -O3 -L$CONDA_PREFIX/lib -I$CONDA_PREFIX/include -o main.exe driver_samples.f90 reactor.o \
	-lstdc++ -lcantera -Bstatic -L../isat/ISAT/lib -I../isat/ISAT/isatab_ser -lisat7_ser -qmkl \
	-Bdynamic -qmkl -Wl,-rpath,$CONDA_PREFIX/lib

# SOURCE VARIABLES AGIAN AND RUN EXE
source ../../intel/oneapi/setvars.sh
./main.exe
