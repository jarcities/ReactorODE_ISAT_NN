#﻿!/bin/bash

#export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
#export CANTERA_INCLUDE_DIR=/mnt/c/complete/cantera/include
#export CANTERA_INCLUDE_DIR=$CONDA_PREFIX/include
#export CANTERA_LIB_DIR=/mnt/c/complete/cantera/build/lib
#export CANTERA_LIB_DIR=$CONDA_PREFIX/lib
#export ISATAB_SERIAL_DIR=/mnt/c/complete/ISAT/isatab_ser
#export ISATAB_SERIAL_DIR=$HOME/Code/isat/ISAT/isatab_ser
#export ISAT_LIB_DIR=/mnt/c/complete/ISAT/lib
#export ISAT_LIB_DIR=$HOME/Code/isat/ISAT/lib

rm *.o
rm *.mod
rm PaSR.exe

#icpx -O3 -L$CANTERA_LIB_DIR -I$CANTERA_INCLUDE_DIR -c reactor.cpp -lcantera
#icpx -O3 -L$CONDA_PREFIX/lib -I$CONDA_PREFIX/include -c reactor.cpp -lcantera

#ifx -O3 -L$CANTERA_LIB_DIR -I$CANTERA_INCLUDE_DIR -o PaSR.exe PaSR.f90 reactor.o -lstdc++ -lcantera -Bstatic -L$ISAT_LIB_DIR -I$ISATAB_SERIAL_DIR -lisat7_ser -L -mkl -Bdynamic -mkl -lstdc++
#ifx -O3 -L$CONDA_PREFIX/lib -I$CONDA_PREFIX/include -o PaSR.exe PaSR.f90 reactor.o -lstdc++ -lcantera -Bstatic -L../isat/ISAT/lib -I../isat/ISAT/isatab_ser -lisat7_ser -qmkl -Bdynamic -qmkl -Wl,-rpath,$CONDA_PREFIX/lib

icpx -O3 -I$CONDA_PREFIX/include -c reactor.cpp
ifx -O3 \
    -L$CONDA_PREFIX/lib \
    -I$CONDA_PREFIX/include \
    -o PaSR.exe \
    PaSR.f90 reactor.o \
    -lstdc++ -lcantera \
    -Bstatic -L../isat/ISAT/lib -I../isat/ISAT/isatab_ser -lisat7_ser \
    -qmkl -Bdynamic -qmkl \
    -Wl,-rpath,$CONDA_PREFIX/lib
