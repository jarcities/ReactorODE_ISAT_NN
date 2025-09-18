#!/bin/bash
rm *.o
rm *.mod
rm *.exe

icpx -O3 -I$CONDA_PREFIX/include -c reactor.cpp 

ifx -O3 \
    -L$CONDA_PREFIX/lib \
    -I$CONDA_PREFIX/include \
    -o PaSR.exe \
    PaSR.f90 reactor.o \
    -lstdc++ -lcantera \
    -Bstatic -L$HOME/code/isat/ISAT/lib -I$HOME/code/isat/ISAT/isatab_ser -lisat7_ser \
    -qmkl -Bdynamic -qmkl \
    -Wl,-rpath,$CONDA_PREFIX/lib