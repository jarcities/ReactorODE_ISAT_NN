#!/bin/bash
CXX=clang++
$CXX -std=c++17 -O2 \
    -I"$CONDA_PREFIX/include" \
    -I"$CONDA_PREFIX/include/cantera" \
    -L"$CONDA_PREFIX/lib" \
    -Wl,-rpath,"$CONDA_PREFIX/lib" \
    -lcantera \
    -lsundials_cvodes \
    -lsundials_nvecserial \
    -lsundials_sunmatrixdense \
    -lsundials_sunlinsoldense \
    -lsundials_core \
    custom_sens.cpp -o custom_sens
./custom_sens