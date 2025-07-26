#!/bin/bash
CXX=clang++
$CXX -std=c++17 -O2 \
    -I"$CONDA_PREFIX/include" -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" -lcantera \
    -I"$CONDA_PREFIX/include/eigen3" \
    custom_jac.cpp -o custom
./custom_jac