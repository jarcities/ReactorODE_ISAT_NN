#!/bin/bash
CXX=clang++
$CXX -std=c++17 -I"$CONDA_PREFIX/include" -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" custom.cpp -lcantera -o custom