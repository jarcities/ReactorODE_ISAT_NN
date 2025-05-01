#!/bin/bash
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
clang++ -std=c++23 -O3 -I"$CONDA_PREFIX/include" -L"$CONDA_PREFIX/lib" simple_state.cpp -o main -lcantera -Wl,-rpath,"$CONDA_PREFIX/lib"
