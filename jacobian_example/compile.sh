# #!/bin/bash
# CXX=clang++
# $CXX -std=c++17 -O2 \
#     -I"$CONDA_PREFIX/include" -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" -lcantera \
#     -I"$CONDA_PREFIX/include/eigen3" \
#     custom__.cpp -o custom
# ./custom














#!/usr/bin/env bash
set -e

CXX=clang++

# Conda‐provided include & lib dirs
INC_DIR="$CONDA_PREFIX/include"
EIGEN_DIR="$INC_DIR/eigen3"
LIB_DIR="$CONDA_PREFIX/lib"

$CXX -std=c++17 -O2 \
    -I"$INC_DIR" \
    -I"$EIGEN_DIR" \
    custom__.cpp -o custom \
    -L"$LIB_DIR" \
    -Wl,-rpath,"$LIB_DIR" \
    -lcantera \
    -lsundials_cvode \
    -lsundials_cvodes \
    -lsundials_nvecserial \
    -lsundials_sunmatrixdense \
    -lsundials_sunlinsoldense


