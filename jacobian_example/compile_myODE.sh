#!/bin/bash
gcc-15 -o myODE myODE.c -O3 -DNDEBUG -I/Users/jay/homebrew/opt/sundials/include -I/Users/jay/homebrew/opt/open-mpi/include -L/Users/jay/homebrew/opt/sundials/lib -L/Users/jay/homebrew/opt/open-mpi/lib -lsundials_cvodes -lsundials_nvecserial -lsundials_nvecmanyvector -lsundials_core -lmpi -lm -Wl,-rpath,/Users/jay/homebrew/opt/sundials/lib -Wl,-rpath,/Users/jay/homebrew/opt/open-mpi/lib
./myODE -sensi sim t
