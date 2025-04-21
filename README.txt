
HOW TO RUN BASH SCRIPT:

1/ Installed miniconda3 and created a virtual env
2/ Installed intel oneapi base and hpc toolkit in home dir (I PUT IT IN $HOME/CODE)
3/ Installed cantera (c++/fortran version) using conda in that env
4/ Set env variables using setvars.sh in intel oneapi dir
5/ Changed paths to link ISAT which shares same parent dir as nn_test
6/ Changed paths to link cantera to conda version instead of ISAT version
7/ Comments are made in build.sh 


LOCAL LINUX INSTRUCTIONS:
1/ `conda activate cantera` #dedicated conda env for cantera c++
2/ `./build.sh` #bash script runs setvars.sh