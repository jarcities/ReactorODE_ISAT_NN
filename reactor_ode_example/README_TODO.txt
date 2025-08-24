.......................
HOW TO RUN BASH SCRIPT:
^^^^^^^^^^^^^^^^^^^^^^^
1/ Installed miniconda3 and created a virtual env
2/ Installed intel oneapi base and hpc toolkit in home dir (I PUT IT IN $HOME/CODE)
3/ Installed cantera (c++/fortran version) using conda in that env
4/ Set env variables using setvars.sh in intel oneapi dir
5/ Changed paths to link ISAT which shares same parent dir as nn_test
6/ Changed paths to link cantera to conda version instead of ISAT version
7/ Comments are made in build.sh 

.........................
LOCAL LINUX INSTRUCTIONS:
^^^^^^^^^^^^^^^^^^^^^^^^^
1/ `conda activate cantera` #dedicated conda env for cantera c++
2/ `./build.sh` #bash script runs setvars.sh


......
TODOS:
^^^^^^
1/ Figure out canteras jacobian.
2/ Figure out how to implement back propagation.
3/ fix cantera conda env
4/ maybe do cantera + intel conda env
5/ finish transferring fnn to codejenn

...........
REFERENCES:
^^^^^^^^^^^
#cantera distro
1/ https://cantera.org/stable/reference/

#how cantera integrates
2/ https://cantera.org/3.1/develop/reactor-integration.html

#sundials
3/ https://sundials.readthedocs.io/en/latest/cvodes/Usage/FSA.html#user-supplied-routines-for-forward-sensitivity-analysis
4/ https://sundials.readthedocs.io/en/latest/cvodes/Mathematics_link.html#forward-sensitivity-analysis

#custom reactor example
5/ https://cantera.org/3.1/examples/cxx/custom.html#sphx-glr-examples-cxx-custom-cpp

#sensitivity function class
6/ https://cantera.org/3.1/cxx/d0/d1c/ReactorNet_8cpp_source.html