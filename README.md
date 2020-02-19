# quTree

A tensor tree and general linear algebra package in C++.

## Getting Started

TODO: Basic example

## Installation
quTree can simply be installed using common package managers. In order to get quTree via apt, just write 
```
apt-get install qutree
``` 
or use a different package managers.

TODO: Instruction on compilation and/or installation through package managers

If you want to enable OpenMP for multi-threading on MacOS, do the following:
1) Install llvm and its version of openmp using homebrew:
    $ brew install llvm
2) In CMakeLists.txt, set(openmp ON)
3) When you run, $ export OMP_NUM_THREADS=\<desired no. threads\>

## Citation

If quTree is useful to your work, please cite the following paper:

TODO