============
Installation
============

Requirements
============

* CMake >= 3.0
* a C++17 compiler
* yaml-cpp (only for mctdh & qvm)
* Boost (only for qvm)

You can install yaml-cpp from source via::

    git clone https://github.com/jbeder/yaml-cpp.git
    cd yaml-cpp
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=/home/[username]/usr/

Compiling from Source
=====================

Compilation follows the standard CMake procedure::

    mkdir build
    cd build
    cmake ../
    make
    make install

For testing, run :code:`make TestQuTree` instead.
Applications are compiled as targets on request. They will NOT be compiled on a standard make call since they
have additional requirements. Run :code:`make mctdh` or :code:`make qvm` to create the targets but note that they require
yaml-cpp (for mctdh & qvm) and Boost (for qvm).

If you work on a system with specialized installation pathways, you can tell cmake to choose the correct
implementations::

    cmake .. -DQuTree_DIR=[path]
    make
    make mctdh
    make qvm

Switch out the :code:`[path]` with the path to your preferred cmake installation directory.
If cmake cannot find required packages, you can select them by adding :code:`-Dyaml-cpp_DIR=[path]/yaml-cpp/`
for yaml-cpp, :code:`-DBOOST_ROOT=[path]/usr/` for boost. If blas/lapack includes files are not found, please check
the system's environment variables.
If the problems occur at the linking step, then the paths set to libraries might not be set correctly.
The library paths for Lapack & Blas can be set by adding :code:`-DLAPACKE_PATH=[path]`, :code:`-DBLAS_LIBRARIES=[path]`
to the cmake call. The paths should be set to directories
that contain the libraries of Blas & Lapack. Make sure that the libraries were compiled using the same compiler, otherwise there
might be missing functions due to differing name decorations.

HomeBrew/LinuxBrew
==================

QuTree is available via HomeBrew on OSX or LinuxBrew on Linux distributions.

Install HomeBrew::

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
and follow any post-installation steps.

Then install the QuTree package is distributed via a HomeBrew tap::

    brew tap sseritan/qu-tree
    brew install qu-tree

On Debian/Ubuntu, you may need the following packages for LinuxBrew::

    apt install build-essential curl file git locales

On RHEL/CentOS, you may need the following packages for LinuxBrew::

    yum groupinstall 'Development Tools'
    yum install curl file git which perl-core

Tested operating systems:

* OS X 10.14
* Ubuntu 18.04 LTS
* CentOS 7

OpenMP Support
--------------
If you want to enable OpenMP for multi-threading on OS X, do the following:
1) Install :code:`llvm` and its version of :code:`openmp` using HomeBrew::

    brew install llvm

2) In CMakeLists.txt, :code:`set(openmp ON)`
3) When you run, :code:`export OMP_NUM_THREADS=<desired no. threads>`

