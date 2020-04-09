============
Installation
============

Prerequisites include:

* CMake >= 3.0
* Eigen3
* UnitTest++ (only for testing)

If installing via HomeBrew, these are automatically handled for you.
Otherwise, QuTree will attempt to use `git` submodules to build them locally.

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

Compiling from Source
=====================

If Eigen or UnitTest++ are not detected on your system, QuTree will use `git` submodules to build them locally.

Compilation follows standard CMake procedure::

    mkdir build
    cd build
    cmake ../
    make
    make install

For testing, run :code:`make TestQuTree` instead.

OpenMP Support
--------------
If you want to enable OpenMP for multi-threading on OS X, do the following:
1) Install :code:`llvm` and its version of :code:`openmp` using HomeBrew::

    brew install llvm

2) In CMakeLists.txt, :code:`set(openmp ON)`
3) When you run, :code:`export OMP_NUM_THREADS=<desired no. threads>`
