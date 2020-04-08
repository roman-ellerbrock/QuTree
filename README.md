# QuTree

A tensor tree linear algebra package in C++ designed for quantum dynamics and machine learning applications.

## Installation

QuTree can simply be installed using HomeBrew on OS X, `apt` on Debian/Ubuntu, and `yum` on RHEL/CentOS. 

```
brew install qutree
apt install qutree
yum install qutree
```

Tested operating systems:
* OS X 10.14
* Ubuntu 16.04 LTS
* CentOS 7

## Compiling From Source

Prerequisites include:
* CMake >= 3.0
* Eigen3
* UnitTest++ (only for testing)

If Eigen or UnitTest++ are not detected on your system, QuTree will use `git` submodules to build them locally.

Compilation follows standard CMake procedure:
```
mkdir build
cd build
cmake ../
make
make install
```
For testing, run `make TestQuTree` instead.

### OpenMP Support
If you want to enable OpenMP for multi-threading on OS X, do the following:
1) Install `llvm` and its version of `openmp` using HomeBrew:
```
brew install llvm
```
2) In CMakeLists.txt, `set(openmp ON)`
3) When you run, `export OMP_NUM_THREADS=<desired no. threads>`

## Getting Started

After installation, QuTree can be easily used in downstream CMake projects.

In your CMake project, 
```
find_package(QuTree REQUIRED)

target_link_libraries(<your executable> QuTree::QuTree)
```
For Makefiles, set your include and library paths to `$PREFIX/include` and `$PREFIX/lib`, respectively.

For detailed examples on how to use the library,
please see the `examples` folder or further documentation here.

TODO: Link to RTD

## Citation

If QuTree is useful to your work, please cite the following paper:

TODO