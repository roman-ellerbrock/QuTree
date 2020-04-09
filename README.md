# QuTree

A tensor tree linear algebra package in C++ designed for quantum dynamics and machine learning applications.

## Installation

QuTree can simply be installed using several package managers.

On OS X, QuTree is available via a HomeBrew tap: 
```
brew tap sseritan/qu-tree
brew install qu-tree
```

TODO: Apt and yum

Tested operating systems:
* OS X 10.14

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

Example CMakeLists.txt:
```
cmake_minimum_required(VERSION 3.0)
project(example CXX)

find_package(QuTree REQUIRED)

add_executable(test test.cpp)
target_link_libraries(test QuTree::QuTree)
```

Example `test.cpp`:
```
#include <Core/Tensor.h>
#include <Core/Matrix.h>

int main()
{
TensorShape tdim({2, 3, 4});
Tensorcd A(tdim);
for (size_t i = 0; i < A.shape().totalDimension(); i++) {
    A(i) = i;
}
Matrixcd w = A.DotProduct(A);
w.print();
}
```

For Makefiles, set your include and library paths to `$PREFIX/include` and `$PREFIX/lib`, respectively.

For detailed examples on how to use the library, please see the `examples` folder.
For more extensive documentation, please see our [ReadTheDocs page](https://qutree.readthedocs.io/en/latest/).


## Citation

If QuTree is useful to your work, please cite the following paper:

TODO