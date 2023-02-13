[![Ubuntu Clang Deployment](https://github.com/roman-ellerbrock/QuTree/actions/workflows/ubuntu-clang.yml/badge.svg)](https://github.com/roman-ellerbrock/QuTree/actions/workflows/ubuntu-clang.yml)

[![Ubuntu GCC Deployment](https://github.com/roman-ellerbrock/QuTree/actions/workflows/ubuntu-gcc.yml/badge.svg)](https://github.com/roman-ellerbrock/QuTree/actions/workflows/ubuntu-gcc.yml)

# QuTree

A tensor tree linear algebra package in C++ designed for quantum dynamics and machine learning applications.

## Getting Started

Installation is easy using HomeBrew (on OS X) or LinuxBrew (on Linux):
```
brew tap sseritan/qu-tree
brew install qu-tree
```

Developers can get access to the cutting-edge package by using the `--HEAD` flag to install directly from the HEAD of master.
Subsequent builds can be upgraded by using the `--fetch-HEAD` flag.
```
brew install --HEAD qu-tree
```
Subsequent builds can be upgraded by using the `--fetch-HEAD` flag.
```
brew upgrade --fetch-HEAD qu-tree
```

For instructions on how to install from source, check out our [documentation](https://qutree.readthedocs.io/en/latest/).

## Using QuTree

QuTree can be used as a library or via it's standalone applications like the multiconfigurational timedependent 
Hartree (mctdh) application or the quantum virtual machine (qvm).

Applications like mctdh & qvm can be found in ${project}/contrib. Examplary inputs are located in 
${project}/examples/${application_name}. Applications are run via
```
./mctdh {input.yaml} > {output.txt}
```

After installation, QuTree can be easily used in downstream CMake projects.
If installed to non-standard locations, make sure to set `QuTree_DIR` to the location of `QuTreeConfig.cmake`.

Example CMakeLists.txt:
```
cmake_minimum_required(VERSION 3.0)
project(example CXX)

find_package(QuTree REQUIRED)

add_executable(qutree_app app.cpp)
target_link_libraries(qutree_app QuTree::QuTree)
```

Example `app.cpp`:
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
Matrixcd w = A.dotProduct(A);
w.print();
}
```

For detailed examples on how to use the library, please see the `examples` folder
 or check out our [documentation](https://qutree.readthedocs.io/en/latest/).

## Citation

If QuTree is useful to your work, please cite the following paper:

R. Ellerbrock, K. G. Johnson, S. Seritan, H. Hoppe, H. J. Zhang, T. Lenzen, T. Weike, U. Manthe, T. J. Martínez,
"QuTree - a Tree Tensor Network package", 2023 (in preparation)
