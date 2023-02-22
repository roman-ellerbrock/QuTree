[![CI](https://github.com/roman-ellerbrock/QuTree/actions/workflows/test.yml/badge.svg)](https://github.com/roman-ellerbrock/QuTree/actions/workflows/test.yml)
[![Docker build](https://github.com/roman-ellerbrock/QuTree/actions/workflows/docker.yml/badge.svg)](https://github.com/roman-ellerbrock/QuTree/actions/workflows/docker.yml)


# QuTree

A tensor tree linear algebra package in C++ designed for quantum dynamics and machine learning applications.

## Getting Started

There are different ways to install QuTree library and corresponding executables.
The easiest way to try out the executables would be to use a Docker container.
You can always install QuTree from source.


### Docker container of the executables

To obtain a Docker container with precompiled executables, type:

```bash
docker pull mtzgroup/qutree:latest
```

To run the executables in Docker:

```bash
# Example 1: running an mctdh example
docker run \
	-it \
	--rm \
	-v $(pwd)/examples/mctdh:/app \
	mtzgroup/qutree:latest \
	mctdh portfolioOptimization.nasdaq.25.yaml

# Example 2: mount your own input dir, then run mctdh interactively
mkdir my-inputs
docker run \
	-it \
	--rm \
	-v $(pwd)/my-inputs:/app \
	mtzgroup/qutree:latest
# You can run `mctdh {input.yaml} > {output.txt}` inside the container.
```

### Homebrew

To install QuTree using HomeBrew (on OS X) or LinuxBrew (on Linux), run:
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

After installation, the QuTree library can be easily used in downstream CMake projects.
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

R. Ellerbrock, K. G. Johnson, S. Seritan, H. Hoppe, J. H. Zhang, T. Lenzen, T. Weike, U. Manthe, T. J. Martínez,
"QuTree - a Tree Tensor Network package", 2023 (in preparation)
