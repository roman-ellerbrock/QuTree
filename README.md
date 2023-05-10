# QuTree

A tensor network package for solving ordinary and partial differential equations.

## Requirements

The package depepends on

* C++-17 compiler (clang)
* cmake
* torch (https://pytorch.org/cppdocs/)
* json-nlohman (https://github.com/nlohmann/json)
* yaml-cpp (https://github.com/jbeder/yaml-cpp)

The code is currently tested on MacOS with clang.

## Getting Started

After installing the required dependencies, the package is installed via

```
mkdir -p build
cd build
cmake ..
make
```

and the executable is found in the build directory.

## Running the code

You can find examplary input files in ./examples.
Execute the examples via

```
./qutree binary.yaml
```

Output is written to std::cout.

## Contributing

We appreciate that you help us writing this code!
If you want to contribute to the project, here are a few steps to follow:

1. Familiarize Yourself with the Code:
   * Read the project documentation, including the README file, to understand the project's purpose, scope, and guidelines.
   * Explore the existing codebase to gain an understanding of its structure and conventions
2. Identify and Discuss the Contribution:
   * Determine the area of the project you want to contribute to, such as fixing a bug, implementing a feature, or improving documentation.
   * Check the project's issue tracker or communication channels (e.g., GitHub issues, mailing list, or chat platform) to see if there are any related discussions or tasks.
   * If necessary, start a discussion with the project maintainers or other contributors to ensure your contribution aligns with the project's goals and guidelines.
3. Create a Branch and Implement the Changes:
   * Create a new branch in git called "feature-*something*" or "bug-fix-*something*", etc. Follow the conventions of the package.
   * Implement your changes following the project's coding style and conventions.
   * Write tests to cover the new code or modifications you've made, ensuring the changes maintain or improve the project's overall code quality.
4. Submit a Pull Request (PR), and merge into main branch:
   * Rebase the project onto the most recent version
   * Submit a pull request to propose your changes for review.
   * Provide a clear and concise description of your changes in the PR, explaining the problem you're addressing and the solution you're proposing.
   * Participate in the review process by addressing comments, making necessary changes, and providing additional context or documentation as requested.
   * Celebrate your contribution! Your changes are now part of the project, and you can move on to the next task or contribute further.

## Goals

New contributions should be checked against our general requirements:

* Is the code well documented?
* Does it follow our style-guidelines?
* Does the new code come with appropriate tests and does it pass all tests?
* Can the code be extendable to new hardwares in the future?
* Are hardware specific pieces separated into modules that can be turned on and off via compile options?
* Can the code easily be shipped to different architechures (HPC clusters, Cloud, workstations, etc.)?
* Does the code achieve high performance and can it be straightforwardly accelerated?

If not mentioned otherwise, we strive for following the standard C++ guidelines:

* see [Bjarne Stroustrup's guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)
* we use the [Google style format for C++](https://google.github.io/styleguide/cppguide.html) (use clang-format or similar software)

## Contributors

Reach out to romanellerbrock@gmail.com for questions.

* Roman Ellerbrock
* K. Grace Johnson
* J. Harry Zhang,
* Hannes Hoppe
* Tim Lenzen
* Thomas Weike
* Uwe Manthe
* Todd J. Martínez

