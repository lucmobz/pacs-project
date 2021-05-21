# Wave Equation Library

## Information
This is a C++ library to solve the scalar linear 3D wave equation problem (for example an acoustic problem) in a domain with suitable Dirichlet or Neumann boundary conditions and Cauchy initial conditions. The numerical method used is the DG method (discontinuous Galerkin) to discretize the problem in space and the finite difference method to discretize in time (there is a built-in solver that does the second order explicit leap-frog scheme). The variant of the method is the modal DG method on polyhedral grids (the basis functions are tensor product Legendre polynomials restriced to the polyhedrons).
The library uses an *expression-templates* system in order to give the mathematical expression of the variational formulation as written on paper as input.

## Folders
The folder structure is explained below:

* `./build` will contain all the built .o and .d files from the compilation
* `./core` contains all the files for the main algorithm and expression-template system
* `./doc` will contain the documentation after doing `make doc`
* `./external` contains the **GetPot** header
* `./functional` contains the  files for the quadrature rules and basis functions implementation
* `./mesh` contains the files for the mesh implementations
* `./tests` contains the files for the executables (`./tests/data` contains the meshfiles)
* `./utilities` contains some helper functions

## Requirements
The library was developed using some **C++17** and **C++20** features:

* The compiler used was `g++-10 (Ubuntu 10.1.0-2ubuntu1~18.04) 10.1.0` with the `-std=c++20` flag
* It uses these **C++20** headers: `<span>` and `<type_traits>` and some of the *constexpr* standard library functions
* It uses the **Eigen** library, version **3.3.7**
* It uses the internal **openmp** support for the **Eigen** library (by using the `-fopenmp` flag)
* It uses the **CGAL** library, version **5.0.1**
* It uses the header `<boost/math/quadrature/gauss.hpp>` from the **boost** library, version **1.72.0**
* Some exectuables use the header **GetPot** (included in `./external`) for command line and file options parsing
* The documentation is built with **Doxygen**

Everything else was as in the toolchain used in the course.

## Compilation
There is a `Makefile` in the main folder. Typing `make help` inside the main folder prints the make commands to compile (or just type `make all` to build all executables, after that run a quick test with `make run` or go inside `./tests` and run them as described below). 

**Inside the `Makefile` some variables must be changed:**

* **CXX** must be set to a compiler that supports the **C++20** features mentioned (only `g++-10` resulting from `sudo apt install g++-10` was tested on the ubuntu VM)
* **THIRD_PARTY_INC** must contain the locations of the **Eigen**, **CGAL** and **boost** libraries (for example `THIRD_PARTY_INC = -I$mkEigenInc -I$mkCgalInc -I$mkBoostInc`)
* (*optional*) setting **OPTCPPFLAGS = -DVERBOSE** prints some additional timings 
* `make doc` builds the documentation in *html* format inside the `./doc` folder
* `make distclean` removes all compiled and additional files
* `make shared` or `make static` will create in the current directory a folder `./wave/lib` with the shared and static version of the library (soname is `libwave.so`), the headers will be copied in `./wave/include`, the global include file is `wave.hpp` (building the libraries is independent from the tests)

## Run
To run the executables enter the `./tests` folder and type:

* `./hconv --file hconv.pot` to run the h-convergence test (3 to 5 minutes)
* `./pconv --file pconv.pot` to run the p-convergence test (2 to 3 minutes)
* `./test ./data/voro1000.voro` to run a sample test on a specific meshfile (other meshes are inside the `./tests/data` folder)