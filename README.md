## Quantum Magnets Large Scale Simulations

This code provides access to the ground state, and to excites states as well, of spin-1/2 systems.
The most general Hamiltonian that can be used is the Heisenberg one with nearest and next-nearest 
neighbors interactions. Anisotropy and disorder along the magnetization axis can also be simulated.

The spins may live in a one-dimensional chain lattice or in a two-dimensional square lattice. However, 
extending the available lattices is possible.

This code has as final purpose the study of frustrated magnetic systems, and for that it includes 
the computation of several parameters: order parameter (magnetization), correlation functions and 
dynamical structure factor.

### Compiling
For this code you need:
 1. [Ninja](https://ninja-build.org/)
 1. [Meson](https://mesonbuild.com/)
 2. [PETSc](https://www.mcs.anl.gov/petsc/) (in its complex versions with 64-bit integers)
 3. [SLEPc](http://slepc.upv.es/)
 4. [Boost](https://www.boost.org/)
 5. an MPI implementation

##### Compile not using Meson

```
$ mpicxx main.cpp ./src/*.cpp ./src/spinOperators/*.cpp ./src/tools/*.cpp `pkg-config --libs --cflags SLEPc` -I./include -I./include/spinOperators -I./include/tools -I/path-to-boost-lib/ -I/path-to-mpi-impl -o main.x
```
If you are not using `pkg-config` then you can link manually with PETSc and SLEPc libraries by adding the include directories and the library directories:
```
-I/path-to-slepc-build/include -I/path-to-petsc-build/include -L/path-to-slepc-build/lib -L/path-to-petsc-build/lib -lslepc -lpetsc
```

##### Compile using Meson
```
$ CXX=<mpi c++ compiler> CXXFLAGS='-I/path-to-boost-root-dir/ -I/path-to-the-mpi.h-file' meson dir/where/to/build
$ cd dir/where/to/build
$ ninja
```

### Running
To run a simple example, with the following setting:
 1. no disorder
 2. with nearest neighbours interactions only
 3. no anisotropy in the z-direction
 2. computing 10 excited states for the dynamical structure factor
 3. for `nspins` spins
 4. with `N` processes
 
```
$ cd dir/where/to/build
$ mpiexec -np N ./main.x nspins -nn j1 1.0 -d1 1.0 -eps_type krylovschur -eps_tol 1e-9
```
