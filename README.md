# nextDFTB

Next-generation **DFTB** (Density Functional Tight-Binding) code, written in modern Fortran (2008).

## Features

- Modular source tree covering SCC, Hamiltonian, Coulomb, forces, dynamics, I/O, Slater–Koster tables, solvers and analysis.
- CMake-based build system with optional OpenMP and MPI parallelism.
- BLAS / LAPACK / HDF5 integration.
- CTest-driven test suite.

## Requirements

- CMake ≥ 3.20
- A Fortran compiler supporting Fortran 2008 (e.g. `gfortran`, `ifort`, `ifx`)
- BLAS and LAPACK
- HDF5 with Fortran bindings (`libhdf5-fortran-dev` on Debian/Ubuntu)
- *(optional)* MPI implementation (OpenMPI, MPICH, …)
- *(optional)* OpenMP-capable compiler

On Debian/Ubuntu:

```bash
sudo apt install cmake gfortran libblas-dev liblapack-dev libhdf5-dev libhdf5-fortran-102
```

## Build

```bash
# Configure
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# Compile
cmake --build build -j
```

### Configuration options

| Option                    | Default | Description                    |
| ------------------------- | ------- | ------------------------------ |
| `CMAKE_BUILD_TYPE`        | Release | `Debug`, `Release`, `RelWithDebInfo` |
| `NEXTDFTB_ENABLE_MPI`     | OFF     | Enable MPI parallelism         |
| `NEXTDFTB_ENABLE_OPENMP`  | OFF     | Enable OpenMP parallelism      |
| `NEXTDFTB_BUILD_TESTS`    | ON      | Build the CTest test suite     |

Example with OpenMP and debug symbols:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DNEXTDFTB_ENABLE_OPENMP=ON
cmake --build build -j
```

## Run

A minimal sanity-check binary is provided:

```bash
./build/src/hello
```

## Tests

```bash
ctest --test-dir build --output-on-failure
```


## License

Released under the [MIT License](LICENSE).
