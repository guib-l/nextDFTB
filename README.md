# nextDFTB

Next-generation **DFTB** (Density Functional Tight-Binding) code, built as a
hybrid **Python / C++20 / Fortran 2018** stack with OpenMP and a deferred
MPI layer.

## Architecture at a glance

```
          ┌─────────────────────┐
          │  Python (pybind11)  │   Ergonomic API. NumPy ↔ F-order auto-copy.
          └──────────┬──────────┘
                     │
          ┌──────────▼──────────┐
          │   C++20 Orchestrator │  Aligned allocation, simulation loop,
          │       (core_cpp)     │  logging sink, error state, exceptions.
          └──────────┬──────────┘
                     │   ← C ABI (abi/)  iso_c_binding types only
          ┌──────────▼──────────┐
          │   Fortran 2018      │   Production compute engine.
          │      (engine)       │   OpenMP-parallel kernels.
          └─────────────────────┘

          ┌─────────────────────┐
          │   CLI (cli/)        │   Standalone nextdftb executable.
          └─────────────────────┘   Links core_cpp, runs without Python.

          ┌─────────────────────┐
          │   HPC / MPI         │   DESIGN ONLY — implementation deferred.
          │   (hpc_mpi/)        │   Cartesian topology, halo exchange,
          └─────────────────────┘   collective reductions.
```

### Source tree

| Directory           | Role                                                        |
|---------------------|-------------------------------------------------------------|
| `abi/`              | Public C ABI headers (`nextdftb_*`)                          |
| `engine/`           | Fortran 2018 compute engine with BIND(C) facade              |
| `core_cpp/`         | C++20 orchestrator, aligned buffers, logger, error state     |
| `bindings/`         | pybind11 extension module `nextdftb._core`                   |
| `cli/`              | Standalone CLI executable `nextdftb`                          |
| `io/`               | Placeholder I/O abstraction (text writer; HDF5 deferred)     |
| `hpc_mpi/`          | Design-only MPI interfaces (stubs throw when linked)         |
| `config/`, `logging/`, `docs/` | Reserved for future artefacts                     |
| `tests/`            | CTest-driven pytest runners                                  |
| `examples/`         | End-user examples (Python)                                   |
| `cmake/`            | Shared CMake modules                                         |

## Requirements

- CMake ≥ 3.20
- A C11 compiler and a C++20 compiler (GCC ≥ 12, Clang ≥ 15, or ICX)
- A Fortran 2018 compiler (`gfortran` ≥ 11, `ifx`, or NVHPC)
- LAPACK (`liblapack-dev` on Debian/Ubuntu)
- Python 3.9+ with `numpy` and `pybind11` (only needed when
  `NEXTDFTB_ENABLE_PYTHON=ON`, the default — see `requirements.txt`)
- *(optional)* An OpenMP-capable compiler
- *(optional)* An MPI implementation — **not required**, the MPI layer is
  design-only in this release

## Install

Three paths are supported depending on what you need.

### 1. `pip install .` (Python module only)

Easiest if you just want `import nextdftb`. A `pyproject.toml` drives CMake
via [scikit-build-core](https://scikit-build-core.readthedocs.io); it
disables the CLI and the test suite and installs `nextdftb` into your
active Python environment.

```bash
sudo apt install cmake gfortran g++ python3-dev liblapack-dev
python3 -m venv .env && source .env/bin/activate
pip install -e .
python -c "import nextdftb; nextdftb.init(); print(nextdftb.test()); nextdftb.finalize()"
```

### 2. Full build (CLI + Python + tests)

```bash
sudo apt install cmake gfortran g++ python3-dev liblapack-dev
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt        # numpy + pybind11 (build deps)

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DPython3_EXECUTABLE="$(which python3)"
cmake --build build -j
cmake --install build                  # CLI → $prefix/bin, module → site-packages
```

### 3. CLI only (no Python)

```bash
sudo apt install cmake gfortran g++ liblapack-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNEXTDFTB_ENABLE_PYTHON=OFF
cmake --build build -j
sudo cmake --install build
```

## Build (without installing)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Options

| Option                      | Default | Description                                  |
|-----------------------------|---------|----------------------------------------------|
| `CMAKE_BUILD_TYPE`          | Release | `Debug`, `Release`, `RelWithDebInfo`         |
| `NEXTDFTB_ENABLE_OPENMP`    | ON      | OpenMP in the Fortran kernels                |
| `NEXTDFTB_ENABLE_MPI`       | OFF     | Build the deferred MPI layer as stubs        |
| `NEXTDFTB_ENABLE_PYTHON`    | ON      | Build the pybind11 module                    |
| `NEXTDFTB_BUILD_CLI`        | ON      | Build the `nextdftb` CLI                     |
| `NEXTDFTB_BUILD_TESTS`      | ON      | Build the CTest suite                        |
| `NEXTDFTB_BUILD_FORTRAN_MAIN` | OFF   | Build the native Fortran driver (dev tool)   |

## Run

### CLI

Runs an end-to-end infrastructure smoke test through every layer
(C++ orchestrator → C ABI → Fortran + OpenMP) and prints a deterministic
reference value (`500500`, the closed-form of `1+2+…+1000`).

```bash
./build/cli/nextdftb --log-stdout
```

Available flags: `--threads N`, `--log PATH`, `--log-level INFO|DEBUG|WARN|ERROR`,
`--log-stdout`, `--output PATH`, `--version`, `--help`.

### Python

```bash
PYTHONPATH=build/python python -c "
import nextdftb
nextdftb.init()
print('test =', nextdftb.test())
print('version:', nextdftb.version())
nextdftb.finalize()
"
```

The public Python surface is currently `init`, `finalize`, `is_initialized`,
`test`, `version`, and the `log` submodule. Numerical kernels (saxpy, L2
norm, …) are not yet exposed — only the infrastructure smoke test is wired
end-to-end at this stage.

## Tests

```bash
ctest --test-dir build --output-on-failure
```

## License

Released under the [MIT License](LICENSE).
