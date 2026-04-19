cmake_minimum_required(VERSION 3.20)

project(nextdftb
    VERSION     0.1.0
    DESCRIPTION "Next-generation DFTB — hybrid Python / C++20 / Fortran 2018 with OpenMP (MPI layer deferred)"
    LANGUAGES   C CXX Fortran
)

# ---------------------------------------------------------------------
# Language standards (Fortran standard is enforced via compiler flags —
# CMAKE_Fortran_STANDARD is not used because CMake 3.20 doesn't know
# Fortran 2018).
# ---------------------------------------------------------------------
set(CMAKE_C_STANDARD            11 )
set(CMAKE_C_STANDARD_REQUIRED   ON )
set(CMAKE_C_EXTENSIONS          OFF)
set(CMAKE_CXX_STANDARD          20 )
set(CMAKE_CXX_STANDARD_REQUIRED ON )
set(CMAKE_CXX_EXTENSIONS        OFF)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# ---------------------------------------------------------------------
# Build type default
# ---------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo)
endif()

# ---------------------------------------------------------------------
# User-facing options
# ---------------------------------------------------------------------
option(NEXTDFTB_ENABLE_OPENMP "Enable OpenMP in the Fortran kernels"              ON )
option(NEXTDFTB_ENABLE_MPI    "Enable the HPC/MPI layer (implementation deferred)" OFF)
option(NEXTDFTB_ENABLE_PYTHON "Build the pybind11 Python module"                   ON )
option(NEXTDFTB_BUILD_CLI     "Build the standalone nextdftb CLI executable"      ON )
option(NEXTDFTB_BUILD_TESTS   "Build the CTest test suite"                        ON )
option(NEXTDFTB_BUILD_FORTRAN_MAIN "Build the native Fortran driver (dev tool)"   OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(compiler_flags)

# ---------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------
find_package(LAPACK REQUIRED)

if(NEXTDFTB_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

if(NEXTDFTB_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
endif()

if(NEXTDFTB_ENABLE_PYTHON)
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 REQUIRED CONFIG)
endif()

# ---------------------------------------------------------------------
# Subdirectories. Order matters: downstream consumers use nextdftb::*
# aliases that must already exist when referenced.
# ---------------------------------------------------------------------
add_subdirectory(abi)                 # C ABI headers (INTERFACE target)
add_subdirectory(kernels_fortran)     # Fortran compute engine
add_subdirectory(io)                  # placeholder I/O abstraction
add_subdirectory(hpc_mpi)             # deferred MPI interfaces
add_subdirectory(core_cpp)            # C++ orchestrator

if(NEXTDFTB_ENABLE_PYTHON)
    add_subdirectory(bindings)
endif()

if(NEXTDFTB_BUILD_CLI)
    add_subdirectory(cli)
endif()

if(NEXTDFTB_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------
# Configuration summary
# ---------------------------------------------------------------------
message(STATUS "")
message(STATUS "nextdftb ${PROJECT_VERSION} — configuration summary")
message(STATUS "  Build type  : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C           : ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "  C++         : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Fortran     : ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  OpenMP      : ${NEXTDFTB_ENABLE_OPENMP}")
message(STATUS "  LAPACK      : ${LAPACK_LIBRARIES}")
message(STATUS "  MPI         : ${NEXTDFTB_ENABLE_MPI} (design only)")
message(STATUS "  Python      : ${NEXTDFTB_ENABLE_PYTHON}")
message(STATUS "  CLI         : ${NEXTDFTB_BUILD_CLI}")
message(STATUS "  Tests       : ${NEXTDFTB_BUILD_TESTS}")
message(STATUS "  F90 driver  : ${NEXTDFTB_BUILD_FORTRAN_MAIN}")
message(STATUS "")
