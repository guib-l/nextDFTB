cmake_minimum_required(VERSION 3.20)

project(nextdftb
    VERSION 0.1.0
    DESCRIPTION "Next-generation DFTB code"
    LANGUAGES Fortran
)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

option(NEXTDFTB_ENABLE_MPI    "Enable MPI parallelism"    OFF)
option(NEXTDFTB_ENABLE_OPENMP "Enable OpenMP parallelism" OFF)
option(NEXTDFTB_BUILD_TESTS   "Enable CTest-driven tests" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(compiler_flags)

find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)
find_package(HDF5   REQUIRED COMPONENTS Fortran HL)

if(NEXTDFTB_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

if(NEXTDFTB_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

add_subdirectory(src)
# add_subdirectory(app)

if(NEXTDFTB_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

message(STATUS "nextdftb configuration summary")
message(STATUS "  Build type  : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fortran     : ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  BLAS libs   : ${BLAS_LIBRARIES}")
message(STATUS "  LAPACK libs : ${LAPACK_LIBRARIES}")
message(STATUS "  HDF5 Fortran: ${HDF5_Fortran_LIBRARIES}")
message(STATUS "  MPI         : ${NEXTDFTB_ENABLE_MPI}")
message(STATUS "  OpenMP      : ${NEXTDFTB_ENABLE_OPENMP}")
