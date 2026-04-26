cmake_minimum_required(VERSION 3.20)

project(nextdftb
    VERSION     0.1.0
    DESCRIPTION "nextDFTB — Fortran 2018 DFTB compute engine"
    LANGUAGES   Fortran
)

# ---------------------------------------------------------------------
# Build type default
# ---------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# ---------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------
option(NEXTDFTB_ENABLE_OPENMP "Enable OpenMP in compute kernels" ON)

# ---------------------------------------------------------------------
# Fortran flags (enforce Fortran 2018; CMAKE_Fortran_STANDARD needs 3.28)
# ---------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-std=f2018>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
    )
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fcheck=all>)
    add_compile_options($<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-g>)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR
       CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-stand=f18>)
endif()

# ---------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------
find_package(LAPACK REQUIRED)

if(NEXTDFTB_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

# ---------------------------------------------------------------------
# Executable
# ---------------------------------------------------------------------
add_subdirectory(src)

# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------
message(STATUS "")
message(STATUS "nextdftb ${PROJECT_VERSION}")
message(STATUS "  Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fortran    : ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  OpenMP     : ${NEXTDFTB_ENABLE_OPENMP}")
message(STATUS "  LAPACK     : ${LAPACK_LIBRARIES}")
message(STATUS "")
