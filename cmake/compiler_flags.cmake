# Per-compiler Fortran flags for nextdftb.
# Applied via add_compile_options on the Fortran language only.

set(_nextdftb_flags_common    "")
set(_nextdftb_flags_debug     "")
set(_nextdftb_flags_release   "")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(_nextdftb_flags_common  -std=f2008 -Wall -Wextra -Wno-maybe-uninitialized)
    set(_nextdftb_flags_debug   -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    set(_nextdftb_flags_release -O3 -march=native -funroll-loops)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(_nextdftb_flags_common  -stand f08 -warn all)
    set(_nextdftb_flags_debug   -O0 -g -check all -traceback -fpe0)
    set(_nextdftb_flags_release -O3 -xHost)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(_nextdftb_flags_common  -Mstandard)
    set(_nextdftb_flags_debug   -O0 -g -Mbounds -traceback)
    set(_nextdftb_flags_release -O3 -fast)
endif()

add_compile_options(
    "$<$<COMPILE_LANGUAGE:Fortran>:${_nextdftb_flags_common}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:${_nextdftb_flags_debug}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:${_nextdftb_flags_release}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RelWithDebInfo>>:${_nextdftb_flags_release};-g>"
)
