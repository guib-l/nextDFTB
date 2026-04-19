# nextDFTB — per-compiler warning & optimisation flags.
#
# Applied via add_compile_options on the relevant language only, using
# generator expressions so each flag lands on the correct files.

# ---- Fortran 2018 -----------------------------------------------------
set(_nextdftb_f_common  "")
set(_nextdftb_f_debug   "")
set(_nextdftb_f_release "")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(_nextdftb_f_common  -std=f2018 -Wall -Wextra -Wno-maybe-uninitialized)
    set(_nextdftb_f_debug   -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    set(_nextdftb_f_release -O3 -march=native -funroll-loops)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(_nextdftb_f_common  -stand f18 -warn all)
    set(_nextdftb_f_debug   -O0 -g -check all -traceback -fpe0)
    set(_nextdftb_f_release -O3 -xHost)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(_nextdftb_f_common  -Mstandard)
    set(_nextdftb_f_debug   -O0 -g -Mbounds -traceback)
    set(_nextdftb_f_release -O3 -fast)
endif()

# ---- C (C11) ----------------------------------------------------------
set(_nextdftb_c_common  "")
set(_nextdftb_c_debug   "")
set(_nextdftb_c_release "")

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(_nextdftb_c_common  -Wall -Wextra -Wpedantic -Wshadow)
    set(_nextdftb_c_debug   -O0 -g -fno-omit-frame-pointer)
    set(_nextdftb_c_release -O3)
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(_nextdftb_c_common  -Wall -Wextra)
    set(_nextdftb_c_debug   -O0 -g)
    set(_nextdftb_c_release -O3 -xHost)
endif()

# ---- C++20 ------------------------------------------------------------
set(_nextdftb_cxx_common  "")
set(_nextdftb_cxx_debug   "")
set(_nextdftb_cxx_release "")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(_nextdftb_cxx_common  -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor)
    set(_nextdftb_cxx_debug   -O0 -g -fno-omit-frame-pointer)
    set(_nextdftb_cxx_release -O3)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(_nextdftb_cxx_common  -Wall -Wextra)
    set(_nextdftb_cxx_debug   -O0 -g)
    set(_nextdftb_cxx_release -O3 -xHost)
endif()

add_compile_options(
    "$<$<COMPILE_LANGUAGE:Fortran>:${_nextdftb_f_common}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:${_nextdftb_f_debug}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:${_nextdftb_f_release}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RelWithDebInfo>>:${_nextdftb_f_release};-g>"

    "$<$<COMPILE_LANGUAGE:C>:${_nextdftb_c_common}>"
    "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:${_nextdftb_c_debug}>"
    "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:${_nextdftb_c_release}>"
    "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:RelWithDebInfo>>:${_nextdftb_c_release};-g>"

    "$<$<COMPILE_LANGUAGE:CXX>:${_nextdftb_cxx_common}>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:${_nextdftb_cxx_debug}>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:${_nextdftb_cxx_release}>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RelWithDebInfo>>:${_nextdftb_cxx_release};-g>"
)
