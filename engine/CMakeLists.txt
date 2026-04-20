# nextDFTB — Fortran 2018 production compute layer.
#
# Build order of modules (resolved automatically by the Fortran compiler via
# CMake's Fortran module dependency tracking):
#   nextdftb_kinds          — precision / kind parameters
#   nextdftb_constants      — physical / numerical constants
#   nextdftb_abi_errors     — bind(C) declarations for nextdftb_set_error
#   nextdftb_abi_logging    — bind(C) declarations for nextdftb_log
#   nextdftb_runtime        — is_initialized flag
#   nextdftb_test           — infrastructure smoke-test routine
#   nextdftb_abi            — BIND(C) facade exporting nextdftb_* symbols
#
# Linked downstream by nextdftb_core_cpp.

set(NEXTDFTB_FORTRAN_SOURCES
    src/nextdftb_kinds.f90
    src/nextdftb_constants.f90
    src/nextdftb_abi_errors.f90
    src/nextdftb_abi_logging.f90
    src/nextdftb_runtime.f90
    src/nextdftb_test.f90
    src/nextdftb_abi.f90
)

add_library(nextdftb_engine STATIC ${NEXTDFTB_FORTRAN_SOURCES})
add_library(nextdftb::engine ALIAS nextdftb_engine)

# Fortran modules emitted by this target end up in a dedicated sub-directory
# so downstream consumers (if any) can locate the .mod files.
set_target_properties(nextdftb_engine PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/engine
    POSITION_INDEPENDENT_CODE ON
)

target_link_libraries(nextdftb_engine
    PUBLIC
        nextdftb::abi_c
        LAPACK::LAPACK
)

if(NEXTDFTB_ENABLE_OPENMP)
    target_link_libraries(nextdftb_engine PUBLIC OpenMP::OpenMP_Fortran)
endif()

# ---------------------------------------------------------------------
# Optional native Fortran driver — bypasses the C ABI entirely.
# Intended for Fortran-only profiling and OpenMP scaling studies.
# ---------------------------------------------------------------------
if(NEXTDFTB_BUILD_FORTRAN_MAIN)
    add_executable(nextdftb_main src/nextdftb_main.f90)
    target_link_libraries(nextdftb_main PRIVATE nextdftb_engine)
    set_target_properties(nextdftb_main PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/engine
    )
endif()
