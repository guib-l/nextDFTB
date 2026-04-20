# nextDFTB — Fortran 2018 production compute layer.
#
# Layout
#   core/   DFTB program itself: kinds, constants, runtime state, kernels,
#           and the optional native Fortran driver. Pure Fortran modules
#           with no iso_c_binding exposure.
#   abi/    Thin BIND(C) facade — translates DFTB kernels into the C ABI
#           symbols declared in abi/include/nextdftb/. No numerical work
#           happens here; the facade only marshals iso_c_binding types
#           and dispatches to core/.
#
# Module resolution is automatic (CMake's Fortran dependency tracking),
# but the logical build order is:
#   nextdftb_kinds          (core/)   precision / kind parameters
#   nextdftb_constants      (core/)   physical / numerical constants
#   nextdftb_runtime        (core/)   is_initialized flag
#   nextdftb_test_kernel    (core/)   infrastructure smoke-test routine
#   nextdftb_abi_errors     (abi/)    bind(C) declarations for nextdftb_set_error
#   nextdftb_abi_logging    (abi/)    bind(C) declarations for nextdftb_log
#   nextdftb_abi            (abi/)    BIND(C) facade exporting nextdftb_* symbols
#
# Linked downstream by nextdftb_core_cpp.

set(NEXTDFTB_FORTRAN_CORE_SOURCES
    core/nextdftb_kinds.f90
    core/nextdftb_constants.f90
    core/nextdftb_runtime.f90
    core/nextdftb_test.f90
)

set(NEXTDFTB_FORTRAN_ABI_SOURCES
    abi/nextdftb_abi_errors.f90
    abi/nextdftb_abi_logging.f90
    abi/nextdftb_abi.f90
)

add_library(nextdftb_engine STATIC
    ${NEXTDFTB_FORTRAN_CORE_SOURCES}
    ${NEXTDFTB_FORTRAN_ABI_SOURCES}
)
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
    add_executable(nextdftb_main core/nextdftb_main.f90)
    target_link_libraries(nextdftb_main PRIVATE nextdftb_engine)
    set_target_properties(nextdftb_main PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/engine
    )
endif()
