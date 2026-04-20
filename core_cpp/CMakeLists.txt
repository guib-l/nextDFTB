# nextDFTB — C++20 orchestrator.
#
# Produces the central shared library that:
#   * Implements the error-state ABI (thread-local)
#   * Implements the logging-sink ABI (thread-safe)
#   * Provides the C++ API used by the CLI, pybind11 bindings and tests
#
# Links against nextdftb_abi_c (headers) and nextdftb_engine
# (the Fortran symbols nextdftb_init / _finalize / _test).

set(NEXTDFTB_CORE_CPP_SOURCES
    src/buffer.cpp
    src/error.cpp
    src/logger.cpp
    src/orchestrator.cpp
    src/version.cpp
)

add_library(nextdftb_core_cpp SHARED ${NEXTDFTB_CORE_CPP_SOURCES})
add_library(nextdftb::core_cpp ALIAS nextdftb_core_cpp)

target_compile_features(nextdftb_core_cpp PUBLIC cxx_std_20)

target_include_directories(nextdftb_core_cpp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_definitions(nextdftb_core_cpp PRIVATE
    NEXTDFTB_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    NEXTDFTB_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    NEXTDFTB_VERSION_PATCH=${PROJECT_VERSION_PATCH}
    NEXTDFTB_VERSION_STRING="${PROJECT_VERSION}"
)

target_link_libraries(nextdftb_core_cpp
    PUBLIC  nextdftb::abi_c
    PRIVATE nextdftb::engine
)

# When linking a C++ shared library against a Fortran static archive, we need
# the Fortran runtime (libgfortran / libifcore / etc.) to resolve Fortran
# intrinsics and OpenMP-Fortran runtime. CMake's Fortran language component
# of the OpenMP target already pulls what it needs; we additionally ensure
# the Fortran linker flags propagate.
set_target_properties(nextdftb_core_cpp PROPERTIES
    LINKER_LANGUAGE CXX
    POSITION_INDEPENDENT_CODE ON
)

install(TARGETS nextdftb_core_cpp
        EXPORT  nextdftbTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
