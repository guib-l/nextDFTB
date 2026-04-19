# nextDFTB — C ABI (header-only interface)
#
# Exposes the public C headers consumed by:
#   - Fortran kernels (via iso_c_binding BIND(C) declarations)
#   - C++ orchestrator (direct #include)
#   - pybind11 module (direct #include)
#   - CLI main.cpp (direct #include)

add_library(nextdftb_abi_c INTERFACE)
add_library(nextdftb::abi_c ALIAS nextdftb_abi_c)

target_include_directories(nextdftb_abi_c INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

install(DIRECTORY include/ DESTINATION include)
