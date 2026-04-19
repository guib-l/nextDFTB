# nextDFTB — standalone CLI executable.
#
# Links against core_cpp (which transitively pulls the ABI and Fortran
# kernels). This binary runs a full simulation with no Python involved.

add_executable(nextdftb_cli src/main.cpp)
add_executable(nextdftb::cli ALIAS nextdftb_cli)

target_link_libraries(nextdftb_cli PRIVATE nextdftb::core_cpp)

# The executable is installed under the plain name "nextdftb".
set_target_properties(nextdftb_cli PROPERTIES
    OUTPUT_NAME nextdftb
)

install(TARGETS nextdftb_cli RUNTIME DESTINATION bin)
