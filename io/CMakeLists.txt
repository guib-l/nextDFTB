# nextDFTB — I/O layer (placeholder).
#
# Only the text backend ships in the foundation phase. An HDF5 backend is
# planned (section 12 of the architecture document).

add_library(nextdftb_io STATIC src/text_writer.cpp)
add_library(nextdftb::io ALIAS nextdftb_io)

target_compile_features(nextdftb_io PUBLIC cxx_std_20)

target_include_directories(nextdftb_io
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(nextdftb_io PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(DIRECTORY include/ DESTINATION include)
