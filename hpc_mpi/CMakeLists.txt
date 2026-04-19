# nextDFTB — HPC/MPI layer (DESIGN ONLY, implementation deferred).
#
# When NEXTDFTB_ENABLE_MPI=OFF (default) we produce only a header-only
# INTERFACE target so that code can #include "nextdftb/hpc/..." without
# forcing an MPI dependency.
#
# When NEXTDFTB_ENABLE_MPI=ON the same target is promoted to a static
# library that compiles stubs.cpp — every entry point throws
# std::logic_error until the real implementation lands.

if(NEXTDFTB_ENABLE_MPI)
    add_library(nextdftb_hpc_mpi STATIC src/stubs.cpp)
    target_compile_features(nextdftb_hpc_mpi PUBLIC cxx_std_20)
    target_link_libraries(nextdftb_hpc_mpi PUBLIC MPI::MPI_CXX)
    set_target_properties(nextdftb_hpc_mpi PROPERTIES POSITION_INDEPENDENT_CODE ON)
else()
    add_library(nextdftb_hpc_mpi INTERFACE)
endif()

add_library(nextdftb::hpc_mpi ALIAS nextdftb_hpc_mpi)

if(TARGET nextdftb_hpc_mpi AND NOT NEXTDFTB_ENABLE_MPI)
    target_include_directories(nextdftb_hpc_mpi INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )
else()
    target_include_directories(nextdftb_hpc_mpi PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )
endif()

install(DIRECTORY include/ DESTINATION include)
