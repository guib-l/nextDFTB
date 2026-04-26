# nextDFTB executable — assembled from the Fortran modules in src/.

add_executable(nextdftb
    kinds.f90
    constants.f90
    main.f90
)

target_link_libraries(nextdftb PRIVATE ${LAPACK_LIBRARIES})

if(NEXTDFTB_ENABLE_OPENMP)
    target_link_libraries(nextdftb PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(TARGETS nextdftb RUNTIME DESTINATION bin)
