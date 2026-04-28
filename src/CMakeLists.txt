# nextDFTB executable — assemblé depuis les modules dans src/.

add_executable(nextdftb
    base/kind.f90
    base/constants.f90
    base/units.f90
    base/default.f90
    base/keywords.f90
    base/globals.f90

    errors/errors.f90

    io/parser/parse_input.f90
    io/parser/parse_param.f90
    io/logger/timer.f90
    io/logger/logger.f90
    io/output/write_matrix.f90
    io/output/write_output.f90
    io/output/dftb_results.f90

    dftb/skf/readskf.f90
    dftb/skf/interp.f90

    dftb/build/dftb_matrix.f90
    dftb/build/gamma.f90
    dftb/build/density.f90

    dftb/compute/diag.f90
    dftb/compute/charges.f90
    dftb/compute/coulomb.f90
    dftb/compute/dftb_energy.f90
    dftb/compute/dftb_grad.f90

    dftb/repulsif/repulsif.f90

    dftb/resolved/linalg.f90
    dftb/resolved/mixer.f90
    dftb/resolved/scc.f90

    dftb/qmmm/qmmm.f90
    dftb/disp/disp.f90
    dftb/pola/pola.f90

    skf/skf_gen.f90
    dft/dft.f90

    driver/opt/opt.f90
    driver/md/md.f90

    dftb/dftb.f90
    driver/driver.f90

    cli.f90
    main.f90
)

target_link_libraries(nextdftb PRIVATE ${LAPACK_LIBRARIES})

if(NEXTDFTB_ENABLE_OPENMP)
    target_link_libraries(nextdftb PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(TARGETS nextdftb RUNTIME DESTINATION bin)
