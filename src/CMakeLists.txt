# nextDFTB executable — assemblé depuis les modules dans src/.

add_executable(nextdftb
    utils/kind.f90
    utils/constants.f90
    utils/units.f90
    base/atoms.f90
    base/molecule.f90
    base/neighbor.f90
    base/structure.f90
    globals.f90
    driver/abstract_calc.f90
    driver/abstract_basis.f90

    errors/errors.f90

    utils/logger/timer.f90
    utils/logger/logger.f90
    utils/output/write_output.f90
    utils/output/write_matrix.f90
    utils/output/output_base.f90

    parser/property.f90
    parser/keywords.f90
    parser/parse_input.f90
    parser/parse_geometry.f90

    skf/slakos.f90
    skf/readskf.f90
    skf/polint.f90
    skf/interp.f90
    skf/electronic.f90
    skf/repulsive.f90
    skf/write_skf.f90
    skf/skf.f90

    dftb/dftbstate.f90
    dftb/write_dftb.f90
    dftb/orbitals.f90
    dftb/core/skrot.f90
    dftb/core/matel.f90
    dftb/core/gamma.f90
    dftb/core/density.f90
    dftb/compute/coulomb.f90
    dftb/compute/charges.f90
    dftb/core/shift.f90
    dftb/compute/dftb_energy.f90
    dftb/compute/dftb_grad.f90
    dftb/solver/linalg.f90
    dftb/solver/mixer/mixer.f90
    dftb/solver/mixer/simple_mixer.f90
    dftb/solver/mixer/random_mixer.f90
    dftb/solver/mixer/broyden_mixer.f90
    dftb/solver/mixer/mixer_factory.f90
    dftb/solver/scc.f90
    dftb/repulsif/repulsif.f90
    dftb/disp/disp.f90
    dftb/pola/pola.f90
    dftb/dftb.f90

    dft/dft.f90
    dft/write_dft.f90

    driver/opt/opt.f90
    driver/md/md.f90
    driver/single_point.f90
    driver/driver.f90

    cli.f90
    main.f90
)

target_link_libraries(nextdftb PRIVATE ${LAPACK_LIBRARIES})

if(NEXTDFTB_ENABLE_OPENMP)
    target_link_libraries(nextdftb PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(TARGETS nextdftb RUNTIME DESTINATION bin)
