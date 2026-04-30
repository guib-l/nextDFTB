!> Tests : conversions d'unités (units.f90).
program test_units_prog
    use kinds,      only: wp
    use units,      only: ang_to_bohr, bohr_to_ang, ha_to_ev, ev_to_ha
    use constants,  only: BOHR_TO_ANGSTROM, HARTREE_TO_EV
    use test_utils, only: assert_close, test_summary
    implicit none

    real(wp), parameter :: TOL = 1.0e-12_wp
    real(wp) :: x

    x = 1.0_wp
    call assert_close(bohr_to_ang(x), BOHR_TO_ANGSTROM, TOL, "bohr_to_ang(1)")
    call assert_close(ang_to_bohr(bohr_to_ang(x)), x, TOL, "ang/bohr round-trip")

    call assert_close(ha_to_ev(x), HARTREE_TO_EV, TOL, "ha_to_ev(1)")
    call assert_close(ev_to_ha(ha_to_ev(x)), x, TOL, "ha/ev round-trip")

    x = 2.5_wp
    call assert_close(bohr_to_ang(x), 2.5_wp * BOHR_TO_ANGSTROM, TOL, "bohr_to_ang(2.5)")

    call test_summary("units")
end program test_units_prog
