!> Tests : objet Atoms (atoms.f90).
program test_atoms_prog
    use kinds,      only: wp
    use atoms_mod,  only: atoms_t, atoms_set
    use test_utils, only: assert_close, assert_eq_int, assert_eq_str, test_summary
    implicit none

    real(wp), parameter :: TOL = 1.0e-14_wp
    type(atoms_t) :: a

    ! valeurs par défaut
    call assert_eq_str(a%symbol, "", "default symbol empty")
    call assert_close(a%charge, 0.0_wp, TOL, "default charge 0")
    call assert_eq_int(a%molecule_id, 1, "default molecule_id 1")

    call atoms_set(a, "O", [0.1_wp, 0.2_wp, 0.3_wp], -0.5_wp, 2)
    call assert_eq_str(a%symbol, "O", "symbol set")
    call assert_close(a%position(1), 0.1_wp, TOL, "pos x")
    call assert_close(a%position(2), 0.2_wp, TOL, "pos y")
    call assert_close(a%position(3), 0.3_wp, TOL, "pos z")
    call assert_close(a%charge, -0.5_wp, TOL, "charge set")
    call assert_eq_int(a%molecule_id, 2, "molecule_id set")

    call test_summary("atoms")
end program test_atoms_prog
