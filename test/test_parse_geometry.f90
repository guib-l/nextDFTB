!> Tests : round-trip écriture/lecture geometry.dat.
program test_parse_geometry_prog
    use kinds,          only: wp
    use units,          only: ang_to_bohr
    use atoms_mod,      only: atoms_t, atoms_set
    use structure_mod,  only: structure_t, structure_init, structure_set_atom
    use parse_geometry, only: read_geometry, write_geometry
    use test_utils,     only: assert_close, assert_eq_int, assert_eq_str, test_summary
    implicit none

    real(wp), parameter :: TOL = 1.0e-6_wp
    character(len=*), parameter :: F = "tmp_test_geometry.dat"
    type(structure_t) :: s_out, s_in
    type(atoms_t)     :: a

    call structure_init(s_out, 2)
    call atoms_set(a, "O", &
        [ang_to_bohr(0.0_wp), ang_to_bohr(0.0_wp), ang_to_bohr(0.0_wp)], &
        0.0_wp, 1)
    call structure_set_atom(s_out, 1, a)
    call atoms_set(a, "H", &
        [ang_to_bohr(1.0_wp), ang_to_bohr(0.0_wp), ang_to_bohr(0.0_wp)], &
        0.1_wp, 1)
    call structure_set_atom(s_out, 2, a)

    call write_geometry(F, s_out, "test")
    call read_geometry(F, s_in)

    call assert_eq_int(s_in%natoms, 2, "natoms after read")
    call assert_eq_str(s_in%atoms(1)%symbol, "O", "atom 1 symbol")
    call assert_eq_str(s_in%atoms(2)%symbol, "H", "atom 2 symbol")
    call assert_close(s_in%atoms(1)%position(1), 0.0_wp, TOL, "atom 1 x")
    call assert_close(s_in%atoms(2)%position(1), ang_to_bohr(1.0_wp), TOL, "atom 2 x (bohr)")
    call assert_close(s_in%atoms(2)%charge, 0.1_wp, TOL, "atom 2 charge")

    open(unit=99, file=F, status='old')
    close(99, status='delete')

    call test_summary("parse_geometry")
end program test_parse_geometry_prog
