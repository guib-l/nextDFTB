!> Tests : objet Structure (structure.f90).
program test_structure_prog
    use kinds,         only: wp
    use atoms_mod,     only: atoms_t, atoms_set
    use structure_mod, only: structure_t, structure_init, structure_set_atom
    use test_utils,    only: assert_eq_int, assert_eq_str, test_summary
    implicit none

    type(structure_t) :: s
    type(atoms_t)     :: a

    call structure_init(s, 3)
    call assert_eq_int(s%natoms, 3, "natoms")
    call assert_eq_int(s%nmolecules, 1, "nmolecules default 1")
    call assert_eq_int(s%molecules(1)%natoms, 0, "molecule empty before set")

    call atoms_set(a, "O", [0.0_wp, 0.0_wp, 0.0_wp], 0.0_wp, 1)
    call structure_set_atom(s, 1, a)
    call atoms_set(a, "H", [1.0_wp, 0.0_wp, 0.0_wp], 0.0_wp, 1)
    call structure_set_atom(s, 2, a)
    call atoms_set(a, "H", [0.0_wp, 1.0_wp, 0.0_wp], 0.0_wp, 1)
    call structure_set_atom(s, 3, a)

    call assert_eq_int(s%molecules(1)%natoms, 3, "molecule 1 has 3 atoms")
    call assert_eq_int(s%molecules(1)%atom_idx(1), 1, "idx 1")
    call assert_eq_int(s%molecules(1)%atom_idx(3), 3, "idx 3")
    call assert_eq_str(s%atoms(1)%symbol, "O", "atom 1 is O")
    call assert_eq_str(s%atoms(2)%symbol, "H", "atom 2 is H")

    call test_summary("structure")
end program test_structure_prog
