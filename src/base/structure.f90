!> Objet Structure : représentation complète du système moléculaire.
!>
!> Conteneur principal du calcul. Stocke la liste de tous les atomes du
!> système et la liste des molécules (par défaut une seule molécule
!> contenant tous les atomes).
module structure_mod
    use kinds,        only: wp
    use atoms_mod,    only: atoms_t
    use molecule_mod, only: molecule_t, molecule_init, molecule_add_atom
    implicit none
    private

    type, public :: structure_t
        integer                       :: natoms      = 0
        type(atoms_t),    allocatable :: atoms(:)            ! (natoms)
        integer                       :: nmolecules  = 0
        type(molecule_t), allocatable :: molecules(:)        ! (nmolecules)
    end type structure_t

    public :: structure_init, structure_set_atom

contains

    subroutine structure_init(s, natoms)
        type(structure_t), intent(out) :: s
        integer,           intent(in)  :: natoms
        s%natoms = natoms
        allocate(s%atoms(natoms))
        s%nmolecules = 1
        allocate(s%molecules(1))
        call molecule_init(s%molecules(1), 1)
    end subroutine structure_init

    !> Affecte un atome (1-based) et l'enregistre dans sa molécule.
    subroutine structure_set_atom(s, idx, a)
        type(structure_t), intent(inout) :: s
        integer,           intent(in)    :: idx
        type(atoms_t),     intent(in)    :: a
        integer :: mid
        s%atoms(idx) = a
        mid = a%molecule_id
        if (mid < 1 .or. mid > s%nmolecules) mid = 1
        call molecule_add_atom(s%molecules(mid), idx)
    end subroutine structure_set_atom
end module structure_mod
