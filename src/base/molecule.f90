!> Objet Molecule : groupement d'atomes appartenant à la même molécule.
!>
!> Une molécule conserve la liste des indices (dans la structure parente)
!> des atomes qu'elle contient et son centre de masse. Par défaut, lors
!> de la lecture d'une géométrie, tous les atomes appartiennent à la
!> molécule numéro 1.
module molecule_mod
    use kinds,     only: wp
    use atoms_mod, only: atoms_t
    implicit none
    private

    type, public :: molecule_t
        integer              :: id        = 1
        integer              :: natoms    = 0
        integer, allocatable :: atom_idx(:)        ! indices (1-based)
        real(wp)             :: com(3)    = 0.0_wp
    contains
        procedure :: set_com => molecule_set_com
    end type molecule_t

    public :: molecule_init, molecule_add_atom

contains

    subroutine molecule_init(m, id)
        type(molecule_t), intent(out) :: m
        integer,          intent(in)  :: id
        m%id     = id
        m%natoms = 0
        if (allocated(m%atom_idx)) deallocate(m%atom_idx)
        allocate(m%atom_idx(0))
        m%com = 0.0_wp
    end subroutine molecule_init

    subroutine molecule_add_atom(m, idx)
        type(molecule_t), intent(inout) :: m
        integer,          intent(in)    :: idx
        integer, allocatable :: tmp(:)
        allocate(tmp(m%natoms + 1))
        if (m%natoms > 0) tmp(1:m%natoms) = m%atom_idx(1:m%natoms)
        tmp(m%natoms + 1) = idx
        call move_alloc(tmp, m%atom_idx)
        m%natoms = m%natoms + 1
    end subroutine molecule_add_atom

    !> Calcule et stocke le centre de masse de la molécule.
    !> `atoms` est la liste complète des atomes de la structure parente ;
    !> seuls les indices `atom_idx` sont utilisés.
    subroutine molecule_set_com(self, atoms)
        class(molecule_t), intent(inout) :: self
        type(atoms_t),     intent(in)    :: atoms(:)
        real(wp) :: total_mass
        integer  :: k, i

        self%com   = 0.0_wp
        total_mass = 0.0_wp
        do k = 1, self%natoms
            i = self%atom_idx(k)
            self%com   = self%com + atoms(i)%mass * atoms(i)%position
            total_mass = total_mass + atoms(i)%mass
        end do
        if (total_mass > 0.0_wp) self%com = self%com / total_mass
    end subroutine molecule_set_com
end module molecule_mod
