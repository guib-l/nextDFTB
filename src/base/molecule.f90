!> Objet Molecule : groupement d'atomes appartenant à la même molécule.
!>
!> Une molécule conserve la liste des indices (dans la structure parente)
!> des atomes qu'elle contient. Par défaut, lors de la lecture d'une
!> géométrie, tous les atomes appartiennent à la molécule numéro 1.
module molecule_mod
    use kinds, only: wp
    implicit none
    private

    type, public :: molecule_t
        integer              :: id        = 1
        integer              :: natoms    = 0
        integer, allocatable :: atom_idx(:)   ! indices des atomes (1-based)
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
end module molecule_mod
