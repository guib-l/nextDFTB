!> Objet Atoms : représentation d'un atome unique.
!>
!> Un atome est défini par son symbole chimique, sa position spatiale
!> (en bohr, unités atomiques), sa charge atomique initiale et son
!> appartenance à une molécule (par défaut 1).
module atoms_mod
    use kinds,     only: wp
    use constants, only: SYMBOL_LEN
    implicit none
    private

    type, public :: atoms_t
        character(len=SYMBOL_LEN) :: symbol      = ""
        real(wp)                  :: position(3) = 0.0_wp
        real(wp)                  :: charge      = 0.0_wp
        integer                   :: molecule_id = 1
    end type atoms_t

    public :: atoms_set
contains

    subroutine atoms_set(a, symbol, position, charge, molecule_id)
        type(atoms_t),    intent(out) :: a
        character(len=*), intent(in)  :: symbol
        real(wp),         intent(in)  :: position(3)
        real(wp),         intent(in)  :: charge
        integer,          intent(in)  :: molecule_id
        a%symbol      = symbol
        a%position    = position
        a%charge      = charge
        a%molecule_id = molecule_id
    end subroutine atoms_set
end module atoms_mod
