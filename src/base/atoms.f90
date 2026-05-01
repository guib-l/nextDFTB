!> Objet Atoms : représentation d'un atome unique.
!>
!> Un atome est défini par son symbole chimique, sa position spatiale
!> (en bohr, unités atomiques), sa charge atomique initiale, sa masse
!> (en u.m.a., lue depuis le SKF homonucléaire), sa vitesse, ses
!> orbitales et leurs occupations, et son appartenance à une molécule
!> (par défaut 1). Les listes `orbitals` et `occupation` sont liées
!> 1:1 (même taille, même indexation).
module atoms_mod
    use kinds,     only: wp
    use constants, only: SYMBOL_LEN
    implicit none
    private

    integer, parameter, public :: ORB_LEN = 4

    type, public :: atoms_t
        character(len=SYMBOL_LEN)              :: symbol      = ""
        real(wp)                               :: position(3) = 0.0_wp
        real(wp)                               :: velocity(3) = 0.0_wp
        real(wp)                               :: charge      = 0.0_wp
        real(wp)                               :: mass        = 0.0_wp
        integer                                :: molecule_id = 1
        character(len=ORB_LEN), allocatable    :: orbitals(:)
        real(wp),               allocatable    :: occupation(:)
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
