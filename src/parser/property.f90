!> Objets `property_*` : propriétés de calcul transmises aux calculateurs.
!>
!> Trois objets, un par méthode supportée :
!>   - `property_basis_t`  → calculateur SKF (et balise BASIS)
!>   - `property_dftb_t`   → calculateur DFTB
!>   - `property_dft_t`    → calculateur DFT (placeholder)
!>
!> Les valeurs par défaut sont définies dans la déclaration des types.
module property
    use kinds,     only: wp
    use constants, only: PATH_LEN, SYMBOL_LEN
    implicit none
    private

    !> Base abstraite pour les propriétés d'un calculateur de méthode
    !> (DFTB / DFT). Permet le passage polymorphe à `method_calc%init`.
    type, abstract, public :: property_method_t
    end type property_method_t

    !> Spécification d'orbitales pour un symbole atomique.
    type, public :: orbital_spec_t
        character(len=SYMBOL_LEN) :: symbol   = ""
        character(len=64)         :: orbitals = ""
    end type orbital_spec_t

    !> Propriétés du calculateur de base (SKF) / balise BASIS.
    type, public :: property_basis_t
        character(len=PATH_LEN)           :: src      = "."
        character(len=16)                 :: ext      = ".skf"
        character(len=4)                  :: sep      = "-"
        character(len=16)                 :: type     = "spd"
        type(orbital_spec_t), allocatable :: orbitals(:)
    end type property_basis_t

    !> Propriétés du calculateur DFTB.
    type, extends(property_method_t), public :: property_dftb_t
        logical  :: scc    = .false.
        integer  :: maxscc = 100
        real(wp) :: tolscc = 1.0e-5_wp
    end type property_dftb_t

    !> Propriétés du calculateur DFT (vide pour l'instant).
    type, extends(property_method_t), public :: property_dft_t
        ! placeholder
    end type property_dft_t
end module property
