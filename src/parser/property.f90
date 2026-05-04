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

    !> Paramètres du mixer de charges SCC.
    type, public :: property_mixer_t
        character(len=16) :: type    = "SIMPLE"   ! SIMPLE | RANDOM | BROYDEN
        real(wp)          :: factor  = 0.1_wp
        integer           :: history = 6
        real(wp)          :: omega0  = 0.01_wp
    end type property_mixer_t

    !> Propriétés du calculateur DFTB.
    type, extends(property_method_t), public :: property_dftb_t
        logical                 :: scc          = .false.
        integer                 :: maxscc       = 100
        real(wp)                :: tolscc       = 1.0e-5_wp
        logical                 :: write_matrix = .false.
        type(property_mixer_t)  :: mixing
    end type property_dftb_t

    !> Propriétés du calculateur DFT (vide pour l'instant).
    type, extends(property_method_t), public :: property_dft_t
        ! placeholder
    end type property_dft_t
end module property
