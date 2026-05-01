!> Mots-clés et objets de balise de l'input utilisateur.
!>
!> Ce module rassemble :
!>   - les constantes de chaînes des mots-clés autorisés ;
!>   - un objet par balise du format d'input, contenant ses options
!>     avec leurs valeurs par défaut directement dans la déclaration ;
!>   - l'objet agrégateur `input_kw_t` qui regroupe toutes les balises
!>     et qui est rempli par le parser.
module keywords
    use kinds,     only: wp
    use constants, only: PATH_LEN
    use property,  only: property_basis_t, property_dftb_t, property_dft_t
    implicit none
    private

    !-- Constantes : noms des balises ----------------------------------
    character(len=*), parameter, public :: KW_GEOMETRY = "GEOMETRY"
    character(len=*), parameter, public :: KW_BASIS    = "BASIS"
    character(len=*), parameter, public :: KW_CALC     = "CALC"
    character(len=*), parameter, public :: KW_DRIVER   = "DRIVER"
    character(len=*), parameter, public :: KW_OUTPUT   = "OUTPUT"
    character(len=*), parameter, public :: KW_SYSTEM   = "SYSTEM"
    character(len=*), parameter, public :: KW_OPTION   = "OPTION"

    !-- Constantes : noms des clés -------------------------------------
    ! GEOMETRY
    character(len=*), parameter, public :: KW_GEO = "GEO"
    ! BASIS
    character(len=*), parameter, public :: KW_SRC      = "SRC"
    character(len=*), parameter, public :: KW_EXT      = "EXT"
    character(len=*), parameter, public :: KW_SEP      = "SEP"
    character(len=*), parameter, public :: KW_TYPE     = "TYPE"
    character(len=*), parameter, public :: KW_ORBITALS = "ORBITALS"
    ! CALC
    character(len=*), parameter, public :: KW_DFTB   = "DFTB"
    character(len=*), parameter, public :: KW_DFT    = "DFT"
    character(len=*), parameter, public :: KW_SKF    = "SKF"
    character(len=*), parameter, public :: KW_SCC    = "SCC"
    character(len=*), parameter, public :: KW_MAXSCC = "MAXSCC"
    character(len=*), parameter, public :: KW_TOLSCC = "TOLSCC"
    ! DRIVER
    character(len=*), parameter, public :: KW_DRV_TYPE = "TYPE"
    ! OUTPUT
    character(len=*), parameter, public :: KW_OUT = "OUT"
    character(len=*), parameter, public :: KW_LOG = "LOG"

    !-- Objets de balise -----------------------------------------------

    !> Balise GEOMETRY.
    type, public :: geometry_kw_t
        character(len=PATH_LEN) :: geo = "geometry.dat"
    end type geometry_kw_t

    !> Balise CALC : sélecteur de méthode + sous-balises typées
    !> (chaque sous-balise est exactement l'objet `property_*` transmis
    !> au calculateur correspondant).
    type, public :: calc_kw_t
        character(len=8)       :: kind = "DFTB"     ! DFTB | DFT | SKF
        type(property_dftb_t)  :: dftb
        type(property_dft_t)   :: dft
        type(property_basis_t) :: skf
    end type calc_kw_t

    !> Balise DRIVER.
    type, public :: driver_kw_t
        character(len=16) :: type       = "SINGLE"
        logical           :: has_driver = .false.
    end type driver_kw_t

    !> Balise OUTPUT.
    type, public :: output_kw_t
        character(len=PATH_LEN) :: out    = "results.out"
        character(len=PATH_LEN) :: log    = "logs"
        logical                 :: log_on = .false.
    end type output_kw_t

    !> Balise SYSTEM (placeholder).
    type, public :: system_kw_t
        ! placeholder
    end type system_kw_t

    !> Balise OPTION (placeholder).
    type, public :: option_kw_t
        ! placeholder
    end type option_kw_t

    !> Agrégateur : rassemble toutes les balises de l'input.
    !> La balise BASIS est exactement un `property_basis_t`.
    type, public :: input_kw_t
        type(geometry_kw_t)    :: geometry
        type(property_basis_t) :: basis
        type(calc_kw_t)        :: calc
        type(driver_kw_t)      :: driver
        type(output_kw_t)      :: output
        type(system_kw_t)      :: system
        type(option_kw_t)      :: option
    end type input_kw_t
end module keywords
