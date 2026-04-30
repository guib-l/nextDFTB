!> Mots-clés autorisés dans l'input utilisateur.
module keywords
    implicit none
    private

    ! Sections (balises)
    character(len=*), parameter, public :: KW_GEOMETRY = "GEOMETRY"
    character(len=*), parameter, public :: KW_BASIS    = "BASIS"
    character(len=*), parameter, public :: KW_CALC     = "CALC"
    character(len=*), parameter, public :: KW_DRIVER   = "DRIVER"
    character(len=*), parameter, public :: KW_OUTPUT   = "OUTPUT"
    character(len=*), parameter, public :: KW_SYSTEM   = "SYSTEM"
    character(len=*), parameter, public :: KW_OPTION   = "OPTION"

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
end module keywords
