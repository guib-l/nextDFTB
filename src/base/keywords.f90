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

    ! GEOMETRY
    character(len=*), parameter, public :: KW_NATOMS = "NATOMS"

    ! BASIS
    character(len=*), parameter, public :: KW_SRC  = "SRC"
    character(len=*), parameter, public :: KW_EXT  = "EXT"
    character(len=*), parameter, public :: KW_SEP  = "SEP"
    character(len=*), parameter, public :: KW_TYPE = "TYPE"

    ! CALC > DFTB
    character(len=*), parameter, public :: KW_DFTB    = "DFTB"
    character(len=*), parameter, public :: KW_DFT     = "DFT"
    character(len=*), parameter, public :: KW_SKF     = "SKF"
    character(len=*), parameter, public :: KW_SCC     = "SCC"
    character(len=*), parameter, public :: KW_MAXSCC  = "MAXSCC"
    character(len=*), parameter, public :: KW_TOLSCC  = "TOLSCC"

    ! OUTPUT
    character(len=*), parameter, public :: KW_OUT = "OUT"
    character(len=*), parameter, public :: KW_LOG = "LOG"
end module keywords
