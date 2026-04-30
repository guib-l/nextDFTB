!> Variables par défaut pour le parsing et le calcul.
module defaults
    use kinds, only: wp
    implicit none
    private

    character(len=*), parameter, public :: DEFAULT_GEO  = "geometry.dat"
    character(len=*), parameter, public :: DEFAULT_EXT  = ".skf"
    character(len=*), parameter, public :: DEFAULT_SEP  = "-"
    character(len=*), parameter, public :: DEFAULT_SRC  = "."
    character(len=*), parameter, public :: DEFAULT_TYPE = "spd"
    character(len=*), parameter, public :: DEFAULT_OUT  = "results.out"
    character(len=*), parameter, public :: DEFAULT_LOG  = "logs"

    integer,  parameter, public :: DEFAULT_MAXSCC = 100
    real(wp), parameter, public :: DEFAULT_TOLSCC = 1.0e-5_wp
    logical,  parameter, public :: DEFAULT_SCC    = .false.
    logical,  parameter, public :: DEFAULT_LOG_ON = .false.

    real(wp), parameter, public :: DEFAULT_MIX_ALPHA = 0.2_wp
    integer,  parameter, public :: DEFAULT_BROYDEN_HISTORY = 6
end module defaults
