!> Constantes physiques et numériques (unités atomiques).
module constants
    use kinds, only: wp
    implicit none
    private

    real(wp), parameter, public :: PI     = 3.141592653589793_wp
    real(wp), parameter, public :: TWO_PI = 2.0_wp * PI
    real(wp), parameter, public :: SQRT_PI = 1.7724538509055159_wp

    real(wp), parameter, public :: BOHR_TO_ANGSTROM = 0.5291772109_wp
    real(wp), parameter, public :: ANGSTROM_TO_BOHR = 1.0_wp / BOHR_TO_ANGSTROM
    real(wp), parameter, public :: HARTREE_TO_EV    = 27.211386245988_wp
    real(wp), parameter, public :: EV_TO_HARTREE    = 1.0_wp / HARTREE_TO_EV

    real(wp), parameter, public :: ZERO  = 0.0_wp
    real(wp), parameter, public :: ONE   = 1.0_wp
    real(wp), parameter, public :: TWO   = 2.0_wp
    real(wp), parameter, public :: HALF  = 0.5_wp
end module constants
