!> nextDFTB — physical and numerical constants (atomic units).
module constants
    use kinds, only: wp
    implicit none
    private

    real(wp), parameter, public :: PI     = 3.141592653589793_wp
    real(wp), parameter, public :: TWO_PI = 2.0_wp * PI

    real(wp), parameter, public :: BOHR_TO_ANGSTROM = 0.5291772109_wp
    real(wp), parameter, public :: ANGSTROM_TO_BOHR = 1.0_wp / BOHR_TO_ANGSTROM
    real(wp), parameter, public :: HARTREE_TO_EV    = 27.211386245988_wp
    real(wp), parameter, public :: EV_TO_HARTREE    = 1.0_wp / HARTREE_TO_EV
end module constants
