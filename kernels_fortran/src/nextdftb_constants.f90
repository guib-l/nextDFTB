!> nextDFTB — physical and numerical constants.
!>
!> Single source of truth for constants used across the compute layer.
!> All values are in SI-derived atomic units unless noted.
module nextdftb_constants
    use nextdftb_kinds, only: wp
    implicit none
    private

    ! --- Mathematical constants ----------------------------------------
    real(wp), parameter, public :: PI     = 3.141592653589793_wp
    real(wp), parameter, public :: TWO_PI = 2.0_wp * PI

    ! --- Unit conversions ----------------------------------------------
    real(wp), parameter, public :: BOHR_TO_ANGSTROM = 0.5291772109_wp
    real(wp), parameter, public :: ANGSTROM_TO_BOHR = 1.0_wp / BOHR_TO_ANGSTROM
    real(wp), parameter, public :: HARTREE_TO_EV    = 27.211386245988_wp
    real(wp), parameter, public :: EV_TO_HARTREE    = 1.0_wp / HARTREE_TO_EV

end module nextdftb_constants
