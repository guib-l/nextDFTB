!> Conversions d'unités courantes.
module units
    use kinds, only: wp
    use constants, only: BOHR_TO_ANGSTROM, ANGSTROM_TO_BOHR, &
                         HARTREE_TO_EV, EV_TO_HARTREE
    implicit none
    private

    public :: ang_to_bohr, bohr_to_ang, ha_to_ev, ev_to_ha

contains

    elemental function ang_to_bohr(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y
        y = x * ANGSTROM_TO_BOHR
    end function ang_to_bohr

    elemental function bohr_to_ang(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y
        y = x * BOHR_TO_ANGSTROM
    end function bohr_to_ang

    elemental function ha_to_ev(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y
        y = x * HARTREE_TO_EV
    end function ha_to_ev

    elemental function ev_to_ha(x) result(y)
        real(wp), intent(in) :: x
        real(wp) :: y
        y = x * EV_TO_HARTREE
    end function ev_to_ha
end module units
