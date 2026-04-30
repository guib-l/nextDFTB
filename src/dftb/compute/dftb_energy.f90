!> Calcul de l'énergie DFTB.
!>
!>   E_total = E_band - E_coulomb + E_repulsive    (Elstner 1998, eq. 24)
!>
!> Pour le mode non-SCC, E_coulomb = 0.
module dftb_energy
    use kinds, only: wp
    implicit none
    private

    public :: band_energy, total_energy

contains

    function band_energy(occ, eig) result(e)
        real(wp), intent(in) :: occ(:), eig(:)
        real(wp) :: e
        e = sum(occ * eig)
    end function band_energy

    function total_energy(e_band, e_coul, e_rep) result(e)
        real(wp), intent(in) :: e_band, e_coul, e_rep
        real(wp) :: e
        e = e_band - e_coul + e_rep
    end function total_energy
end module dftb_energy
