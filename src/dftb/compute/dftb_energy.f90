!> Calcul de l'énergie DFTB.
!>
!>   E_total = E_band + E_coulomb + E_repulsive
!>
!> où E_band = sum_i n_i ε_i, ε_i étant les valeurs propres du H SCC.
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
        ! Note: pour SCC-DFTB, l'énergie de bande contient déjà le terme
        ! (1/2) Σ S P γ Δq, on retire la moitié de E_coul pour éviter le
        ! double comptage : E = E_band - E_coul + E_rep n'est PAS la forme
        ! standard. Forme retenue (Elstner 1998, eq. 24) :
        !   E = sum_i n_i ε_i  -  (1/2) Σ_{AB} γ_AB Δq_A Δq_B  +  E_rep
        e = e_band - e_coul + e_rep
    end function total_energy
end module dftb_energy
