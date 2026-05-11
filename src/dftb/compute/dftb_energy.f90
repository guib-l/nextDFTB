!> Calcul de l'énergie DFTB.
!>
!>   E_total = E_elec + E_scc + E_repulsive            (Elstner 1998, eq. 24)
!>   E_elec  = Σ_i occ_i Σ_{A,B} Σ_{μ,ν} c^i_μ c^i_ν H^0_{μν}
!>   E_scc   = (1/2) Σ_{A,B} γ_AB Δq_A Δq_B
!>
!> Pour le mode non-SCC, E_scc = 0.
module dftb_energy
    use kinds, only: wp
    implicit none
    private

    public :: band_energy, scc_energy, electronic_energy, total_energy

contains

    function band_energy(occ, eig) result(e)
        real(wp), intent(in) :: occ(:), eig(:)
        real(wp) :: e
        e = sum(occ * eig)
    end function band_energy

    function scc_energy(gamma, dq) result(e)
        real(wp), intent(in) :: gamma(:,:), dq(:)
        real(wp) :: e
        integer  :: i, j, n
        n = size(dq)
        e = 0.0_wp
        do i = 1, n
            do j = 1, n
                e = e + gamma(i, j) * dq(i) * dq(j)
            end do
        end do
        e = 0.5_wp * e
    end function scc_energy

    function electronic_energy(occ, C, H) result(e)
        real(wp), intent(in) :: occ(:)        ! (norb)
        real(wp), intent(in) :: C(:,:)        ! (norb, norb) coefficients MO
        real(wp), intent(in) :: H(:,:)       ! (norb, norb) Hamiltonien H^0
        real(wp) :: e, ei
        integer  :: i, mu, nu, norb
        norb = size(occ)
        e = 0.0_wp
        do i = 1, norb
            ei = 0.0_wp
            do nu = 1, norb
                do mu = 1, norb
                    ei = ei + C(mu, i) * C(nu, i) * H(mu, nu)
                end do
            end do
            e = e + occ(i) * ei
        end do
    end function electronic_energy

    function total_energy(e_elec, e_coul, e_rep) result(e)
        real(wp), intent(in) :: e_elec, e_coul, e_rep
        real(wp) :: e
        e = e_elec - e_coul + e_rep
    end function total_energy
end module dftb_energy
