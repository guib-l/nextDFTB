!> Calcul de l'énergie DFTB.
!>
!>   E_total = E_elec + E_repulsive                       (Elstner 1998, eq. 24)
!>   E_elec  = E_H0 + E_scc
!>   E_H0    = Σ_i occ_i Σ_{A,B} Σ_{μ,ν} c^i_μ c^i_ν H^0_{μν}
!>   E_scc   = (1/2) Σ_{A,B} γ_AB Δq_A Δq_B
!>
!> Pour le mode BASIC, E_scc = 0 et donc E_elec = E_H0.
module dftb_energy
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use dftbstate,     only: dftbstate_t
    use repulsif,      only: repulsive_energy
    use property,      only: SCHEME_BASIC, SCHEME_NOSCC, SCHEME_SCC
    implicit none
    private

    public :: band_energy, scc_energy, electronic_energy, total_energy, &
              compute_energy, scc_loop_energy, coulomb_energy

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

    function coulomb_energy(gamma, dq) result(e)
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
    end function coulomb_energy

    subroutine scc_loop_energy(st, ecit,escc)
        type(dftbstate_t), intent(inout) :: st
        real(wp), intent(out)  :: ecit,escc
        st%e_band = sum(st%occ * st%eig)
        ecit = coulomb_energy(st%gamma, st%dq)
        escc = st%e_band + ecit
    end subroutine scc_loop_energy

    
    function total_energy(e_elec, e_rep) result(e)
        real(wp), intent(in) :: e_elec, e_rep
        real(wp) :: e
        e = e_elec + e_rep
    end function total_energy

    !> Calcule toutes les énergies DFTB et les stocke dans `st`.
    !> Le schéma à utiliser est lu depuis `st%scheme` (positionné par
    !> `solve_scc` au début de la résolution).
    subroutine compute_energy(struct, st)
        type(structure_t), intent(in)    :: struct
        type(dftbstate_t), intent(inout) :: st

        st%e_band = band_energy(st%occ, st%eig)
        st%e_H0   = electronic_energy(st%occ, st%C, st%H0)

        select case (st%scheme)
        case (SCHEME_BASIC)
            st%e_scc = 0.0_wp
        case (SCHEME_NOSCC, SCHEME_SCC)
            st%e_scc = scc_energy(st%gamma, st%dq)
        case default
            st%e_scc = 0.0_wp
        end select

        st%e_coul = st%e_scc
        st%e_elec = st%e_H0 + st%e_scc
        st%e_rep  = repulsive_energy(struct)
        st%e_total = total_energy(st%e_elec, st%e_rep)
    end subroutine compute_energy

end module dftb_energy
