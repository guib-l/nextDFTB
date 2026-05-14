!> Calcul de la matrice densité et de la matrice densité pondérée par
!> les énergies propres.
!>
!>   P_μν = Σ_i  occ_i  C_μi  C_νi
!>   W_μν = Σ_i  occ_i  eig_i  C_μi  C_νi
module density
    use kinds, only: wp
    implicit none
    private

    public :: build_density, build_dens_block, build_edens_block

contains

    !> Construit la matrice densité totale P et la matrice densité
    !> pondérée par les énergies propres W (norb × norb).
    subroutine build_density(C, occ, eig, P, W)
        real(wp), intent(in)  :: C(:,:)
        real(wp), intent(in)  :: occ(:)
        real(wp), intent(in)  :: eig(:)
        real(wp), intent(out) :: P(:,:)
        real(wp), intent(out) :: W(:,:)

        integer  :: norb, i, mu, nu
        real(wp) :: w_i

        norb = size(C, 1)
        P = 0.0_wp
        W = 0.0_wp
        do i = 1, norb
            if (occ(i) <= 0.0_wp) cycle
            w_i = occ(i) * eig(i)
            do nu = 1, norb
                do mu = 1, norb
                    P(mu, nu) = P(mu, nu) + occ(i) * C(mu, i) * C(nu, i)
                    W(mu, nu) = W(mu, nu) + w_i  * C(mu, i) * C(nu, i)
                end do
            end do
        end do
    end subroutine build_density


    !> Bloc densité P_{μν} = Σ_i occ_i C_{μi} C_{νi} restreint à
    !> μ ∈ [oa, oa+na-1], ν ∈ [ob, ob+nb-1].
    subroutine build_dens_block(C, occ, oa, na, ob, nb, P_block)
        real(wp), intent(in)  :: C(:,:)
        real(wp), intent(in)  :: occ(:)
        integer,  intent(in)  :: oa, na, ob, nb
        real(wp), intent(out) :: P_block(:,:)

        integer :: nstates, i, mu, nu

        nstates = size(C, 2)
        P_block = 0.0_wp
        do i = 1, nstates
            if (occ(i) <= 0.0_wp) cycle
            do nu = 1, nb
                do mu = 1, na
                    P_block(mu, nu) = P_block(mu, nu) &
                        + occ(i) * C(oa + mu - 1, i) * C(ob + nu - 1, i)
                end do
            end do
        end do
    end subroutine build_dens_block


    !> Bloc densité pondérée W_{μν} = Σ_i occ_i eig_i C_{μi} C_{νi}
    !> restreint à μ ∈ [oa, oa+na-1], ν ∈ [ob, ob+nb-1].
    subroutine build_edens_block(C, occ, eig, oa, na, ob, nb, W_block)
        real(wp), intent(in)  :: C(:,:)
        real(wp), intent(in)  :: occ(:)
        real(wp), intent(in)  :: eig(:)
        integer,  intent(in)  :: oa, na, ob, nb
        real(wp), intent(out) :: W_block(:,:)

        integer  :: nstates, i, mu, nu
        real(wp) :: w_i

        nstates = size(C, 2)
        W_block = 0.0_wp
        do i = 1, nstates
            if (occ(i) <= 0.0_wp) cycle
            w_i = occ(i) * eig(i)
            do nu = 1, nb
                do mu = 1, na
                    W_block(mu, nu) = W_block(mu, nu) &
                        + w_i * C(oa + mu - 1, i) * C(ob + nu - 1, i)
                end do
            end do
        end do
    end subroutine build_edens_block

end module density
