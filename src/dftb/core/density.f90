!> Calcul de la matrice densité.
!>
!>   P_μν = Σ_i  occ_i  C_μi  C_νi
module density
    use kinds, only: wp
    implicit none
    private

    public :: build_density

contains

    subroutine build_density(C, occ, P)
        real(wp), intent(in)  :: C(:,:)
        real(wp), intent(in)  :: occ(:)
        real(wp), intent(out) :: P(:,:)

        integer :: norb, i, mu, nu

        norb = size(C, 1)
        P = 0.0_wp
        do i = 1, norb
            if (occ(i) <= 0.0_wp) cycle
            do nu = 1, norb
                do mu = 1, norb
                    P(mu, nu) = P(mu, nu) + occ(i) * C(mu, i) * C(nu, i)
                end do
            end do
        end do
    end subroutine build_density
end module density
