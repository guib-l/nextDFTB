!> Algèbre linéaire (diagonalisation) — délégation à diag.
module linalg
    use kinds, only: wp
    use diag,  only: gen_eig
    implicit none
    private

    public :: solve_gen_eig

contains

    subroutine solve_gen_eig(H, S, eig, C)
        real(wp), intent(in)  :: H(:,:), S(:,:)
        real(wp), intent(out) :: eig(:), C(:,:)
        call gen_eig(H, S, eig, C)
    end subroutine solve_gen_eig
end module linalg
