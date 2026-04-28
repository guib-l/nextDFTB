!> Diagonalisation généralisée H C = S C E (LAPACK dsygv).
module diag
    use kinds,  only: wp
    use errors, only: fatal
    implicit none
    private

    public :: gen_eig

contains

    subroutine gen_eig(H_in, S_in, eig, C)
        real(wp), intent(in)  :: H_in(:,:), S_in(:,:)
        real(wp), intent(out) :: eig(:)
        real(wp), intent(out) :: C(:,:)

        integer :: n, info, lwork
        real(wp), allocatable :: A(:,:), B(:,:), work(:)
        real(wp) :: wkopt(1)

        n = size(H_in, 1)
        allocate(A(n, n), B(n, n))
        A = H_in
        B = S_in

        ! workspace query
        lwork = -1
        call dsygv(1, 'V', 'U', n, A, n, B, n, eig, wkopt, lwork, info)
        if (info /= 0) call fatal("diag", "dsygv workspace query failed")
        lwork = nint(wkopt(1))
        allocate(work(lwork))

        call dsygv(1, 'V', 'U', n, A, n, B, n, eig, work, lwork, info)
        if (info /= 0) call fatal("diag", "dsygv failed")

        C = A
        deallocate(A, B, work)
    end subroutine gen_eig
end module diag
