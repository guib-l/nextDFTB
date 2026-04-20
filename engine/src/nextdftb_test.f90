!> nextDFTB — infrastructure smoke test.
!>
!> Single deterministic routine used to exercise the full stack
!> (C++ allocation -> C ABI -> Fortran + OpenMP -> ABI status return).
!> No physics, no caller-owned memory. Returns the closed-form value
!> of the sum 1 + 2 + ... + N (= N*(N+1)/2) for N = 1000 = 500500.
module nextdftb_test_kernel
    use nextdftb_kinds, only: wp, ip
    implicit none
    private

    public :: test_compute

contains

    subroutine test_compute(out)
        real(wp), intent(out) :: out

        integer(ip), parameter :: N = 1000_ip
        integer(ip) :: i
        real(wp)    :: s

        s = 0.0_wp
        !$omp parallel do default(none) private(i) &
        !$omp     reduction(+:s) schedule(static)
        do i = 1_ip, N
            s = s + real(i, wp)
        end do
        !$omp end parallel do

        out = s
    end subroutine test_compute

end module nextdftb_test_kernel
