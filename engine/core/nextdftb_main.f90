!> nextDFTB — standalone Fortran driver.
!>
!> Minimal scaffolding for the native Fortran entry point. Exposes two
!> virtual extension points (function_A, function_B) implemented as
!> procedure pointers bound to an abstract interface. Concrete physics
!> routines can replace the default stubs by reassigning the pointers.
program nextdftb_main
    use nextdftb_kinds, only: wp, ip
    !$ use omp_lib, only: omp_get_max_threads, omp_get_wtime

    implicit none

    abstract interface
        function virtual_func_iface(x, n) result(y)
            import :: wp, ip
            real(wp),    intent(in) :: x(*)
            integer(ip), intent(in) :: n
            real(wp)                :: y
        end function virtual_func_iface
    end interface

    procedure(virtual_func_iface), pointer :: function_A => null()
    procedure(virtual_func_iface), pointer :: function_B => null()

    integer(ip), parameter :: n_default = 1024_ip

    integer(ip)           :: n, i
    integer               :: argc, ios
    character(len=32)     :: arg
    real(wp), allocatable :: x(:)
    real(wp)              :: ra, rb, t0, t1

    ! --- Bind virtual slots to default stubs ---------------------------
    function_A => default_A
    function_B => default_B

    ! --- CLI parsing : `nextdftb_main [n]` -----------------------------
    n    = n_default
    argc = command_argument_count()
    if (argc >= 1) then
        call get_command_argument(1, arg)
        read(arg, *, iostat=ios) n
        if (ios /= 0 .or. n <= 0_ip) then
            write(*,'(A)') "nextdftb_main: fatal: invalid problem size"
            stop 1
        end if
    end if

    ! --- Banner ---------------------------------------------------------
    write(*,'(A)')          "nextDFTB - Fortran native driver"
    write(*,'(A, I0)')      "  problem size : ", n
    !$ write(*,'(A, I0)')   "  omp threads  : ", omp_get_max_threads()

    ! --- Allocation & initialization -----------------------------------
    allocate(x(n), stat=ios)
    if (ios /= 0) then
        write(*,'(A)') "nextdftb_main: fatal: allocation failed"
        stop 1
    end if

    !$omp parallel do default(none) shared(x, n) private(i) schedule(static)
    do i = 1_ip, n
        x(i) = real(i, wp)
    end do
    !$omp end parallel do

    ! --- Virtual pipeline ----------------------------------------------
    t0 = wall_time()
    ra = function_A(x, n)
    rb = function_B(x, n)
    t1 = wall_time()

    write(*,'(A, ES14.6)')  "  function_A   = ", ra
    write(*,'(A, ES14.6)')  "  function_B   = ", rb
    write(*,'(A, F9.3, A)') "  elapsed      = ", t1 - t0, " s"

    deallocate(x)

contains

    !> Default virtual stub for function_A. Replace by reassigning the
    !> module-level procedure pointer to a real implementation.
    function default_A(x, n) result(y)
        real(wp),    intent(in) :: x(*)
        integer(ip), intent(in) :: n
        real(wp)                :: y
        write(*,'(A)') "  [function_A] default stub (not implemented)"
        y = 0.0_wp
        if (n > 0_ip) y = x(1)
    end function default_A

    function default_B(x, n) result(y)
        real(wp),    intent(in) :: x(*)
        integer(ip), intent(in) :: n
        real(wp)                :: y
        write(*,'(A)') "  [function_B] default stub (not implemented)"
        y = 0.0_wp
        if (n > 0_ip) y = x(1)
    end function default_B

    !> Wall-clock seconds. Uses omp_get_wtime when available, falls back
    !> to system_clock otherwise.
    function wall_time() result(t)
        real(wp) :: t
        integer(ip) :: count, rate
        !$ t = real(omp_get_wtime(), wp)
        !$ return
        call system_clock(count, rate)
        t = real(count, wp) / real(max(rate, 1_ip), wp)
    end function wall_time

end program nextdftb_main
