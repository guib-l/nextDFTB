!> Mesure de temps d'exécution simple (CPU + wall).
module timer
    use kinds, only: wp
    implicit none
    private

    type, public :: timer_t
        real(wp) :: t_start = 0.0_wp
        real(wp) :: t_stop  = 0.0_wp
    end type timer_t

    public :: tic, toc, elapsed

contains

    subroutine tic(tm)
        type(timer_t), intent(out) :: tm
        call cpu_time(tm%t_start)
    end subroutine tic

    subroutine toc(tm)
        type(timer_t), intent(inout) :: tm
        call cpu_time(tm%t_stop)
    end subroutine toc

    function elapsed(tm) result(dt)
        type(timer_t), intent(in) :: tm
        real(wp) :: dt
        dt = tm%t_stop - tm%t_start
    end function elapsed
end module timer
