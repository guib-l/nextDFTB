!> Mesure de temps d'exécution simple + registre nommé.
!>
!> API :
!>   - tic / toc / elapsed     : timer local (type timer_t)
!>   - timer_record(name, tm)  : enregistre un timing dans le registre global
!>   - timer_count / timer_get : accès au registre (nom + durée)
!>   - timer_reset             : vide le registre
module timer
    use kinds, only: wp
    implicit none
    private

    type, public :: timer_t
        real(wp) :: t_start = 0.0_wp
        real(wp) :: t_stop  = 0.0_wp
    end type timer_t

    integer, parameter :: NAME_LEN  = 32
    integer, parameter :: REG_MAX   = 64

    type :: timer_entry_t
        character(len=NAME_LEN) :: name = ""
        real(wp) :: dt = 0.0_wp
    end type timer_entry_t

    type(timer_entry_t), save :: g_reg(REG_MAX)
    integer,             save :: g_n = 0

    public :: tic, toc, elapsed
    public :: timer_record, timer_count, timer_get, timer_reset

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

    subroutine timer_record(name, tm)
        character(len=*), intent(in) :: name
        type(timer_t),    intent(in) :: tm
        if (g_n >= REG_MAX) return
        g_n = g_n + 1
        g_reg(g_n)%name = name
        g_reg(g_n)%dt   = tm%t_stop - tm%t_start
    end subroutine timer_record

    function timer_count() result(n)
        integer :: n
        n = g_n
    end function timer_count

    subroutine timer_get(i, name, dt)
        integer,                intent(in)  :: i
        character(len=NAME_LEN), intent(out) :: name
        real(wp),               intent(out) :: dt
        if (i < 1 .or. i > g_n) then
            name = ""
            dt   = 0.0_wp
            return
        end if
        name = g_reg(i)%name
        dt   = g_reg(i)%dt
    end subroutine timer_get

    subroutine timer_reset()
        g_n = 0
    end subroutine timer_reset
end module timer
