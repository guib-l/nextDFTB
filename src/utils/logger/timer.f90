!> Mesure de temps d'exécution : registre nommé.
!>
!> API publique :
!>   - start_timer(name)   : démarre le chrono (création ou refresh)
!>   - stop_timer(name)    : accumule le delta et incrémente niter
!>   - write_timer()       : écrit le tableau récapitulatif sur l'output
!>   - timer_reset()       : vide le registre
module timer
    use kinds,        only: wp
    use write_output, only: output_unit_id, output_is_open
    implicit none
    private

    integer, parameter :: NAME_LEN = 32
    integer, parameter :: REG_MAX  = 64

    character(len=*), parameter :: TOTAL_NAME = "TOTAL"

    type :: timer_entry_t
        character(len=NAME_LEN) :: name     = ""
        real(wp)                :: dt_cumul = 0.0_wp
        real(wp)                :: t_start  = 0.0_wp
        integer                 :: niter    = 0
    end type timer_entry_t

    type(timer_entry_t), save :: g_reg(REG_MAX)
    integer,             save :: g_n = 0

    character(len=19), save :: g_start_stamp = ""
    character(len=19), save :: g_stop_stamp  = ""

    public :: start_timer, stop_timer, write_timer, timer_reset

contains

    subroutine start_timer(name)
        character(len=*), intent(in) :: name
        integer  :: idx
        real(wp) :: t

        idx = find_entry(name)
        if (idx == 0) idx = create_entry(name)
        if (idx == 0) return

        call cpu_time(t)
        g_reg(idx)%t_start = t

        if (trim(name) == TOTAL_NAME) g_start_stamp = now_stamp()
    end subroutine start_timer

    subroutine stop_timer(name)
        character(len=*), intent(in) :: name
        integer  :: idx
        real(wp) :: t

        idx = find_entry(name)
        if (idx == 0) return

        call cpu_time(t)
        g_reg(idx)%dt_cumul = g_reg(idx)%dt_cumul + (t - g_reg(idx)%t_start)
        g_reg(idx)%niter    = g_reg(idx)%niter + 1

        if (trim(name) == TOTAL_NAME) g_stop_stamp = now_stamp()
    end subroutine stop_timer

    subroutine write_timer()
        integer :: u, i, idx_total
        real(wp) :: dt_total, pct
        character(len=43) :: name_pad
        character(len=128) :: buf

        if (.not. output_is_open()) return
        if (g_n == 0) return

        u = output_unit_id()

        idx_total = find_entry(TOTAL_NAME)
        if (idx_total == 0) then
            dt_total = 0.0_wp
        else
            dt_total = g_reg(idx_total)%dt_cumul
        end if

        write(u, '(a)') repeat('*', 79)
        write(u, '(a)') " TIMER ::"
        write(u, '(a)') " "//hbar_top()
        write(u, '(a)') " |           Name                            |  Loop  |"// &
                        "  Time [s] | Relat [%] |"
        write(u, '(a)') " "//hbar_mid()

        do i = 1, g_n
            if (i == idx_total) cycle
            name_pad = pad_name(g_reg(i)%name)
            if (dt_total > 0.0_wp) then
                pct = 100.0_wp * g_reg(i)%dt_cumul / dt_total
            else
                pct = 0.0_wp
            end if
            write(buf, '(a,a,a,i6,a,f9.4,a,f8.3,a)') &
                " | ", name_pad, " |", g_reg(i)%niter, " |", &
                g_reg(i)%dt_cumul, " |", pct, "    |"
            write(u, '(a)') trim(buf)
        end do

        if (idx_total /= 0) then
            write(u, '(a)') " "//hbar_mid()
            name_pad = pad_name("TOTAL")
            write(buf, '(a,a,a,i6,a,f9.4,a,a)') &
                " | ", name_pad, " |", g_reg(idx_total)%niter, " |", &
                g_reg(idx_total)%dt_cumul, " |    100    |"
            write(u, '(a)') trim(buf)
        end if

        write(u, '(a)') " "//hbar_bot()
        write(u, '(a)') ""
        write(u, '(a,a)') " > Start at : ", g_start_stamp
        write(u, '(a,a)') " > Stop at  : ", g_stop_stamp
        write(u, '(a)') ""
        write(u, '(a)') repeat('*', 79)
    end subroutine write_timer

    subroutine timer_reset()
        integer :: i
        do i = 1, REG_MAX
            g_reg(i)%name     = ""
            g_reg(i)%dt_cumul = 0.0_wp
            g_reg(i)%t_start  = 0.0_wp
            g_reg(i)%niter    = 0
        end do
        g_n = 0
        g_start_stamp = ""
        g_stop_stamp  = ""
    end subroutine timer_reset

    !-- helpers privés -------------------------------------------------

    pure function find_entry(name) result(idx)
        character(len=*), intent(in) :: name
        integer :: idx, i
        idx = 0
        do i = 1, g_n
            if (trim(g_reg(i)%name) == trim(name)) then
                idx = i
                return
            end if
        end do
    end function find_entry

    function create_entry(name) result(idx)
        character(len=*), intent(in) :: name
        integer :: idx
        if (g_n >= REG_MAX) then
            idx = 0
            return
        end if
        g_n = g_n + 1
        idx = g_n
        g_reg(idx)%name     = name
        g_reg(idx)%dt_cumul = 0.0_wp
        g_reg(idx)%t_start  = 0.0_wp
        g_reg(idx)%niter    = 0
    end function create_entry

    function now_stamp() result(s)
        character(len=19) :: s
        character(len=8)  :: date
        character(len=10) :: time
        call date_and_time(date=date, time=time)
        s = date(1:4)//"-"//date(5:6)//"-"//date(7:8)//" "// &
            time(1:2)//":"//time(3:4)//":"//time(5:6)
    end function now_stamp

    pure function pad_name(s) result(o)
        character(len=*), intent(in) :: s
        character(len=43) :: o
        integer :: n
        o = repeat(' ', 43)
        n = min(len_trim(s), 43)
        if (n > 0) o(1:n) = s(1:n)
    end function pad_name

    pure function hbar_top() result(s)
        character(len=:), allocatable :: s
        s = char_box_top()
    end function hbar_top

    pure function hbar_mid() result(s)
        character(len=:), allocatable :: s
        s = char_box_mid()
    end function hbar_mid

    pure function hbar_bot() result(s)
        character(len=:), allocatable :: s
        s = char_box_bot()
    end function hbar_bot

    pure function char_box_top() result(s)
        character(len=:), allocatable :: s
        s = "┌"//repeat("─",43)//"┬"//repeat("─",8)//"┬"// &
            repeat("─",11)//"┬"//repeat("─",11)//"┐"
    end function char_box_top

    pure function char_box_mid() result(s)
        character(len=:), allocatable :: s
        s = "├"//repeat("─",43)//"┼"//repeat("─",8)//"┼"// &
            repeat("─",11)//"┼"//repeat("─",11)//"┤"
    end function char_box_mid

    pure function char_box_bot() result(s)
        character(len=:), allocatable :: s
        s = "└"//repeat("─",43)//"┴"//repeat("─",8)//"┴"// &
            repeat("─",11)//"┴"//repeat("─",11)//"┘"
    end function char_box_bot
end module timer
