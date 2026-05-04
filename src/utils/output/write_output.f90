!> Output généraliste : entête, sections, lignes formatées, timings, footer.
!>
!> Structure d'un fichier d'output produit par nextDFTB :
!>
!>   1. write_header(title)        — bannière + date/heure
!>   2. section/subsection + kv_*  — propriétés du calcul, état
!>   3. (sortie spécifique au calculateur via write_dftb / write_skf)
!>   4. write_footer               — clôture
module write_output
    use, intrinsic :: iso_fortran_env, only: output_unit
    use kinds,  only: wp
    implicit none
    private

    integer, save :: out_unit = -1
    logical, save :: out_open = .false.

    public :: open_output, close_output, output_unit_id, output_is_open
    public :: banner, section, subsection, line
    public :: kv_str, kv_int, kv_real, kv_real_es, kv_log
    public :: write_header, write_footer

contains

    subroutine open_output(filename)
        character(len=*), intent(in) :: filename
        integer :: ios
        open(newunit=out_unit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) out_unit = output_unit
        out_open = .true.
    end subroutine open_output

    subroutine close_output()
        if (out_open .and. out_unit /= output_unit) close(out_unit)
        out_unit = -1
        out_open = .false.
    end subroutine close_output

    function output_unit_id() result(u)
        integer :: u
        u = out_unit
    end function output_unit_id

    function output_is_open() result(b)
        logical :: b
        b = out_open
    end function output_is_open

    !-- bannières -------------------------------------------------------

    subroutine banner(title)
        character(len=*), intent(in) :: title
        if (.not. out_open) return
        write(out_unit, '(a)') repeat('=', 80)
        write(out_unit, '(a)') trim(title)
        write(out_unit, '(a)') repeat('=', 80)
    end subroutine banner

    subroutine section(title)
        character(len=*), intent(in) :: title
        if (.not. out_open) return
        write(out_unit, '(a)') ""
        write(out_unit, '(a)') repeat('-', 80)
        write(out_unit, '(a)') " " // trim(title)
        write(out_unit, '(a)') repeat('-', 80)
    end subroutine section

    subroutine subsection(title)
        character(len=*), intent(in) :: title
        if (.not. out_open) return
        write(out_unit, '(a)') ""
        write(out_unit, '(a)') "  > " // trim(title)
    end subroutine subsection

    subroutine line(s)
        character(len=*), intent(in) :: s
        if (.not. out_open) return
        write(out_unit, '(a)') trim(s)
    end subroutine line

    !-- key / value (clés alignées sur 22 colonnes) ---------------------

    subroutine kv_str(key, val)
        character(len=*), intent(in) :: key, val
        if (.not. out_open) return
        write(out_unit, '(a,a,a)') "  ", pad22(key), trim(val)
    end subroutine kv_str

    subroutine kv_int(key, val)
        character(len=*), intent(in) :: key
        integer,          intent(in) :: val
        if (.not. out_open) return
        write(out_unit, '(a,a,i0)') "  ", pad22(key), val
    end subroutine kv_int

    subroutine kv_real(key, val)
        character(len=*), intent(in) :: key
        real(wp),         intent(in) :: val
        if (.not. out_open) return
        write(out_unit, '(a,a,f20.10)') "  ", pad22(key), val
    end subroutine kv_real

    subroutine kv_real_es(key, val)
        character(len=*), intent(in) :: key
        real(wp),         intent(in) :: val
        if (.not. out_open) return
        write(out_unit, '(a,a,es20.10)') "  ", pad22(key), val
    end subroutine kv_real_es

    subroutine kv_log(key, val)
        character(len=*), intent(in) :: key
        logical,          intent(in) :: val
        if (.not. out_open) return
        if (val) then
            write(out_unit, '(a,a,a)') "  ", pad22(key), "True"
        else
            write(out_unit, '(a,a,a)') "  ", pad22(key), "False"
        end if
    end subroutine kv_log

    !-- header / footer / timings ---------------------------------------

    subroutine write_header(title, version)
        character(len=*), intent(in) :: title, version
        character(len=8)  :: date
        character(len=10) :: time
        character(len=32) :: stamp
        if (.not. out_open) return
        call date_and_time(date=date, time=time)
        stamp = date(1:4)//"-"//date(5:6)//"-"//date(7:8)//"  "// &
                time(1:2)//":"//time(3:4)//":"//time(5:6)
        call banner(trim(title))
        call kv_str("version", version)
        call kv_str("started", trim(stamp))
    end subroutine write_header

    subroutine write_footer()
        character(len=8)  :: date
        character(len=10) :: time
        character(len=32) :: stamp
        if (.not. out_open) return
        call date_and_time(date=date, time=time)
        stamp = date(1:4)//"-"//date(5:6)//"-"//date(7:8)//"  "// &
                time(1:2)//":"//time(3:4)//":"//time(5:6)
        write(out_unit, '(a)') ""
        write(out_unit, '(a)') repeat('=', 80)
        write(out_unit, '(a,a)') " finished : ", trim(stamp)
        write(out_unit, '(a)') repeat('=', 80)
    end subroutine write_footer

    !-- helper interne --------------------------------------------------

    pure function pad22(s) result(o)
        character(len=*), intent(in) :: s
        character(len=22) :: o
        integer :: n
        o = repeat(' ', 22)
        n = min(len_trim(s), 19)
        o(1:n)        = s(1:n)
        o(20:22)      = " : "
    end function pad22
end module write_output
