!> Output généraliste — bannières, sections, lignes formatées.
module write_output
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    private

    public :: open_output, close_output, banner, section, line

    integer, save :: out_unit = -1

contains

    subroutine open_output(filename)
        character(len=*), intent(in) :: filename
        integer :: ios
        open(newunit=out_unit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) out_unit = output_unit
    end subroutine open_output

    subroutine close_output()
        if (out_unit /= output_unit .and. out_unit >= 0) close(out_unit)
        out_unit = -1
    end subroutine close_output

    subroutine banner(title)
        character(len=*), intent(in) :: title
        if (out_unit < 0) return
        write(out_unit, '(a)') repeat('=', 60)
        write(out_unit, '(a)') trim(title)
        write(out_unit, '(a)') repeat('=', 60)
    end subroutine banner

    subroutine section(title)
        character(len=*), intent(in) :: title
        if (out_unit < 0) return
        write(out_unit, '(a)') ""
        write(out_unit, '(a)') "--- " // trim(title) // " ---"
    end subroutine section

    subroutine line(s)
        character(len=*), intent(in) :: s
        if (out_unit < 0) return
        write(out_unit, '(a)') trim(s)
    end subroutine line
end module write_output
