!> Driver par défaut : sélectionne le pipeline en fonction de l'input.
!>
!> Pour le calculateur DFTB/DFT, le driver par défaut est `single_point`.
!> Sinon, le driver termine simplement le programme.
module driver
    use kinds,           only: wp
    use structure_mod,   only: structure_t
    use parse_input,     only: input_t
    use single_point,    only: run_single_point
    use write_output,    only: open_output, close_output, &
                                write_header, write_footer, write_timings
    use timer,           only: timer_t, tic, toc, timer_record, timer_reset
    use logger,          only: log_open, log_close
    use errors,          only: warn
    implicit none
    private

    character(len=*), parameter :: PROG_TITLE   = "nextDFTB — single point"
    character(len=*), parameter :: PROG_VERSION = "0.1.0"

    public :: run_default

contains

    subroutine run_default(struct, inp)
        type(structure_t), intent(in) :: struct
        type(input_t),     intent(in) :: inp
        type(timer_t) :: t_total

        call timer_reset()
        call tic(t_total)

        if (inp%out%log_on) call log_open(trim(inp%out%log))
        call open_output(trim(inp%out%out))
        call write_header(PROG_TITLE, PROG_VERSION)

        select case (trim(inp%calc%kind))
        case ("DFTB", "DFT")
            call run_single_point(struct, inp)
        case default
            call warn("driver", "no driver path for calc kind: "//trim(inp%calc%kind))
        end select

        call toc(t_total)
        call timer_record("total", t_total)

        call write_timings()
        call write_footer()

        call close_output()
        call log_close()
    end subroutine run_default
end module driver
