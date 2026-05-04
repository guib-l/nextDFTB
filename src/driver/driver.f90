!> Driver par défaut : sélectionne le pipeline en fonction de l'input.
!>
!> Pour le calculateur DFTB/DFT, le driver par défaut est `single_point`.
!> Sinon, le driver termine simplement le programme.
module driver
    use kinds,           only: wp
    use structure_mod,   only: structure_t
    use keywords,        only: input_kw_t
    use single_point,    only: run_single_point
    use write_output,    only: open_output, close_output, &
                                write_header, write_footer
    use timer,           only: start_timer, stop_timer, write_timer, timer_reset
    use logger,          only: log_open, log_close
    use errors,          only: warn
    implicit none
    private

    character(len=*), parameter :: PROG_TITLE   = "nextDFTB — Calculation "
    character(len=*), parameter :: PROG_VERSION = "0.1.0"

    public :: run_default

contains

    subroutine run_default(struct, inp)
        type(structure_t), intent(in) :: struct
        type(input_kw_t),  intent(in) :: inp

        call timer_reset()
        call start_timer("TOTAL")

        if (inp%output%log_on) call log_open(trim(inp%output%log))
        call open_output(trim(inp%output%out))
        call write_header(PROG_TITLE, PROG_VERSION)

        select case (trim(inp%calc%kind))
        case ("DFTB", "DFT")
            call run_single_point(struct, inp)
        case default
            call warn("driver", "no driver path for calc kind: "//trim(inp%calc%kind))
        end select

        call stop_timer("TOTAL")
        call write_timer()
        call write_footer()

        call close_output()
        call log_close()
    end subroutine run_default
end module driver
