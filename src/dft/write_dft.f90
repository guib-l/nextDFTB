!> Output spécifique aux calculs DFT (placeholder).
!>
!> Le calculateur DFT n'est pas encore implémenté ; les méthodes
!> écrivent un message indiquant que la sortie n'est pas implémentée.
module write_dft
    use write_output, only: section, line
    use output_base,  only: output_base_t
    implicit none
    private

    type, extends(output_base_t), public :: output_dft_t
    contains
        procedure :: write_process => dft_write_process
        procedure :: write_result  => dft_write_result
    end type output_dft_t

contains

    subroutine dft_write_process(self)
        class(output_dft_t), intent(inout) :: self
        if (self%verbose >= 1) then
            call section("DFT process")
            call line("  (not implemented)")
        end if
    end subroutine dft_write_process

    subroutine dft_write_result(self)
        class(output_dft_t), intent(inout) :: self
        if (self%verbose >= 1) then
            call section("DFT result")
            call line("  (not implemented)")
        end if
    end subroutine dft_write_result
end module write_dft
