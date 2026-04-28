!> Logger simple : route les messages vers stdout et/ou un fichier.
!> Optionnel — n'a aucun effet sur le calcul.
module logger
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    private

    integer, save :: log_unit = -1
    logical, save :: log_active = .false.

    public :: log_open, log_close, log_msg, log_is_active

contains

    subroutine log_open(filename)
        character(len=*), intent(in) :: filename
        integer :: ios
        if (log_active) return
        open(newunit=log_unit, file=filename, status='replace', action='write', iostat=ios)
        if (ios == 0) log_active = .true.
    end subroutine log_open

    subroutine log_close()
        if (log_active) then
            close(log_unit)
            log_active = .false.
            log_unit   = -1
        end if
    end subroutine log_close

    subroutine log_msg(msg)
        character(len=*), intent(in) :: msg
        write(output_unit, '(a)') trim(msg)
        if (log_active) write(log_unit, '(a)') trim(msg)
    end subroutine log_msg

    function log_is_active() result(b)
        logical :: b
        b = log_active
    end function log_is_active
end module logger
