!> Gestion centralisée des erreurs nextDFTB.
!>
!> Aucun STOP brutal sans message explicite. Tous les modules doivent
!> utiliser fatal/warn/info de ce module.
module errors
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    implicit none
    private

    public :: fatal, warn, info, assert

contains

    subroutine fatal(where, msg)
        character(len=*), intent(in) :: where, msg
        write(error_unit, '(a)') ""
        write(error_unit, '(a,a,a,a)') "[FATAL] ", trim(where), ": ", trim(msg)
        write(error_unit, '(a)') ""
        error stop 1
    end subroutine fatal

    subroutine warn(where, msg)
        character(len=*), intent(in) :: where, msg
        write(error_unit, '(a,a,a,a)') "[WARN ] ", trim(where), ": ", trim(msg)
    end subroutine warn

    subroutine info(where, msg)
        character(len=*), intent(in) :: where, msg
        write(output_unit, '(a,a,a,a)') "[INFO ] ", trim(where), ": ", trim(msg)
    end subroutine info

    subroutine assert(cond, where, msg)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: where, msg
        if (.not. cond) call fatal(where, msg)
    end subroutine assert
end module errors
