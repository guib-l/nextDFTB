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

    subroutine fatal(where_, msg)
        character(len=*), intent(in) :: where_, msg
        write(error_unit, '(a)') ""
        write(error_unit, '(a,a,a,a)') "[FATAL] ", trim(where_), ": ", trim(msg)
        write(error_unit, '(a)') ""
        error stop 1
    end subroutine fatal

    subroutine warn(where_, msg)
        character(len=*), intent(in) :: where_, msg
        write(error_unit, '(a,a,a,a)') "[WARN ] ", trim(where_), ": ", trim(msg)
    end subroutine warn

    subroutine info(where_, msg)
        character(len=*), intent(in) :: where_, msg
        write(output_unit, '(a,a,a,a)') "[INFO ] ", trim(where_), ": ", trim(msg)
    end subroutine info

    subroutine assert(cond, where_, msg)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: where_, msg
        if (.not. cond) call fatal(where_, msg)
    end subroutine assert
end module errors
