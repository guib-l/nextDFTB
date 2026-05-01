!> Output spécifique au module SKF.
!>
!> Pour l'instant : résumé compact du magasin SKF (éléments chargés et
!> nombre de paires).
module write_skf
    use kinds,        only: wp
    use constants,    only: SYMBOL_LEN
    use slakos,       only: skf_store_t
    use write_output, only: section, line
    use output_base,  only: output_base_t
    implicit none
    private

    public :: write_skf_summary

    !> Objet de sortie SKF. Stocke la liste des symboles et le magasin
    !> à écrire. `write_result` délègue à `write_skf_summary`.
    type, extends(output_base_t), public :: output_skf_t
        character(len=SYMBOL_LEN), allocatable :: symbols(:)
        type(skf_store_t)                      :: store
    contains
        procedure :: write_process => skf_write_process
        procedure :: write_result  => skf_write_result
    end type output_skf_t

contains

    subroutine write_skf_summary(symbols, store)
        character(len=*),  intent(in) :: symbols(:)
        type(skf_store_t), intent(in) :: store
        character(len=128) :: buf
        integer :: i

        call section("SKF store")
        write(buf, '(a, i0)') "  elements loaded = ", store%nelem
        call line(buf)
        do i = 1, size(symbols)
            write(buf, '(a, a, a, f8.4)') "  ", trim(symbols(i)), &
                "   U_s = ", store%pair(i, i)%hubbard(3)
            call line(buf)
        end do
    end subroutine write_skf_summary


    !-- output_skf_t (impl. de output_base_t) --------------------------

    subroutine skf_write_process(self)
        class(output_skf_t), intent(inout) :: self
        if (self%verbose >= 2) call section("SKF process")
    end subroutine skf_write_process

    subroutine skf_write_result(self)
        class(output_skf_t), intent(inout) :: self
        if (.not. allocated(self%symbols)) return
        call write_skf_summary(self%symbols, self%store)
    end subroutine skf_write_result
end module write_skf
