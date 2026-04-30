!> Output spécifique au module SKF.
!>
!> Pour l'instant : résumé compact du magasin SKF (éléments chargés et
!> nombre de paires).
module write_skf
    use kinds,        only: wp
    use slakos,       only: skf_store_t
    use write_output, only: section, line
    implicit none
    private

    public :: write_skf_summary

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
end module write_skf
