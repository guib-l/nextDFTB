!> Output spécifique aux calculs DFTB.
module dftb_results
    use kinds,        only: wp
    use write_output, only: section, line
    implicit none
    private

    public :: write_dftb_summary

contains

    subroutine write_dftb_summary(e_total, e_band, e_coulomb, e_rep, charges)
        real(wp), intent(in) :: e_total, e_band, e_coulomb, e_rep
        real(wp), intent(in) :: charges(:)
        character(len=128) :: buf
        integer :: i

        call section("DFTB results (Hartree)")
        write(buf, '(a,es20.10)') "  E_band     = ", e_band;     call line(buf)
        write(buf, '(a,es20.10)') "  E_coulomb  = ", e_coulomb;  call line(buf)
        write(buf, '(a,es20.10)') "  E_repulsive= ", e_rep;      call line(buf)
        write(buf, '(a,es20.10)') "  E_total    = ", e_total;    call line(buf)

        call section("Mulliken charges (q_atom)")
        do i = 1, size(charges)
            write(buf, '(a,i5,a,f12.6)') "  atom ", i, " : ", charges(i)
            call line(buf)
        end do
    end subroutine write_dftb_summary
end module dftb_results
