!> Output spécifique aux calculs DFTB.
!>
!> Couvre :
!>   - en-tête de la boucle SCF
!>   - une ligne par itération SCF (it, max|Δq|)
!>   - statut de convergence
!>   - résultat final (énergies, charges Mulliken)
module write_dftb
    use kinds,        only: wp
    use write_output, only: section, subsection, line, &
                            kv_int, kv_log, kv_real, kv_real_es, &
                            output_unit_id, output_is_open
    implicit none
    private

    public :: write_dftb_scc_header, write_dftb_scc_iter
    public :: write_dftb_scc_status, write_dftb_final

contains

    subroutine write_dftb_scc_header(do_scc, maxscc, tolscc)
        logical,  intent(in) :: do_scc
        integer,  intent(in) :: maxscc
        real(wp), intent(in) :: tolscc

        call section("SCF cycle")
        call kv_log("scc_enabled", do_scc)
        if (.not. do_scc) return
        call kv_int("max_iterations", maxscc)
        call kv_real_es("tolerance", tolscc)
        call line("")
        call line("    iter        max|dq|")
        call line("  --------  --------------")
    end subroutine write_dftb_scc_header

    subroutine write_dftb_scc_iter(it, max_diff)
        integer,  intent(in) :: it
        real(wp), intent(in) :: max_diff
        if (.not. output_is_open()) return
        write(output_unit_id(), '(a,i6,a,es16.6)') "  ", it, "    ", max_diff
    end subroutine write_dftb_scc_iter

    subroutine write_dftb_scc_status(converged, niter)
        logical, intent(in) :: converged
        integer, intent(in) :: niter
        call line("")
        if (converged) then
            call kv_int("converged in", niter)
        else
            call kv_int("NOT converged after", niter)
        end if
    end subroutine write_dftb_scc_status

    subroutine write_dftb_final(e_total, e_band, e_coulomb, e_rep, charges_)
        real(wp), intent(in) :: e_total, e_band, e_coulomb, e_rep
        real(wp), intent(in) :: charges_(:)
        character(len=128) :: buf
        integer :: i

        call section("DFTB final result")

        call subsection("Energies (Hartree)")
        call kv_real_es("E_band",      e_band)
        call kv_real_es("E_coulomb",   e_coulomb)
        call kv_real_es("E_repulsive", e_rep)
        call kv_real_es("E_total",     e_total)

        call subsection("Mulliken charges")
        do i = 1, size(charges_)
            write(buf, '(a,i5,a,f12.6)') "  atom ", i, "  q = ", charges_(i)
            call line(buf)
        end do
    end subroutine write_dftb_final
end module write_dftb
