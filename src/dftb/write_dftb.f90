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
    use write_matrix, only: print_matrix
    use output_base,  only: output_base_t
    implicit none
    private

    public :: write_dftb_scc_header, write_dftb_scc_iter
    public :: write_dftb_scc_status, write_dftb_final
    public :: write_dftb_matrices

    !> Objet de sortie DFTB. Stocke les résultats essentiels et délègue
    !> à `write_dftb_final` pour l'écriture finale.
    type, extends(output_base_t), public :: output_dftb_t
        real(wp)              :: e_total = 0.0_wp
        real(wp)              :: e_band  = 0.0_wp
        real(wp)              :: e_coul  = 0.0_wp
        real(wp)              :: e_rep   = 0.0_wp
        real(wp), allocatable :: charges(:)
    contains
        procedure :: write_process => dftb_write_process
        procedure :: write_result  => dftb_write_result
    end type output_dftb_t

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
        call line("    iter         dE_elec         max|dq|")
        call line("  --------  ----------------  ----------------")
    end subroutine write_dftb_scc_header

    subroutine write_dftb_scc_iter(it, dE_elec, max_diff, has_dE)
        integer,  intent(in) :: it
        real(wp), intent(in) :: dE_elec, max_diff
        logical,  intent(in) :: has_dE
        if (.not. output_is_open()) return
        if (has_dE) then
            write(output_unit_id(), '(a,i6,a,es18.6,a,es18.6)') &
                "  ", it, "  ", dE_elec, "  ", max_diff
        else
            write(output_unit_id(), '(a,i6,a,a18,a,es18.6)') &
                "  ", it, "  ", "      ---         ", "  ", max_diff
        end if
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


    subroutine write_dftb_matrices(H, S, P, C)
        real(wp), intent(in) :: H(:,:), S(:,:), P(:,:), C(:,:)
        integer :: u, ios

        open(newunit=u, file="hamiltonian.dat", status='replace', &
             action='write', iostat=ios)
        if (ios == 0) then
            call print_matrix(u, "H0", H)
            close(u)
        end if

        open(newunit=u, file="overlaps.dat", status='replace', &
             action='write', iostat=ios)
        if (ios == 0) then
            call print_matrix(u, "S", S)
            close(u)
        end if

        open(newunit=u, file="density.dat", status='replace', &
             action='write', iostat=ios)
        if (ios == 0) then
            call print_matrix(u, "P", P)
            close(u)
        end if

        open(newunit=u, file="coeff.dat", status='replace', &
             action='write', iostat=ios)
        if (ios == 0) then
            call print_matrix(u, "C", C)
            close(u)
        end if
    end subroutine write_dftb_matrices


    !-- output_dftb_t (impl. de output_base_t) -------------------------

    subroutine dftb_write_process(self)
        class(output_dftb_t), intent(inout) :: self
        if (self%verbose >= 2) call section("DFTB process")
    end subroutine dftb_write_process

    subroutine dftb_write_result(self)
        class(output_dftb_t), intent(inout) :: self
        real(wp), allocatable :: q(:)
        if (allocated(self%charges)) then
            q = self%charges
        else
            allocate(q(0))
        end if
        call write_dftb_final(self%e_total, self%e_band, self%e_coul, &
                              self%e_rep, q)
    end subroutine dftb_write_result
end module write_dftb
