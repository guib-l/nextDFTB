!> Output spécifique aux calculs DFTB.
!>
!> Couvre :
!>   - en-tête de la boucle SCF
!>   - une ligne par itération SCF (it, max|Δq|)
!>   - statut de convergence
!>   - résultat final (énergies)
!>   - charges et populations (write_dftb_population)
module write_dftb
    use kinds,        only: wp
    use dftbstate,    only: dftbstate_t
    use write_output, only: section, subsection, line, &
                            kv_int, kv_log, kv_real_es, &
                            output_unit_id, output_is_open
    use write_matrix, only: print_matrix
    use output_base,  only: output_base_t
    implicit none
    private

    public :: write_dftb_scc_header, write_dftb_scc_iter
    public :: write_dftb_scc_status, write_dftb_final
    public :: write_dftb_matrices, write_dftb_population
    public :: write_dftb_gradient

    !> Objet de sortie DFTB. Stocke les résultats essentiels et délègue
    !> à `write_dftb_final` pour l'écriture finale.
    type, extends(output_base_t), public :: output_dftb_t
        real(wp)              :: e_total = 0.0_wp
        real(wp)              :: e_band  = 0.0_wp
        real(wp)              :: e_elec  = 0.0_wp
        real(wp)              :: e_scc   = 0.0_wp
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

    subroutine write_dftb_final(e_total, e_band, e_elec, e_scc, e_rep)
        real(wp), intent(in) :: e_total, e_band, e_elec, e_scc, e_rep

        call section("DFTB final result")

        call subsection("Energies (Hartree)")
        call kv_real_es("E_band",        e_band)
        call kv_real_es("E_H0",          e_elec)
        call kv_real_es("E_electronic",  e_elec+e_scc)
        call kv_real_es("E_scc",         e_scc)
        call kv_real_es("E_repulsive",   e_rep)
        call kv_real_es("E_total",       e_total)
    end subroutine write_dftb_final


    !> Affiche le gradient (3, natoms) en Hartree/Bohr.
    subroutine write_dftb_gradient(grad)
        real(wp), intent(in) :: grad(:,:)
        integer :: ia, natoms
        character(len=120) :: buf

        if (.not. output_is_open()) return

        natoms = size(grad, 2)
        call line("")
        call subsection("Gradient (Hartree/Bohr)")
        call line(" Atom         dE/dx              dE/dy              dE/dz")
        call line(" -----------------------------------------------------------------")
        do ia = 1, natoms
            write(buf, '(i5,3x,es18.8,1x,es18.8,1x,es18.8)') &
                ia, grad(1, ia), grad(2, ia), grad(3, ia)
            call line(trim(buf))
        end do
    end subroutine write_dftb_gradient


    !> Affiche les charges atomiques (dq) puis les populations par
    !> orbitale (m-shell) avec, en regard, la population de la l-shell
    !> et la population atomique brute (q).
    subroutine write_dftb_population(state)
        type(dftbstate_t), intent(in) :: state

        integer :: ia, ils, ks, mu_lo, mu_hi, mu, l_val, m_val
        integer :: o0, ls0, natoms
        real(wp) :: total_charge, lshell_val, atomic_val
        character(len=120) :: buf

        if (.not. output_is_open()) return

        natoms = size(state%bas%atom_norb)
        total_charge = 0.0_wp
        do ia = 1, natoms
            total_charge = total_charge + state%dq(ia)
        end do

        call line("")
        write(buf, '(a,f12.8)') " Total charge: ", total_charge
        call line(trim(buf))

        call line("")
        call line("> Atomic charges ")
        call line("")
        call line(" Atom   Population")
        call line(" ------------------")
        do ia = 1, natoms
            write(buf, '(i5,3x,f12.8)') ia, state%dq(ia)
            call line(trim(buf))
        end do

        call line("")
        call line("> Orbital populations ")
        call line("")
        call line(" Atom Sh.   l   m   Population  l-shell     Atomic")
        call line(" -----------------------------------------------------")

        do ia = 1, natoms
            o0  = state%bas%atom_orb_start(ia)
            ls0 = state%bas%atom_lshell_start(ia)
            atomic_val = state%q(ia)
            do ils = 1, state%bas%atom_nlshell(ia)
                call lshell_orb_range_loc(ils, mu_lo, mu_hi)
                lshell_val = state%lshell_q(ls0 + ils - 1)
                l_val = ils - 1
                ! Ordre m affiché : s (l=0) → 0 ; p (l=1) → 1,-1,0
                ! (px,py,pz) ; d (l=2) → -2,-1,0,1,2 (xy,yz,z²,xz,x²-y²).
                do ks = 1, mu_hi - mu_lo + 1
                    mu = o0 + (mu_lo - 1) + ks - 1
                    m_val = m_of_local_orbital(ils, ks)
                    if (ks == 1) then
                        if (ils == 1) then
                            write(buf, '(i5,i4,i4,i4,3x,f10.8,2x,f10.8,2x,f10.8)') &
                                ia, ils, l_val, m_val, state%mshell_q(mu), &
                                lshell_val, atomic_val
                        else
                            write(buf, '(i5,i4,i4,i4,3x,f10.8,2x,f10.8,2x,a)') &
                                ia, ils, l_val, m_val, state%mshell_q(mu), &
                                lshell_val, "---       "
                        end if
                    else
                        write(buf, '(i5,i4,i4,i4,3x,f10.8,2x,a,2x,a)') &
                            ia, ils, l_val, m_val, state%mshell_q(mu), &
                            "---       ", "---       "
                    end if
                    call line(trim(buf))
                end do
            end do
        end do
    end subroutine write_dftb_population


    subroutine write_dftb_matrices(H, S, P, C, G)
        real(wp), intent(in) :: H(:,:), S(:,:), P(:,:), C(:,:), G(:,:)
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

        open(newunit=u, file="gamma.dat", status='replace', &
             action='write', iostat=ios)
        if (ios == 0) then
            call print_matrix(u, "Gamma", G)
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
        call write_dftb_final(self%e_total, self%e_band, self%e_elec, &
                              self%e_scc, self%e_rep)
    end subroutine dftb_write_result


    !-- helpers privés -------------------------------------------------

    !> Plage locale d'orbitales de la sous-couche `ils` au sein d'un
    !> atome (1-based dans le bloc atomique). 1=s, 2..4=p, 5..9=d.
    pure subroutine lshell_orb_range_loc(ils, mu_lo, mu_hi)
        integer, intent(in)  :: ils
        integer, intent(out) :: mu_lo, mu_hi
        select case (ils)
        case (1); mu_lo = 1; mu_hi = 1
        case (2); mu_lo = 2; mu_hi = 4
        case (3); mu_lo = 5; mu_hi = 9
        case default; mu_lo = 1; mu_hi = 0
        end select
    end subroutine lshell_orb_range_loc


    !> Convention harmoniques sphériques réelles → m :
    !>   s : m=0
    !>   p : px ↔ +1, py ↔ -1, pz ↔ 0   (ordre matel : 1=px,2=py,3=pz)
    !>   d : dxy↔-2, dyz↔-1, dzx↔+1, dx²-y²↔+2, d3z²-r²↔0
    !>      (ordre matel : 1=xy,2=yz,3=zx,4=x²-y²,5=3z²-r²)
    pure function m_of_local_orbital(ils, k) result(m)
        integer, intent(in) :: ils, k
        integer :: m
        m = 0
        select case (ils)
        case (1)
            m = 0
        case (2)
            select case (k)
            case (1); m =  1
            case (2); m = -1
            case (3); m =  0
            end select
        case (3)
            select case (k)
            case (1); m = -2
            case (2); m = -1
            case (3); m =  1
            case (4); m =  2
            case (5); m =  0
            end select
        end select
    end function m_of_local_orbital
end module write_dftb
