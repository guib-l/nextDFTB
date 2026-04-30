!> Résolution SCC-DFTB (et chemin non-SCC : un seul passage).
!>
!> Étapes :
!>   1. construire H0, S, gamma (une fois)
!>   2. Δq = 0 initial
!>   3. répéter :
!>        H_μν = H0_μν + (1/2) S_μν (V_A + V_B), V_A = Σ_C γ_AC Δq_C
!>        diagonaliser H C = S C E
!>        occupations Aufbau closed-shell
!>        P, q (Mulliken), Δq_new
!>        si max|Δq_new - Δq| < tol → fin
!>        mixer Δq
module scc
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use dftbstate,     only: dftbstate_t, basis_system_t
    use matel,         only: build_hs
    use gamma_mod,     only: build_gamma
    use density,       only: build_density
    use linalg,        only: solve_gen_eig
    use charges,       only: mulliken_charges, delta_charges
    use coulomb,       only: coulomb_potential
    use mixer,         only: mixer_t, init_mixer, free_mixer, mix_charges
    use write_dftb,    only: write_dftb_scc_header, write_dftb_scc_iter, &
                              write_dftb_scc_status
    implicit none
    private

    public :: solve_scc

contains

    subroutine solve_scc(struct, do_scc, maxit, tol, st)
        type(structure_t),  intent(in)    :: struct
        logical,            intent(in)    :: do_scc
        integer,            intent(in)    :: maxit
        real(wp),           intent(in)    :: tol
        type(dftbstate_t),  intent(inout) :: st

        integer  :: norb, natoms, it, mu, nu, ia, ib
        real(wp) :: max_diff, V_mu, V_nu
        real(wp), allocatable :: H0(:,:), V_atom(:), dq_new(:)
        type(mixer_t) :: mx

        norb   = st%bas%norb_total
        natoms = struct%natoms

        if (.not. allocated(st%H))     allocate(st%H(norb, norb))
        if (.not. allocated(st%S))     allocate(st%S(norb, norb))
        if (.not. allocated(st%C))     allocate(st%C(norb, norb))
        if (.not. allocated(st%eig))   allocate(st%eig(norb))
        if (.not. allocated(st%P))     allocate(st%P(norb, norb))
        if (.not. allocated(st%occ))   allocate(st%occ(norb))
        if (.not. allocated(st%q))     allocate(st%q(natoms))
        if (.not. allocated(st%dq))    allocate(st%dq(natoms))
        if (.not. allocated(st%gamma)) allocate(st%gamma(natoms, natoms))

        allocate(H0(norb, norb), V_atom(natoms), dq_new(natoms))

        call build_hs(struct, st%bas, H0, st%S)
        call build_gamma(struct, st%bas, st%gamma)

        st%dq = 0.0_wp
        call set_occ(st%bas%nelec, norb, st%occ)

        st%converged = .false.
        st%niter = 0

        call write_dftb_scc_header(do_scc, maxit, tol)

        do it = 1, maxit
            call coulomb_potential(st%gamma, st%dq, V_atom)

            st%H = H0
            if (do_scc) then
                do mu = 1, norb
                    ia = atom_of_orbital(st%bas, mu)
                    V_mu = V_atom(ia)
                    do nu = 1, norb
                        ib = atom_of_orbital(st%bas, nu)
                        V_nu = V_atom(ib)
                        st%H(mu, nu) = st%H(mu, nu) &
                            + 0.5_wp * st%S(mu, nu) * (V_mu + V_nu)
                    end do
                end do
            end if

            call solve_gen_eig(st%H, st%S, st%eig, st%C)

            call build_density(st%C, st%occ, st%P)
            call mulliken_charges(st%P, st%S, st%bas, st%q)
            call delta_charges(st%q, st%bas, dq_new)

            st%e_band = sum(st%occ * st%eig)

            if (.not. do_scc) then
                st%dq = dq_new
                st%niter = 1
                st%converged = .true.
                exit
            end if

            max_diff = maxval(abs(dq_new - st%dq))
            call write_dftb_scc_iter(it, max_diff)

            if (max_diff < tol) then
                st%dq = dq_new
                st%niter = it
                st%converged = .true.
                exit
            end if

            if (it == 1) call init_mixer(mx, natoms, 0)
            block
                real(wp) :: dq_mixed(natoms)
                call mix_charges(mx, st%dq, dq_new, dq_mixed)
                st%dq = dq_mixed
            end block
        end do

        if (.not. st%converged) st%niter = maxit

        if (do_scc) call write_dftb_scc_status(st%converged, st%niter)

        if (do_scc) call free_mixer(mx)
        deallocate(H0, V_atom, dq_new)
    end subroutine solve_scc


    pure function atom_of_orbital(bas, iorb) result(ia)
        type(basis_system_t), intent(in) :: bas
        integer,              intent(in) :: iorb
        integer :: ia, k
        ia = 1
        do k = 1, size(bas%atom_orb_start)
            if (iorb >= bas%atom_orb_start(k) .and. &
                iorb <  bas%atom_orb_start(k) + bas%atom_norb(k)) then
                ia = k; return
            end if
        end do
    end function atom_of_orbital


    !> Aufbau closed-shell : double occupation par orbitale.
    subroutine set_occ(nelec, norb, occ)
        integer,  intent(in)  :: nelec, norb
        real(wp), intent(out) :: occ(:)
        integer :: i, n_doubly, n_left
        occ = 0.0_wp
        n_doubly = nelec / 2
        n_left   = nelec - 2 * n_doubly
        do i = 1, min(n_doubly, norb)
            occ(i) = 2.0_wp
        end do
        if (n_left > 0 .and. n_doubly + 1 <= norb) occ(n_doubly + 1) = real(n_left, wp)
    end subroutine set_occ
end module scc
