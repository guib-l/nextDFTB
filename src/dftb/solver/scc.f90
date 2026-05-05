!> Résolution DFTB selon trois schémas distincts :
!>
!>   - BASIC : pas de terme dépendant des charges, une seule
!>             diagonalisation H = H0.
!>   - NOSCC : une diagonalisation avec V calculé à partir des charges
!>             fournies dans l'input de géométrie (pas d'auto-cohérence).
!>   - SCC   : cycle SCC complet jusqu'à convergence sur Δq.
!>
!> Les helpers `build_shift` et `build_ham` factorisent la construction
!> de l'hamiltonien shifté H_μν = H0_μν + (1/2) S_μν (V_A + V_B).
module scc
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use dftbstate,     only: dftbstate_t, basis_system_t
    use matel,         only: build_hs
    use gamma_mod,     only: build_gamma
    use density,       only: build_density
    use linalg,        only: solve_gen_eig
    use charges,       only: mulliken_charges, delta_charges
    use coulomb,       only: coulomb_potential, coulomb_energy
    use mixer,         only: mixer_t
    use mixer_factory, only: make_mixer
    use property,      only: property_mixer_t,                          &
                              SCHEME_BASIC, SCHEME_NOSCC, SCHEME_SCC
    use write_dftb,    only: write_dftb_scc_header, write_dftb_scc_iter, &
                              write_dftb_scc_status, write_dftb_matrices
    use errors,        only: fatal
    implicit none
    private

    public :: solve_scc

contains

    !> Dispatcher : exécute le schéma de calcul demandé.
    subroutine solve_scc(struct, scheme, gamma_kind, maxit, tol,        &
                         write_matrix, mixing, st)
        type(structure_t),      intent(in)    :: struct
        integer,                intent(in)    :: scheme
        integer,                intent(in)    :: gamma_kind
        integer,                intent(in)    :: maxit
        real(wp),               intent(in)    :: tol
        logical,                intent(in)    :: write_matrix
        type(property_mixer_t), intent(in)    :: mixing
        type(dftbstate_t),      intent(inout) :: st

        call ensure_allocated(struct, st, scheme)

        select case (scheme)
        case (SCHEME_BASIC)
            call solve_basic(struct, write_matrix, st)
        case (SCHEME_NOSCC)
            call solve_noscc(struct, gamma_kind, write_matrix, st)
        case (SCHEME_SCC)
            call solve_scc_loop(struct, gamma_kind, maxit, tol,         &
                                write_matrix, mixing, st)
        case default
            call fatal("scc", "unknown DFTB scheme")
        end select
    end subroutine solve_scc


    !-- Schéma BASIC ----------------------------------------------------

    subroutine solve_basic(struct, write_matrix, st)
        type(structure_t), intent(in)    :: struct
        logical,           intent(in)    :: write_matrix
        type(dftbstate_t), intent(inout) :: st

        call build_hs(struct, st%bas, st%H, st%S)
        call set_occ(st%bas%nelec, st%bas%norb_total, st%occ)

        call solve_gen_eig(st%H, st%S, st%eig, st%C)
        call build_density(st%C, st%occ, st%P)
        call mulliken_charges(st%P, st%S, st%bas, st%q)
        call delta_charges(st%q, st%bas, st%dq)

        if (write_matrix) call write_dftb_matrices(st%H, st%S, st%P, st%C)

        st%e_band    = sum(st%occ * st%eig)
        st%niter     = 1
        st%converged = .true.
    end subroutine solve_basic


    !-- Schéma NOSCC ----------------------------------------------------

    subroutine solve_noscc(struct, gamma_kind, write_matrix, st)
        type(structure_t), intent(in)    :: struct
        integer,           intent(in)    :: gamma_kind
        logical,           intent(in)    :: write_matrix
        type(dftbstate_t), intent(inout) :: st

        integer  :: ia, natoms
        real(wp), allocatable :: H0(:,:), S_shift(:,:), V_atom(:)

        natoms = struct%natoms
        allocate(H0(st%bas%norb_total, st%bas%norb_total))
        allocate(S_shift(st%bas%norb_total, st%bas%norb_total))
        allocate(V_atom(natoms))

        call build_hs(struct, st%bas, H0, st%S)
        call build_gamma(struct, st%bas, st%gamma, gamma_kind)

        ! Δq initial = - charge formelle fournie par l'input géométrie.
        ! Convention : charge>0 (cation) ⇒ moins d'électrons ⇒ Δq<0.
        do ia = 1, natoms
            st%dq(ia) = -struct%atoms(ia)%charge
        end do

        call coulomb_potential(st%gamma, st%dq, V_atom)
        call build_shift(st%S, V_atom, st%bas, S_shift)
        call build_ham(H0, S_shift, st%H)

        call set_occ(st%bas%nelec, st%bas%norb_total, st%occ)
        call solve_gen_eig(st%H, st%S, st%eig, st%C)
        call build_density(st%C, st%occ, st%P)
        call mulliken_charges(st%P, st%S, st%bas, st%q)

        if (write_matrix) call write_dftb_matrices(H0, st%S, st%P, st%C)

        st%e_band    = sum(st%occ * st%eig)
        st%niter     = 1
        st%converged = .true.

        deallocate(H0, S_shift, V_atom)
    end subroutine solve_noscc


    !-- Schéma SCC complet ---------------------------------------------

    subroutine solve_scc_loop(struct, gamma_kind, maxit, tol,           &
                              write_matrix, mixing, st)
        type(structure_t),      intent(in)    :: struct
        integer,                intent(in)    :: gamma_kind
        integer,                intent(in)    :: maxit
        real(wp),               intent(in)    :: tol
        logical,                intent(in)    :: write_matrix
        type(property_mixer_t), intent(in)    :: mixing
        type(dftbstate_t),      intent(inout) :: st

        integer  :: it, natoms
        real(wp) :: max_diff, e_scc, e_scc_prev, dE_elec, e_coul_it
        real(wp), allocatable :: H0(:,:), S_shift(:,:), V_atom(:), dq_new(:)
        class(mixer_t), allocatable :: mx
        logical :: has_dE

        natoms = struct%natoms
        allocate(H0(st%bas%norb_total, st%bas%norb_total))
        allocate(S_shift(st%bas%norb_total, st%bas%norb_total))
        allocate(V_atom(natoms))
        allocate(dq_new(natoms))

        call build_hs(struct, st%bas, H0, st%S)
        call build_gamma(struct, st%bas, st%gamma, gamma_kind)

        st%dq = 0.0_wp
        call set_occ(st%bas%nelec, st%bas%norb_total, st%occ)

        st%converged = .false.
        st%niter     = 0
        e_scc_prev   = 0.0_wp
        has_dE       = .false.

        call write_dftb_scc_header(.true., maxit, tol)

        do it = 1, maxit
            call coulomb_potential(st%gamma, st%dq, V_atom)
            call build_shift(st%S, V_atom, st%bas, S_shift)
            call build_ham(H0, S_shift, st%H)

            call solve_gen_eig(st%H, st%S, st%eig, st%C)
            call build_density(st%C, st%occ, st%P)
            call mulliken_charges(st%P, st%S, st%bas, st%q)
            call delta_charges(st%q, st%bas, dq_new)

            if (write_matrix) call write_dftb_matrices(H0, st%S, st%P, st%C)

            st%e_band = sum(st%occ * st%eig)
            e_coul_it = coulomb_energy(st%gamma, dq_new)
            e_scc     = st%e_band + e_coul_it

            if (has_dE) then
                dE_elec = abs(e_scc - e_scc_prev)
            else
                dE_elec = 0.0_wp
            end if

            max_diff = maxval(abs(dq_new - st%dq))
            call write_dftb_scc_iter(it, dE_elec, max_diff, has_dE)

            if (max_diff < tol) then
                st%dq        = dq_new
                st%niter     = it
                st%converged = .true.
                exit
            end if

            if (it == 1) call make_mixer(mixing, natoms, mx)
            block
                real(wp) :: dq_mixed(natoms)
                call mx%mix(st%dq, dq_new, dq_mixed)
                st%dq = dq_mixed
            end block

            e_scc_prev = e_scc
            has_dE     = .true.
        end do

        if (.not. st%converged) st%niter = maxit

        call write_dftb_scc_status(st%converged, st%niter)

        if (allocated(mx)) then
            call mx%free()
            deallocate(mx)
        end if
        deallocate(H0, S_shift, V_atom, dq_new)
    end subroutine solve_scc_loop


    !-- Helpers --------------------------------------------------------

    !> S_shift_μν = (1/2) S_μν (V_A(μ) + V_B(ν))
    subroutine build_shift(S, V_atom, bas, S_shift)
        real(wp),             intent(in)  :: S(:,:)
        real(wp),             intent(in)  :: V_atom(:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: S_shift(:,:)

        integer  :: mu, nu, ia, ib, norb
        real(wp) :: V_mu, V_nu

        norb = size(S, 1)
        do mu = 1, norb
            ia   = atom_of_orbital(bas, mu)
            V_mu = V_atom(ia)
            do nu = 1, norb
                ib   = atom_of_orbital(bas, nu)
                V_nu = V_atom(ib)
                S_shift(mu, nu) = 0.5_wp * S(mu, nu) * (V_mu + V_nu)
            end do
        end do
    end subroutine build_shift


    !> H_μν = H0_μν + S_shift_μν
    subroutine build_ham(H0, S_shift, H)
        real(wp), intent(in)  :: H0(:,:), S_shift(:,:)
        real(wp), intent(out) :: H(:,:)
        H = H0 + S_shift
    end subroutine build_ham


    !-- Allocation et utilitaires --------------------------------------

    subroutine ensure_allocated(struct, st, scheme)
        type(structure_t), intent(in)    :: struct
        type(dftbstate_t), intent(inout) :: st
        integer,           intent(in)    :: scheme

        integer :: norb, natoms
        norb   = st%bas%norb_total
        natoms = struct%natoms

        if (.not. allocated(st%H))   allocate(st%H(norb, norb))
        if (.not. allocated(st%S))   allocate(st%S(norb, norb))
        if (.not. allocated(st%C))   allocate(st%C(norb, norb))
        if (.not. allocated(st%eig)) allocate(st%eig(norb))
        if (.not. allocated(st%P))   allocate(st%P(norb, norb))
        if (.not. allocated(st%occ)) allocate(st%occ(norb))
        if (.not. allocated(st%q))   allocate(st%q(natoms))
        if (.not. allocated(st%dq))  allocate(st%dq(natoms))
        if (scheme /= SCHEME_BASIC) then
            if (.not. allocated(st%gamma)) allocate(st%gamma(natoms, natoms))
        end if
    end subroutine ensure_allocated


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
