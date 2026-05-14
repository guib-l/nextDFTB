!> Résolution DFTB selon trois schémas distincts :
!>
!>   - BASIC : pas de terme dépendant des charges, une seule
!>             diagonalisation H = H0.
!>   - NOSCC : une diagonalisation avec V calculé à partir des charges
!>             fournies dans l'input de géométrie (pas d'auto-cohérence).
!>   - SCC   : cycle SCC complet jusqu'à convergence sur Δq.
!>
!> Le shift H1_{μν} = (1/2) S_{μν} (V_A + V_B) est délégué au module
!> `shift`. Le helper `build_ham` ajoute simplement H = H0 + H1.
module scc
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use dftbstate,     only: dftbstate_t, basis_system_t
    use matel,         only: build_hs
    use gamma_mod,     only: build_gamma
    use density,       only: build_density
    use linalg,        only: solve_gen_eig
    use charges,       only: build_charges
    use coulomb,       only: coulomb_energy
    use shift,         only: build_shift, SHIFT_DQ
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
        st%H0 = st%H
        call set_occ(st%bas%nelec, st%bas%norb_total, st%occ)

        call solve_gen_eig(st%H, st%S, st%eig, st%C)
        call build_density(st%C, st%occ, st%P)
        call build_charges(st%P, st%S, st)

        if (write_matrix) call write_dftb_matrices(st%H, st%S, st%P, st%C, st%gamma)

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
        real(wp), allocatable :: S_shift(:,:)

        natoms = struct%natoms
        allocate(S_shift(st%bas%norb_total, st%bas%norb_total))

        call build_hs(struct, st%bas, st%H0, st%S)
        call build_gamma(struct, st%bas, st%gamma, gamma_kind)

        ! Δq initial = charge formelle fournie par l'input géométrie.
        ! Convention : charge>0 (cation) ⇒ moins d'électrons ⇒ Δq>0.
        do ia = 1, natoms
            st%dq(ia) = struct%atoms(ia)%charge
        end do

        ! Variante atomique imposée : aucun lshell_q n'est calculé
        ! avant la première (et unique) diagonalisation.
        call build_shift(st%S, st%gamma, st%dq, st%bas, S_shift, kind=SHIFT_DQ)
        call build_ham(st%H0, S_shift, st%H)

        call set_occ(st%bas%nelec, st%bas%norb_total, st%occ)
        call solve_gen_eig(st%H, st%S, st%eig, st%C)
        call build_density(st%C, st%occ, st%P)
        call build_charges(st%P, st%S, st)

        if (write_matrix) call write_dftb_matrices(st%H0, st%S, st%P, st%C, st%gamma)

        st%e_band    = sum(st%occ * st%eig)
        st%niter     = 1
        st%converged = .true.

        deallocate(S_shift)
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

        integer  :: it, ia, natoms
        real(wp) :: max_diff, e_scc, e_scc_prev, dE_elec, e_coul_it
        real(wp), allocatable :: S_shift(:,:), dq_old(:)
        class(mixer_t), allocatable :: mx
        logical :: has_dE

        natoms = struct%natoms
        allocate(S_shift(st%bas%norb_total, st%bas%norb_total))
        allocate(dq_old(natoms))

        call build_hs(struct, st%bas, st%H0, st%S)
        call build_gamma(struct, st%bas, st%gamma, gamma_kind)

        ! Δq initial = charge formelle fournie par l'input géométrie
        ! (0 par défaut). Convention : charge>0 ⇒ Δq>0.
        do ia = 1, natoms
            st%dq(ia) = struct%atoms(ia)%charge
        end do
        st%lshell_q  = 0.0_wp
        st%lshell_dq = 0.0_wp
        call set_occ(st%bas%nelec, st%bas%norb_total, st%occ)

        st%converged = .false.
        st%niter     = 0
        e_scc_prev   = 0.0_wp
        has_dE       = .false.

        call write_dftb_scc_header(.true., maxit, tol)

        do it = 1, maxit
            ! Shift atomique : V_A = -Σ_K γ_AK dq_K. L'auto-cohérence
            ! porte exclusivement sur les charges atomiques dq.
            call build_shift(st%S, st%gamma, st%dq, st%bas, S_shift, &
                              kind=SHIFT_DQ)
            call build_ham(st%H0, S_shift, st%H)

            call solve_gen_eig(st%H, st%S, st%eig, st%C)
            call build_density(st%C, st%occ, st%P)

            ! Mémoriser dq d'entrée avant écrasement par les nouvelles
            ! charges Mulliken : sert au test de convergence et au mixage.
            dq_old = st%dq
            call build_charges(st%P, st%S, st)

            if (write_matrix) call write_dftb_matrices(st%H0, st%S, st%P, st%C, st%gamma)

            st%e_band = sum(st%occ * st%eig)
            e_coul_it = coulomb_energy(st%gamma, st%dq)
            e_scc     = st%e_band + e_coul_it

            if (has_dE) then
                dE_elec = abs(e_scc - e_scc_prev)
            else
                dE_elec = 0.0_wp
            end if

            max_diff = maxval(abs(st%dq - dq_old))
            call write_dftb_scc_iter(it, dE_elec, max_diff, has_dE)

            if (max_diff < tol) then
                st%niter     = it
                st%converged = .true.
                exit
            end if

            ! Mixage des charges atomiques : c'est la seule quantité
            ! pilotant build_shift_dq, donc la grandeur naturelle de la
            ! boucle SCC.
            if (it == 1) call make_mixer(mixing, natoms, mx)
            block
                real(wp) :: dq_mixed(natoms)
                call mx%mix(dq_old, st%dq, dq_mixed)
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
        deallocate(S_shift, dq_old)
    end subroutine solve_scc_loop


    !-- Helpers --------------------------------------------------------

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

        integer :: norb, natoms, nls
        norb   = st%bas%norb_total
        natoms = struct%natoms
        nls    = st%bas%lshell_orbs

        if (.not. allocated(st%H))         allocate(st%H(norb, norb))
        if (.not. allocated(st%H0))        allocate(st%H0(norb, norb))
        if (.not. allocated(st%S))         allocate(st%S(norb, norb))
        if (.not. allocated(st%C))         allocate(st%C(norb, norb))
        if (.not. allocated(st%eig))       allocate(st%eig(norb))
        if (.not. allocated(st%P))         allocate(st%P(norb, norb))
        if (.not. allocated(st%occ))       allocate(st%occ(norb))
        if (.not. allocated(st%q))         allocate(st%q(natoms))
        if (.not. allocated(st%dq))        allocate(st%dq(natoms))
        if (.not. allocated(st%lshell_q))  allocate(st%lshell_q(nls))
        if (.not. allocated(st%lshell_dq)) allocate(st%lshell_dq(nls))
        if (.not. allocated(st%mshell_q))  allocate(st%mshell_q(norb))
        if (scheme /= SCHEME_BASIC) then
            if (.not. allocated(st%gamma)) allocate(st%gamma(natoms, natoms))
        end if
    end subroutine ensure_allocated


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
