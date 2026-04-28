!> Résolution SCC-DFTB.
!>
!> Boucle :
!>   1. construire H0, S, gamma (une fois)
!>   2. Δq = 0 initial
!>   3. répéter :
!>        H_μν = H0_μν + (1/2) S_μν (V_A + V_B)   où V_A = Σ_C γ_AC Δq_C
!>        diagonaliser H C = S C E
!>        remplir occ (Aufbau, closed-shell)
!>        P, q (Mulliken), Δq_new
!>        si max|Δq_new - Δq| < tol → fin
!>        mixer Δq
module scc
    use kinds,        only: wp
    use globals,      only: geometry_t, basis_system_t
    use readskf,      only: skf_store_t
    use dftb_matrix,  only: build_hs
    use gamma_mod,    only: build_gamma
    use density,      only: build_density
    use diag,         only: gen_eig
    use charges,      only: mulliken_charges, delta_charges
    use coulomb,      only: coulomb_potential
    use mixer,        only: mixer_t, init_mixer, free_mixer, mix_charges
    use logger,       only: log_msg
    implicit none
    private

    type, public :: scc_result_t
        real(wp), allocatable :: H(:,:)          ! Hamiltonien final
        real(wp), allocatable :: S(:,:)
        real(wp), allocatable :: C(:,:)          ! coefficients MO
        real(wp), allocatable :: eig(:)
        real(wp), allocatable :: P(:,:)
        real(wp), allocatable :: occ(:)
        real(wp), allocatable :: q(:)            ! charges Mulliken (atomes)
        real(wp), allocatable :: dq(:)
        real(wp), allocatable :: gamma(:,:)
        real(wp) :: e_band   = 0.0_wp
        integer  :: niter    = 0
        logical  :: converged = .false.
    end type scc_result_t

    public :: solve_scc

contains

    subroutine solve_scc(geom, bas, store, do_scc, maxit, tol, res)
        type(geometry_t),     intent(in)  :: geom
        type(basis_system_t), intent(in)  :: bas
        type(skf_store_t),    intent(in)  :: store
        logical,              intent(in)  :: do_scc
        integer,              intent(in)  :: maxit
        real(wp),             intent(in)  :: tol
        type(scc_result_t),   intent(out) :: res

        integer  :: norb, natoms, it, mu, nu, ia, ib
        real(wp) :: max_diff, V_mu, V_nu
        real(wp), allocatable :: H0(:,:), V_atom(:), dq_new(:)
        type(mixer_t) :: mx
        character(len=128) :: buf

        norb   = bas%norb_total
        natoms = geom%natoms

        allocate(res%H(norb, norb), res%S(norb, norb), res%C(norb, norb))
        allocate(res%eig(norb), res%P(norb, norb), res%occ(norb))
        allocate(res%q(natoms), res%dq(natoms), res%gamma(natoms, natoms))
        allocate(H0(norb, norb), V_atom(natoms), dq_new(natoms))

        call build_hs(geom, bas, store, H0, res%S)
        call build_gamma(geom, bas, res%gamma)

        res%dq = 0.0_wp

        ! occupations (closed-shell, Aufbau) — calculées une fois (nelec fixe)
        call set_occ(bas%nelec, norb, res%occ)

        do it = 1, maxit
            ! 1. potentiel atomique V_A = Σ_B γ_AB Δq_B
            call coulomb_potential(res%gamma, res%dq, V_atom)

            ! 2. H_μν = H0_μν + 0.5 S_μν (V_A + V_B)
            res%H = H0
            if (do_scc) then
                do mu = 1, norb
                    ia = atom_of_orbital(bas, mu)
                    V_mu = V_atom(ia)
                    do nu = 1, norb
                        ib = atom_of_orbital(bas, nu)
                        V_nu = V_atom(ib)
                        res%H(mu, nu) = res%H(mu, nu) &
                            + 0.5_wp * res%S(mu, nu) * (V_mu + V_nu)
                    end do
                end do
            end if

            ! 3. diagonalisation
            call gen_eig(res%H, res%S, res%eig, res%C)

            ! 4. P, q, Δq
            call build_density(res%C, res%occ, res%P)
            call mulliken_charges(res%P, res%S, bas, res%q)
            call delta_charges(res%q, bas, dq_new)

            res%e_band = sum(res%occ * res%eig)

            if (.not. do_scc) then
                res%dq = dq_new
                res%niter = 1
                res%converged = .true.
                exit
            end if

            max_diff = maxval(abs(dq_new - res%dq))
            write(buf, '(a,i4,a,es12.4)') "  scc iter ", it, "  max|Δq| diff = ", max_diff
            call log_msg(buf)

            if (max_diff < tol) then
                res%dq = dq_new
                res%niter = it
                res%converged = .true.
                exit
            end if

            if (it == 1) call init_mixer(mx, natoms, 0)
            block
                real(wp) :: dq_mixed(natoms)
                call mix_charges(mx, res%dq, dq_new, dq_mixed)
                res%dq = dq_mixed
            end block
        end do

        if (.not. res%converged) then
            res%niter = maxit
            res%converged = .false.
        end if

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
                iorb < bas%atom_orb_start(k) + bas%atom_norb(k)) then
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
