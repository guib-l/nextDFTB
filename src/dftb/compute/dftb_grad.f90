!> Gradient DFTB (forces sur chaque atome).
!>
!> Décomposition standard SCC-DFTB :
!>   F_{kx} = F^{rep}_{kx} + F^{H0}_{kx} + F^{norm}_{kx} + F^{scc}_{kx}
!>
!> - F^{rep}    : dérivée du potentiel répulsif (paires SK).
!> - F^{H0}     : dérivée du terme bandes via H0 ;
!>                  -2 Σ_{a≠k} Σ_{μ∈a, ν∈k} P_{μν} dH0_{μν}/dr_{kx}
!> - F^{norm}   : dérivée de la contrainte de normalisation S ;
!>                  +2 Σ_{a≠k} Σ_{μ∈a, ν∈k} W_{μν} dS_{μν}/dr_{kx}
!>                où W = Σ_i occ_i eig_i C_μi C_νi
!> - F^{scc}    : - Σ_{μν} 2 P_{μν} (1/2)(V_a + V_b) dS_{μν}/dr_{kx}
!>                - dq_k Σ_{a≠k} dq_a dγ_{ka}/dr_{kx}
module dftb_grad
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use dftbstate,     only: dftbstate_t
    use property,      only: property_dftb_t
    use matel,         only: build_off_dblock
    use gamma_mod,     only: build_dgamma
    use repulsif,     only: repulsive_grad
    implicit none
    private

    logical, parameter, public :: SPLIT_FORCE = .true.

    public :: compute_gradient
    public :: grad_erep, grad_ham, grad_scc, grad_ecoulomb

contains

    !> Calcule le gradient total et le stocke dans `state%grad(3,natoms)`.
    subroutine compute_gradient(struct, state, prop)
        type(structure_t),     intent(in)    :: struct
        type(dftbstate_t),     intent(inout) :: state
        type(property_dftb_t), intent(in)    :: prop

        real(wp), allocatable :: g_rep(:,:), g_ham(:,:)
        real(wp), allocatable :: g_norm(:,:), g_scc_shift(:,:), g_coul(:,:)
        integer :: natoms

        natoms = struct%natoms

        if (SPLIT_FORCE) then
            allocate(g_rep(3, natoms), g_ham(3, natoms))
            allocate(g_norm(3, natoms), g_scc_shift(3, natoms))
            allocate(g_coul(3, natoms))

            g_rep       = grad_erep(struct)
            g_ham       = grad_ham(struct, state)
            g_norm      = grad_norm(struct, state)
            g_scc_shift = grad_scc(struct, state)
            g_coul      = grad_ecoulomb(struct, state, prop%gamma_kind)

            state%grad = g_rep + g_ham + g_norm + g_scc_shift + g_coul

            deallocate(g_rep, g_ham, g_norm, g_scc_shift, g_coul)
        else
            ! Mode regroupé (non distinct ici) : même résultat numérique.
            allocate(g_rep(3, natoms), g_ham(3, natoms))
            allocate(g_norm(3, natoms), g_scc_shift(3, natoms))
            allocate(g_coul(3, natoms))

            g_rep       = grad_erep(struct)
            g_ham       = grad_ham(struct, state)
            g_norm      = grad_norm(struct, state)
            g_scc_shift = grad_scc(struct, state)
            g_coul      = grad_ecoulomb(struct, state, prop%gamma_kind)

            state%grad = g_rep + g_ham + g_norm + g_scc_shift + g_coul

            deallocate(g_rep, g_ham, g_norm, g_scc_shift, g_coul)
        end if
    end subroutine compute_gradient


    !> Gradient répulsif : g(:, k) = -dE_rep/dr_k.
    function grad_erep(struct) result(g)
        type(structure_t), intent(in) :: struct
        real(wp), allocatable :: g(:,:)
        integer :: k, d

        allocate(g(3, struct%natoms))
        g = 0.0_wp
        do k = 1, struct%natoms
            do d = 1, 3
                g(d, k) = -repulsive_grad(struct, k, d)
            end do
        end do
    end function grad_erep


    !> Gradient H0 :  g(d, k) = -2 Σ_{a≠k} Σ_{μ∈a, ν∈k} P_{μν} dH_{μν}/dr_{kx}
    !> Le facteur 2 provient de la contribution symétrique a↔k.
    function grad_ham(struct, state) result(g)
        type(structure_t), intent(in) :: struct
        type(dftbstate_t), intent(in) :: state
        real(wp), allocatable :: g(:,:)
        integer :: k, j, d, oi, oj, ni, nj
        real(wp), allocatable :: dhblk(:,:), dsblk(:,:)
        real(wp) :: acc

        allocate(g(3, struct%natoms))
        g = 0.0_wp

        do k = 1, struct%natoms
            ni = state%bas%atom_norb(k)
            oi = state%bas%atom_orb_start(k)
            do j = 1, struct%natoms
                if (j == k) cycle
                nj = state%bas%atom_norb(j)
                oj = state%bas%atom_orb_start(j)
                allocate(dhblk(ni, nj), dsblk(ni, nj))
                do d = 1, 3
                    call build_off_dblock(struct, state%bas, k, j, k, d, &
                                          dhblk, dsblk)
                    acc = block_trace(state%P, oi, ni, oj, nj, dhblk)
                    g(d, k) = g(d, k) - 2.0_wp * acc
                end do
                deallocate(dhblk, dsblk)
            end do
        end do
    end function grad_ham


    !> Gradient du terme de normalisation (énergie pondérée par eig) :
    !>   g(d, k) = +2 Σ_{a≠k} Σ_{μ∈a, ν∈k} W_{μν} dS_{μν}/dr_{kx}
    !> Sign opposé de F^{H0} car cette pièce vient de la contrainte
    !> d'orthonormalité (terme -2 eig_i dS dans la formule générale).
    function grad_norm(struct, state) result(g)
        type(structure_t), intent(in) :: struct
        type(dftbstate_t), intent(in) :: state
        real(wp), allocatable :: g(:,:)
        integer :: k, j, d, oi, oj, ni, nj
        real(wp), allocatable :: dhblk(:,:), dsblk(:,:)
        real(wp) :: acc

        allocate(g(3, struct%natoms))
        g = 0.0_wp

        do k = 1, struct%natoms
            ni = state%bas%atom_norb(k)
            oi = state%bas%atom_orb_start(k)
            do j = 1, struct%natoms
                if (j == k) cycle
                nj = state%bas%atom_norb(j)
                oj = state%bas%atom_orb_start(j)
                allocate(dhblk(ni, nj), dsblk(ni, nj))
                do d = 1, 3
                    call build_off_dblock(struct, state%bas, k, j, k, d, &
                                          dhblk, dsblk)
                    acc = block_trace(state%W, oi, ni, oj, nj, dsblk)
                    g(d, k) = g(d, k) + 2.0_wp * acc
                end do
                deallocate(dhblk, dsblk)
            end do
        end do
    end function grad_norm


    !> Partie h_shift × dS du F^{scc} :
    !>   g(d, k) = - Σ_{a≠k} Σ_{μ∈a, ν∈k} 2 P_{μν} (1/2)(V_a + V_k)
    !>                                       × dS_{μν}/dr_{kx}
    !> Le facteur 1/2 vient de la convention shift (S_shift =
    !> (1/2) S (V_a+V_b)). Le facteur 2 vient de la symétrie a↔k.
    function grad_scc(struct, state) result(g)
        type(structure_t), intent(in) :: struct
        type(dftbstate_t), intent(in) :: state
        real(wp), allocatable :: g(:,:)
        real(wp), allocatable :: V_atom(:)
        integer  :: k, j, d, oi, oj, ni, nj, natoms, kk
        real(wp), allocatable :: dhblk(:,:), dsblk(:,:)
        real(wp) :: acc, v_sum

        natoms = struct%natoms
        allocate(g(3, natoms))
        g = 0.0_wp

        if (.not. allocated(state%gamma)) return

        ! V_a = -Σ_K γ_aK dq_K (convention chimie).
        allocate(V_atom(natoms))
        do k = 1, natoms
            v_sum = 0.0_wp
            do kk = 1, natoms
                v_sum = v_sum - state%gamma(k, kk) * state%dq(kk)
            end do
            V_atom(k) = v_sum
        end do

        do k = 1, natoms
            ni = state%bas%atom_norb(k)
            oi = state%bas%atom_orb_start(k)
            do j = 1, natoms
                if (j == k) cycle
                nj = state%bas%atom_norb(j)
                oj = state%bas%atom_orb_start(j)
                allocate(dhblk(ni, nj), dsblk(ni, nj))
                do d = 1, 3
                    call build_off_dblock(struct, state%bas, k, j, k, d, &
                                          dhblk, dsblk)
                    acc = block_trace(state%P, oi, ni, oj, nj, dsblk)
                    g(d, k) = g(d, k) - (V_atom(k) + V_atom(j)) * acc
                end do
                deallocate(dhblk, dsblk)
            end do
        end do

        deallocate(V_atom)
    end function grad_scc


    !> Gradient énergie SCC coulombienne :
    !>   g(d, k) = - dq_k Σ_{a≠k} dq_a dγ_{ka}/dr_{kx}
    !> Note : le facteur global 1/2 de E_coul est compensé par les deux
    !> termes symétriques (a,k) et (k,a) qui contribuent à la dérivée
    !> par rapport à r_k.
    function grad_ecoulomb(struct, state, kind) result(g)
        type(structure_t), intent(in) :: struct
        type(dftbstate_t), intent(in) :: state
        integer,           intent(in) :: kind
        real(wp), allocatable :: g(:,:)
        real(wp), allocatable :: dgamma(:,:)
        integer  :: k, j, d, natoms
        real(wp) :: acc

        natoms = struct%natoms
        allocate(g(3, natoms))
        g = 0.0_wp

        if (.not. allocated(state%gamma) .or. .not. allocated(state%dq)) return

        allocate(dgamma(natoms, natoms))
        do k = 1, natoms
            do d = 1, 3
                call build_dgamma(struct, state%bas, k, d, dgamma, kind)
                acc = 0.0_wp
                do j = 1, natoms
                    if (j == k) cycle
                    acc = acc + state%dq(j) * dgamma(k, j)
                end do
                g(d, k) = g(d, k) - state%dq(k) * acc
            end do
        end do
        deallocate(dgamma)
    end function grad_ecoulomb


    !-- Helpers internes -----------------------------------------------

    !> Trace partielle d'un bloc (ni × nj) : Σ_{μ,ν} M_{oi+μ-1, oj+ν-1} B_{μν}.
    pure function block_trace(M, oi, ni, oj, nj, B) result(s)
        real(wp), intent(in) :: M(:,:)
        integer,  intent(in) :: oi, ni, oj, nj
        real(wp), intent(in) :: B(:,:)
        real(wp) :: s
        integer  :: mu, nu

        s = 0.0_wp
        do nu = 1, nj
            do mu = 1, ni
                s = s + M(oi + mu - 1, oj + nu - 1) * B(mu, nu)
            end do
        end do
    end function block_trace

end module dftb_grad
