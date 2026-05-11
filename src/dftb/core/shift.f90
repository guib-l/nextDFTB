!> Construction de la matrice de shift H1 (terme SCC du DFTB) :
!>
!>   H1_{μν} = (1/2) S_{μν} (V_μ + V_ν)        avec μ ∈ A, ν ∈ B
!>
!> Le potentiel V_μ vu par l'orbitale μ est défini selon le choix :
!>   - SHIFT_DQ      : V_A = Σ_K γ_AK (q_K - q_neutral_K) = -Σ_K γ_AK dq_K
!>                     (convention chimie : dq_K = q_neutral_K - q_K).
!>   - SHIFT_LSHELL  : V_μ = -Σ_K γ_AK dq_K  pour les contributions
!>                     non-locales K ≠ A, et la charge locale dq_A est
!>                     redistribuée par l-shell sur les n_orb(A,l_μ)
!>                     orbitales du shell de μ.
module shift
    use kinds,     only: wp
    use dftbstate, only: basis_system_t
    implicit none
    private

    integer, parameter, public :: SHIFT_DQ     = 1
    integer, parameter, public :: SHIFT_LSHELL = 2

    public :: build_shift, build_shift_dq, build_shift_lshell

contains

    !> Sélecteur. Par défaut : SHIFT_LSHELL (fallback sur SHIFT_DQ
    !> si `lshell_q` n'est pas fourni).
    subroutine build_shift(S, gamma, dq, bas, S_shift, kind, lshell_q)
        real(wp),             intent(in)  :: S(:,:)
        real(wp),             intent(in)  :: gamma(:,:)
        real(wp),             intent(in)  :: dq(:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: S_shift(:,:)
        integer,  optional,   intent(in)  :: kind
        real(wp), optional,   intent(in)  :: lshell_q(:)

        integer :: k

        k = SHIFT_LSHELL
        if (present(kind)) k = kind

        select case (k)
        case (SHIFT_LSHELL)
            if (.not. present(lshell_q)) then
                call build_shift_dq(S, gamma, dq, bas, S_shift)
            else
                call build_shift_lshell(S, gamma, lshell_q, bas, S_shift)
            end if
        case default
            call build_shift_dq(S, gamma, dq, bas, S_shift)
        end select
    end subroutine build_shift


    !> H1_{μν} = (1/2) S_{μν} (V_A + V_B), V_A = -Σ_K γ_AK dq_K
    !> (convention chimie : dq_K = q_neutral_K - q_K).
    subroutine build_shift_dq(S, gamma, dq, bas, S_shift)
        real(wp),             intent(in)  :: S(:,:)
        real(wp),             intent(in)  :: gamma(:,:)
        real(wp),             intent(in)  :: dq(:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: S_shift(:,:)

        integer  :: mu, nu, ia, ib, kk, natoms, norb
        real(wp), allocatable :: V_atom(:)
        real(wp) :: acc

        natoms = size(dq)
        norb   = size(S, 1)
        allocate(V_atom(natoms))
        do ia = 1, natoms
            acc = 0.0_wp
            do kk = 1, natoms
                acc = acc - gamma(ia, kk) * dq(kk)
            end do
            V_atom(ia) = acc
        end do

        do mu = 1, norb
            ia = atom_of_orbital(bas, mu)
            do nu = 1, norb
                ib = atom_of_orbital(bas, nu)
                S_shift(mu, nu) = 0.5_wp * S(mu, nu) * (V_atom(ia) + V_atom(ib))
            end do
        end do

        deallocate(V_atom)
    end subroutine build_shift_dq


    !> Variante l-shell : la charge atomique vue par μ est redistribuée
    !> par sous-couche du site local. On note dq_{A,l} = q_neutral_{A,l}
    !> - lshell_q_{A,l} (convention chimie). La contribution locale
    !> (K = A) est dq_{A,l_μ}/n_orb(A,l_μ), tandis que les sites K ≠ A
    !> contribuent avec leur dq_K total. Idem pour ν. Signe global -1
    !> car V = -Σ γ dq dans cette convention.
    subroutine build_shift_lshell(S, gamma, lshell_q, bas, S_shift)
        real(wp),             intent(in)  :: S(:,:)
        real(wp),             intent(in)  :: gamma(:,:)
        real(wp),             intent(in)  :: lshell_q(:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: S_shift(:,:)

        integer  :: mu, nu, ia, ib, kk, natoms, norb, l_mu, l_nu, n_loc
        real(wp), allocatable :: dq_atom(:), dq_lshell(:)
        real(wp) :: V_mu, V_nu, dq_local

        natoms = size(bas%atom_norb)
        norb   = size(S, 1)
        allocate(dq_atom(natoms))
        allocate(dq_lshell(size(lshell_q)))

        call lshell_dq(bas, lshell_q, dq_lshell)
        call sum_dq_per_atom(bas, dq_lshell, dq_atom)

        do mu = 1, norb
            ia   = atom_of_orbital(bas, mu)
            l_mu = lshell_of_orbital(bas, mu)
            do nu = 1, norb
                ib   = atom_of_orbital(bas, nu)
                l_nu = lshell_of_orbital(bas, nu)

                ! V_μ : somme sur tous les sites K, charge atomique sauf
                ! pour le site local A où l'on prend dq_{A,l_μ}/n_orb.
                V_mu = 0.0_wp
                do kk = 1, natoms
                    if (kk == ia) then
                        n_loc = lshell_norb(bas, ia, l_mu)
                        dq_local = dq_lshell(bas%atom_lshell_start(ia) + l_mu - 1)
                        V_mu = V_mu - gamma(ia, kk) * (dq_local / real(n_loc, wp))
                    else
                        V_mu = V_mu - gamma(ia, kk) * dq_atom(kk)
                    end if
                end do

                V_nu = 0.0_wp
                do kk = 1, natoms
                    if (kk == ib) then
                        n_loc = lshell_norb(bas, ib, l_nu)
                        dq_local = dq_lshell(bas%atom_lshell_start(ib) + l_nu - 1)
                        V_nu = V_nu - gamma(ib, kk) * (dq_local / real(n_loc, wp))
                    else
                        V_nu = V_nu - gamma(ib, kk) * dq_atom(kk)
                    end if
                end do

                S_shift(mu, nu) = 0.5_wp * S(mu, nu) * (V_mu + V_nu)
            end do
        end do

        deallocate(dq_atom, dq_lshell)
    end subroutine build_shift_lshell


    !-- Utilitaires ----------------------------------------------------

    !> dq par sous-couche : dq_{A,l} = q_neutral_{A,l} - lshell_q_{A,l}
    !> (convention chimie : charge partielle).
    subroutine lshell_dq(bas, lshell_q, dq_lshell)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(in)  :: lshell_q(:)
        real(wp),             intent(out) :: dq_lshell(:)

        integer  :: ia, ils, ie, ls0, natoms
        real(wp) :: q_neutral_l

        natoms = size(bas%atom_norb)
        do ia = 1, natoms
            ie  = bas%atom_elem(ia)
            ls0 = bas%atom_lshell_start(ia)
            do ils = 1, bas%atom_nlshell(ia)
                select case (ils)
                case (1); q_neutral_l = bas%elems(ie)%occ_s
                case (2); q_neutral_l = bas%elems(ie)%occ_p
                case (3); q_neutral_l = bas%elems(ie)%occ_d
                case default; q_neutral_l = 0.0_wp
                end select
                dq_lshell(ls0 + ils - 1) = q_neutral_l - lshell_q(ls0 + ils - 1)
            end do
        end do
    end subroutine lshell_dq


    !> Somme des dq_lshell sur les sous-couches d'un atome → dq_atom.
    subroutine sum_dq_per_atom(bas, dq_lshell, dq_atom)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(in)  :: dq_lshell(:)
        real(wp),             intent(out) :: dq_atom(:)

        integer  :: ia, ils, ls0, natoms
        real(wp) :: acc

        natoms = size(bas%atom_norb)
        do ia = 1, natoms
            ls0 = bas%atom_lshell_start(ia)
            acc = 0.0_wp
            do ils = 1, bas%atom_nlshell(ia)
                acc = acc + dq_lshell(ls0 + ils - 1)
            end do
            dq_atom(ia) = acc
        end do
    end subroutine sum_dq_per_atom


    !> Atome auquel appartient l'orbitale globale `iorb`.
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


    !> Indice local de l-shell (1=s, 2=p, 3=d) auquel appartient
    !> l'orbitale globale `iorb`. Suit l'ordre de matel.
    pure function lshell_of_orbital(bas, iorb) result(ils)
        type(basis_system_t), intent(in) :: bas
        integer,              intent(in) :: iorb
        integer :: ils, ia, local
        ia    = atom_of_orbital(bas, iorb)
        local = iorb - bas%atom_orb_start(ia) + 1
        if (local == 1) then
            ils = 1
        else if (local <= 4) then
            ils = 2
        else
            ils = 3
        end if
    end function lshell_of_orbital


    !> Nombre d'orbitales du shell `ils` (1=s, 2=p, 3=d) sur l'atome ia.
    pure function lshell_norb(bas, ia, ils) result(n)
        type(basis_system_t), intent(in) :: bas
        integer,              intent(in) :: ia, ils
        integer :: n
        select case (ils)
        case (1); n = 1
        case (2); n = 3
        case (3); n = 5
        case default; n = 1
        end select
        if (n > bas%atom_norb(ia)) n = bas%atom_norb(ia)
    end function lshell_norb
end module shift
