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
    use kinds,        only: wp
    use dftbstate,    only: basis_system_t
    use orbitals_mod, only: atom_of_orbital, lshell_of_orbital, lshell_norb
    implicit none
    private

    integer, parameter, public :: SHIFT_DQ     = 1
    integer, parameter, public :: SHIFT_LSHELL = 2

    public :: build_shift, build_shift_dq, build_shift_lshell

contains

    !> Sélecteur. Par défaut : SHIFT_LSHELL (fallback sur SHIFT_DQ
    !> si `lshell_dq` n'est pas fourni).
    subroutine build_shift(S, gamma, dq, bas, S_shift, kind, lshell_dq)
        real(wp),             intent(in)  :: S(:,:)
        real(wp),             intent(in)  :: gamma(:,:)
        real(wp),             intent(in)  :: dq(:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: S_shift(:,:)
        integer,  optional,   intent(in)  :: kind
        real(wp), optional,   intent(in)  :: lshell_dq(:)

        integer :: k

        k = SHIFT_LSHELL
        if (present(kind)) k = kind

        select case (k)
        case (SHIFT_LSHELL)
            if (.not. present(lshell_dq)) then
                call build_shift_dq(S, gamma, dq, bas, S_shift)
            else
                call build_shift_lshell(S, gamma, dq, lshell_dq, bas, S_shift)
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
    !> par sous-couche du site local. La contribution locale (K = A) est
    !> dq_{A,l_μ}/n_orb(A,l_μ), tandis que les sites K ≠ A contribuent
    !> avec leur dq_K total. Idem pour ν. Signe global -1 car
    !> V = -Σ γ dq dans cette convention. `dq_atom` et `lshell_dq` sont
    !> attendus précalculés (convention chimie).
    subroutine build_shift_lshell(S, gamma, dq_atom, lshell_dq, bas, S_shift)
        real(wp),             intent(in)  :: S(:,:)
        real(wp),             intent(in)  :: gamma(:,:)
        real(wp),             intent(in)  :: dq_atom(:)
        real(wp),             intent(in)  :: lshell_dq(:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: S_shift(:,:)

        integer  :: mu, nu, ia, ib, kk, natoms, norb, l_mu, l_nu, n_loc
        real(wp) :: V_mu, V_nu, dq_local

        natoms = size(bas%atom_norb)
        norb   = size(S, 1)

        do mu = 1, norb
            ia   = atom_of_orbital(bas, mu)
            l_mu = lshell_of_orbital(bas, mu)
            do nu = 1, norb
                ib   = atom_of_orbital(bas, nu)
                l_nu = lshell_of_orbital(bas, nu)

                V_mu = 0.0_wp
                do kk = 1, natoms
                    if (kk == ia) then
                        n_loc = lshell_norb(bas, ia, l_mu)
                        dq_local = lshell_dq(bas%atom_lshell_start(ia) + l_mu - 1)
                        V_mu = V_mu - gamma(ia, kk) * (dq_local / real(n_loc, wp))
                    else
                        V_mu = V_mu - gamma(ia, kk) * dq_atom(kk)
                    end if
                end do

                V_nu = 0.0_wp
                do kk = 1, natoms
                    if (kk == ib) then
                        n_loc = lshell_norb(bas, ib, l_nu)
                        dq_local = lshell_dq(bas%atom_lshell_start(ib) + l_nu - 1)
                        V_nu = V_nu - gamma(ib, kk) * (dq_local / real(n_loc, wp))
                    else
                        V_nu = V_nu - gamma(ib, kk) * dq_atom(kk)
                    end if
                end do

                S_shift(mu, nu) = 0.5_wp * S(mu, nu) * (V_mu + V_nu)
            end do
        end do
    end subroutine build_shift_lshell

end module shift
