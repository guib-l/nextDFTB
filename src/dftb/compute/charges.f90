!> Charges et populations Mulliken.
!>
!>   q_A         = Σ_{μ ∈ A} Σ_ν P_{μν} S_{νμ}
!>   dq_A        = q_neutral_A - q_A         (charge partielle, chimie)
!>   q_{A,l}     = Σ_{μ ∈ shell l de A} Σ_ν P_{μν} S_{νμ}
!>   q_{A,μ}    = Σ_ν P_{μν} S_{νμ}
!>
!> Ordre des orbitales par atome (cf. matel) :
!>   1=s, 2=px, 3=py, 4=pz,
!>   5=dxy, 6=dyz, 7=dzx, 8=dx2−y2, 9=d3z2−r2
module charges
    use kinds,     only: wp
    use dftbstate, only: basis_system_t
    implicit none
    private

    public :: atomic_population, partial_atomic_dq
    public :: lshell_population, mshell_population

contains

    !> Population atomique de Mulliken : q_A = Σ_{μ∈A} (P S)_{μμ}.
    subroutine atomic_population(P, S, bas, q)
        real(wp),             intent(in)  :: P(:,:), S(:,:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: q(:)

        integer  :: ia, mu, nu, o0, n
        real(wp) :: ps_diag

        do ia = 1, size(q)
            o0 = bas%atom_orb_start(ia)
            n  = bas%atom_norb(ia)
            ps_diag = 0.0_wp
            do mu = o0, o0 + n - 1
                do nu = 1, size(P, 1)
                    ps_diag = ps_diag + P(mu, nu) * S(nu, mu)
                end do
            end do
            q(ia) = ps_diag
        end do
    end subroutine atomic_population


    !> Charges partielles dq_A = q_neutral_A - q_A (convention chimie).
    subroutine partial_atomic_dq(q, bas, dq)
        real(wp),             intent(in)  :: q(:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: dq(:)
        integer :: ia, ie
        do ia = 1, size(q)
            ie = bas%atom_elem(ia)
            dq(ia) = bas%elems(ie)%q_neutral - q(ia)
        end do
    end subroutine partial_atomic_dq


    !> Population par sous-couche (l-shell). Stockage 1D plat de taille
    !> bas%lshell_orbs ; index de la 1ère sous-couche de A est
    !> bas%atom_lshell_start(A).
    !> Convention : indices locaux 1..nls = (s, p, d) selon l_max.
    subroutine lshell_population(P, S, bas, lshell_q)
        real(wp),             intent(in)  :: P(:,:), S(:,:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: lshell_q(:)

        integer  :: ia, ils, ls0, o0, mu_lo, mu_hi, mu, nu, natoms
        real(wp) :: acc

        lshell_q = 0.0_wp
        natoms = size(bas%atom_norb)
        do ia = 1, natoms
            o0  = bas%atom_orb_start(ia)
            ls0 = bas%atom_lshell_start(ia)
            do ils = 1, bas%atom_nlshell(ia)
                call lshell_orb_range(ils, o0, mu_lo, mu_hi)
                acc = 0.0_wp
                do mu = mu_lo, mu_hi
                    do nu = 1, size(P, 1)
                        acc = acc + P(mu, nu) * S(nu, mu)
                    end do
                end do
                lshell_q(ls0 + ils - 1) = acc
            end do
        end do
    end subroutine lshell_population


    !> Population par orbitale individuelle (m-shell). Stockage 1D plat
    !> de taille bas%norb_total ; un élément par orbitale μ.
    subroutine mshell_population(P, S, bas, mshell_q)
        real(wp),             intent(in)  :: P(:,:), S(:,:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: mshell_q(:)

        integer  :: mu, nu, norb
        real(wp) :: acc

        norb = size(P, 1)
        do mu = 1, norb
            acc = 0.0_wp
            do nu = 1, norb
                acc = acc + P(mu, nu) * S(nu, mu)
            end do
            mshell_q(mu) = acc
        end do
    end subroutine mshell_population


    !> Plage d'orbitales [mu_lo, mu_hi] correspondant à la l-shell
    !> locale `ils` (1=s, 2=p, 3=d) d'un atome dont la 1ère orbitale
    !> est `o0`. Suit l'ordre de matel (1=s ; 2..4=p ; 5..9=d).
    pure subroutine lshell_orb_range(ils, o0, mu_lo, mu_hi)
        integer, intent(in)  :: ils, o0
        integer, intent(out) :: mu_lo, mu_hi
        select case (ils)
        case (1)
            mu_lo = o0;     mu_hi = o0
        case (2)
            mu_lo = o0 + 1; mu_hi = o0 + 3
        case (3)
            mu_lo = o0 + 4; mu_hi = o0 + 8
        case default
            mu_lo = o0;     mu_hi = o0 - 1
        end select
    end subroutine lshell_orb_range
end module charges
