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
    use dftbstate, only: basis_system_t, dftbstate_t
    implicit none
    private

    public :: partial_atomic_dq, partial_lshell_dq
    public :: mulliken_population, build_charges

contains

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


    !> Charges partielles par sous-couche :
    !>   dq_{A,l} = q_neutral_{A,l} - lshell_q_{A,l}.
    subroutine partial_lshell_dq(bas, lshell_q, lshell_dq)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(in)  :: lshell_q(:)
        real(wp),             intent(out) :: lshell_dq(:)

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
                lshell_dq(ls0 + ils - 1) = q_neutral_l - lshell_q(ls0 + ils - 1)
            end do
        end do
    end subroutine partial_lshell_dq


    !> Populations Mulliken en une seule passe :
    !>   mshell_q_μ = Σ_ν P_{μν} S_{νμ}
    !>   lshell_q_{A,l} = Σ_{μ ∈ shell l de A} mshell_q_μ
    !>   q_A           = Σ_l lshell_q_{A,l}
    subroutine mulliken_population(P, S, bas, q, lshell_q, mshell_q)
        real(wp),             intent(in)  :: P(:,:), S(:,:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: q(:)
        real(wp),             intent(out) :: lshell_q(:)
        real(wp),             intent(out) :: mshell_q(:)

        integer  :: ia, ils, mu, nu, o0, n, ls0, local, norb, natoms
        real(wp) :: acc

        norb   = size(P, 1)
        natoms = size(bas%atom_norb)

        do mu = 1, norb
            acc = 0.0_wp
            do nu = 1, norb
                acc = acc + P(mu, nu) * S(nu, mu)
            end do
            mshell_q(mu) = acc
        end do

        lshell_q = 0.0_wp
        do ia = 1, natoms
            o0  = bas%atom_orb_start(ia)
            n   = bas%atom_norb(ia)
            ls0 = bas%atom_lshell_start(ia)
            acc = 0.0_wp
            do mu = o0, o0 + n - 1
                local = mu - o0 + 1
                if (local == 1) then
                    ils = 1
                else if (local <= 4) then
                    ils = 2
                else
                    ils = 3
                end if
                lshell_q(ls0 + ils - 1) = lshell_q(ls0 + ils - 1) + mshell_q(mu)
                acc = acc + mshell_q(mu)
            end do
            q(ia) = acc
        end do
    end subroutine mulliken_population


    !> Wrapper : remplit l'état DFTB (q, dq, lshell_q, lshell_dq,
    !> mshell_q) à partir de P et S.
    subroutine build_charges(P, S, st)
        real(wp),          intent(in)    :: P(:,:), S(:,:)
        type(dftbstate_t), intent(inout) :: st

        call mulliken_population(P, S, st%bas, st%q, st%lshell_q, st%mshell_q)
        call partial_atomic_dq(st%q, st%bas, st%dq)
        call partial_lshell_dq(st%bas, st%lshell_q, st%lshell_dq)
    end subroutine build_charges

end module charges
