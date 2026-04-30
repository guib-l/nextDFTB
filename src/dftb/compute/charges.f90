!> Calcul des charges Mulliken par atome.
!>
!>   q_A = Σ_{μ ∈ A} Σ_ν P_{μν} S_{νμ}
module charges
    use kinds,     only: wp
    use dftbstate, only: basis_system_t
    implicit none
    private

    public :: mulliken_charges, delta_charges

contains

    subroutine mulliken_charges(P, S, bas, q)
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
    end subroutine mulliken_charges

    subroutine delta_charges(q, bas, dq)
        real(wp),             intent(in)  :: q(:)
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: dq(:)
        integer :: ia, ie
        do ia = 1, size(q)
            ie = bas%atom_elem(ia)
            dq(ia) = q(ia) - bas%elems(ie)%q_neutral
        end do
    end subroutine delta_charges
end module charges
