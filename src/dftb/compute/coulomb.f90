!> Énergie coulombienne 2nd ordre DFTB :
!>      E_coul = (1/2) Σ_{A,B} γ_AB Δq_A Δq_B
module coulomb
    use kinds, only: wp
    implicit none
    private

    public :: coulomb_energy, coulomb_potential

contains

    function coulomb_energy(gamma, dq) result(e)
        real(wp), intent(in) :: gamma(:,:), dq(:)
        real(wp) :: e
        integer  :: i, j, n
        n = size(dq)
        e = 0.0_wp
        do i = 1, n
            do j = 1, n
                e = e + gamma(i, j) * dq(i) * dq(j)
            end do
        end do
        e = 0.5_wp * e
    end function coulomb_energy

    !> Potentiel atomique : V_A = Σ_B γ_AB Δq_B
    subroutine coulomb_potential(gamma, dq, V)
        real(wp), intent(in)  :: gamma(:,:), dq(:)
        real(wp), intent(out) :: V(:)
        integer :: i, j, n
        n = size(dq)
        do i = 1, n
            V(i) = 0.0_wp
            do j = 1, n
                V(i) = V(i) + gamma(i, j) * dq(j)
            end do
        end do
    end subroutine coulomb_potential
end module coulomb
