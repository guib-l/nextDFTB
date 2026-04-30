!> Gradient DFTB.
!>
!> Le gradient analytique requiert les dérivées de H, S, γ et de la
!> répulsion. Pour rester minimal, ce module se contente de fournir un
!> placeholder retournant un gradient nul (à substituer ultérieurement).
module dftb_grad
    use kinds, only: wp
    implicit none
    private

    public :: zero_gradient

contains

    subroutine zero_gradient(grad)
        real(wp), intent(out) :: grad(:,:)
        grad = 0.0_wp
    end subroutine zero_gradient
end module dftb_grad
