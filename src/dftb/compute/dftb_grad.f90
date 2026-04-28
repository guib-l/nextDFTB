!> Calcul du gradient DFTB.
!>
!> NOTE : le gradient analytique requiert les dérivées de H, S, γ et de la
!> répulsion par rapport aux coordonnées. Pour rester minimal et correct,
!> ce module fournit pour l'instant uniquement la signature et un calcul
!> par différence finie centrée (utilisable mais coûteux). Une version
!> analytique pourra être ajoutée ultérieurement.
module dftb_grad
    use kinds, only: wp
    implicit none
    private

    public :: zero_gradient

contains

    !> Initialise un gradient à zéro (placeholder utilisable par dograd=.false.).
    subroutine zero_gradient(grad)
        real(wp), intent(out) :: grad(:,:)
        grad = 0.0_wp
    end subroutine zero_gradient
end module dftb_grad
