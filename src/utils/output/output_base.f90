!> Objet abstrait des modules d'écriture de résultats.
!>
!> Chaque module `write_*.f90` (dftb, dft, skf) expose un objet qui
!> hérite de `output_base_t` et fournit deux méthodes :
!>   - write_process() : écrit l'ensemble du déroulement du calcul
!>   - write_result()  : écrit les résultats essentiels (sortie minimale)
!>
!> Le niveau de verbosité (`verbose`) se lit comme :
!>   1 = peu verbeux ;  2 = standard (défaut) ;  3 = très verbeux.
module output_base
    implicit none
    private

    integer, parameter, public :: VERBOSE_LOW    = 1
    integer, parameter, public :: VERBOSE_NORMAL = 2
    integer, parameter, public :: VERBOSE_HIGH   = 3

    type, abstract, public :: output_base_t
        integer :: verbose = VERBOSE_NORMAL
    contains
        procedure(write_iface), deferred :: write_process
        procedure(write_iface), deferred :: write_result
    end type output_base_t

    abstract interface
        subroutine write_iface(self)
            import :: output_base_t
            class(output_base_t), intent(inout) :: self
        end subroutine write_iface
    end interface
end module output_base
