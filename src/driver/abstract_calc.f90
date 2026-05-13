!> Interface abstraite des calculateurs de méthode (DFTB, DFT).
!>
!> Tout calculateur de méthode doit implémenter au minimum :
!>   - init(struct, basis, method)  : initialise le calcul
!>   - build()                      : construit les objets nécessaires
!>   - execute()                    : exécute le calcul
!>   - get_total_energy()           : énergie totale
!>   - get_total_gradient()         : gradient
!>   - write_output()               : écriture des résultats
!>
!> Les calculateurs concrets peuvent exposer des méthodes additionnelles
!> spécifiques (énergies décomposées, charges, etc.).
module method_calc
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use property,      only: property_basis_t, property_method_t
    implicit none
    private

    type, abstract, public :: method_calc_t
    contains
        procedure(init_iface),     deferred :: init
        procedure(build_iface),    deferred :: build
        procedure(execute_iface),  deferred :: execute
        procedure(energy_iface),   deferred :: get_total_energy
        procedure(grad_iface),     deferred :: get_total_gradient
        procedure(write_iface),    deferred :: write_output
    end type method_calc_t

    abstract interface
        subroutine init_iface(self, struct, basis, method)
            import :: method_calc_t, structure_t, property_basis_t, property_method_t
            class(method_calc_t),     intent(inout) :: self
            type(structure_t),        intent(in)    :: struct
            type(property_basis_t),   intent(in)    :: basis
            class(property_method_t), intent(in)    :: method
        end subroutine init_iface

        subroutine build_iface(self)
            import :: method_calc_t
            class(method_calc_t), intent(inout) :: self
        end subroutine build_iface

        subroutine execute_iface(self)
            import :: method_calc_t
            class(method_calc_t), intent(inout) :: self
        end subroutine execute_iface

        function energy_iface(self) result(e)
            import :: method_calc_t, wp
            class(method_calc_t), intent(in) :: self
            real(wp) :: e
        end function energy_iface

        function grad_iface(self) result(g)
            import :: method_calc_t, wp
            class(method_calc_t), intent(in) :: self
            real(wp), allocatable :: g(:,:)
        end function grad_iface

        subroutine write_iface(self)
            import :: method_calc_t
            class(method_calc_t), intent(inout) :: self
        end subroutine write_iface
    end interface
end module method_calc
