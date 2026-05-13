!> Interface abstraite des calculateurs de base (SKF, ...).
!>
!> Tout calculateur de base doit implémenter au minimum :
!>   - init(struct, basis)  : initialise le calcul
!>   - build()              : construit les objets nécessaires
!>   - execute()            : exécute le calcul
!>   - write_output()       : écriture des résultats
!>
!> Les calculateurs concrets peuvent exposer des méthodes additionnelles
!> spécifiques (accès aux intégrales, à la table H/S, etc.).
module method_basis
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use property,      only: property_basis_t
    implicit none
    private

    type, abstract, public :: method_basis_t
    contains
        procedure(init_iface),     deferred :: init
        procedure(build_iface),    deferred :: build
        procedure(execute_iface),  deferred :: execute
        procedure(write_iface),    deferred :: write_output
    end type method_basis_t

    abstract interface
        subroutine init_iface(self, struct, basis)
            import :: method_basis_t, structure_t, property_basis_t
            class(method_basis_t),  intent(inout) :: self
            type(structure_t),      intent(in)    :: struct
            type(property_basis_t), intent(in)    :: basis
        end subroutine init_iface

        subroutine build_iface(self)
            import :: method_basis_t
            class(method_basis_t), intent(inout) :: self
        end subroutine build_iface

        subroutine execute_iface(self)
            import :: method_basis_t
            class(method_basis_t), intent(inout) :: self
        end subroutine execute_iface

        subroutine write_iface(self)
            import :: method_basis_t
            class(method_basis_t), intent(inout) :: self
        end subroutine write_iface
    end interface
end module method_basis
