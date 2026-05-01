!> Calculateur DFT (squelette).
!>
!> Implémente strictement l'interface abstraite `method_calc_t`. Aucune
!> méthode additionnelle. Les corps des routines lèvent `fatal` tant
!> que la méthode n'est pas implémentée.
module dft
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use property,      only: property_basis_t, property_dft_t, &
                              property_method_t
    use method_calc,   only: method_calc_t
    use errors,        only: fatal
    implicit none
    private

    type, extends(method_calc_t), public :: dft_calc_t
        type(structure_t)      :: struct
        type(property_basis_t) :: basis
        type(property_dft_t)   :: prop
        logical                :: initialized = .false.
    contains
        procedure :: init               => dft_init
        procedure :: build              => dft_build
        procedure :: execute            => dft_execute
        procedure :: get_total_energy   => dft_get_total_energy
        procedure :: get_total_gradient => dft_get_total_gradient
        procedure :: write_output       => dft_write_output
    end type dft_calc_t

contains

    subroutine dft_init(self, struct, basis, method)
        class(dft_calc_t),        intent(inout) :: self
        type(structure_t),        intent(in)    :: struct
        type(property_basis_t),   intent(in)    :: basis
        class(property_method_t), intent(in)    :: method
        self%struct = struct
        self%basis  = basis
        select type (method)
        type is (property_dft_t)
            self%prop = method
        class default
            call fatal("dft", "init: expected property_dft_t")
        end select
        self%initialized = .true.
    end subroutine dft_init

    subroutine dft_build(self)
        class(dft_calc_t), intent(inout) :: self
        call fatal("dft", "build: not implemented")
    end subroutine dft_build

    subroutine dft_execute(self)
        class(dft_calc_t), intent(inout) :: self
        call fatal("dft", "execute: not implemented")
    end subroutine dft_execute

    function dft_get_total_energy(self) result(e)
        class(dft_calc_t), intent(in) :: self
        real(wp) :: e
        e = 0.0_wp
        call fatal("dft", "get_total_energy: not implemented")
    end function dft_get_total_energy

    function dft_get_total_gradient(self) result(g)
        class(dft_calc_t), intent(in) :: self
        real(wp), allocatable :: g(:,:)
        allocate(g(3, self%struct%natoms))
        g = 0.0_wp
        call fatal("dft", "get_total_gradient: not implemented")
    end function dft_get_total_gradient

    subroutine dft_write_output(self)
        class(dft_calc_t), intent(inout) :: self
        call fatal("dft", "write_output: not implemented")
    end subroutine dft_write_output
end module dft
