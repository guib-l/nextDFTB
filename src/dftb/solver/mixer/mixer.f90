!> Interface abstraite pour les mixers de charges SCC.
!>
!> Toute implémentation concrète étend `mixer_t` et fournit `init`,
!> `mix` et `free`. La construction polymorphe se fait via
!> `mixer_factory`.
module mixer
    use kinds, only: wp
    implicit none
    private

    type, abstract, public :: mixer_t
        integer  :: n       = 0
        integer  :: it      = 0
        real(wp) :: factor  = 0.1_wp
    contains
        procedure(mixer_init), deferred :: init
        procedure(mixer_mix),  deferred :: mix
        procedure(mixer_free), deferred :: free
    end type mixer_t

    abstract interface
        subroutine mixer_init(self, n, factor, history, omega0)
            import :: mixer_t, wp
            class(mixer_t), intent(inout) :: self
            integer,        intent(in)    :: n
            real(wp),       intent(in)    :: factor
            integer,        intent(in)    :: history
            real(wp),       intent(in)    :: omega0
        end subroutine mixer_init

        subroutine mixer_mix(self, x_in, x_out, x_new)
            import :: mixer_t, wp
            class(mixer_t), intent(inout) :: self
            real(wp),       intent(in)    :: x_in(:)
            real(wp),       intent(in)    :: x_out(:)
            real(wp),       intent(out)   :: x_new(:)
        end subroutine mixer_mix

        subroutine mixer_free(self)
            import :: mixer_t
            class(mixer_t), intent(inout) :: self
        end subroutine mixer_free
    end interface
end module mixer
