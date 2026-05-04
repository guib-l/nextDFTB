!> Mixer linéaire simple : x_new = x_in + factor * (x_out - x_in).
module simple_mixer
    use kinds, only: wp
    use mixer, only: mixer_t
    implicit none
    private

    type, extends(mixer_t), public :: simpleMixer_t
    contains
        procedure :: init => simple_init
        procedure :: mix  => simple_mix
        procedure :: free => simple_free
    end type simpleMixer_t

contains

    subroutine simple_init(self, n, factor, history, omega0)
        class(simpleMixer_t), intent(inout) :: self
        integer,              intent(in)    :: n
        real(wp),             intent(in)    :: factor
        integer,              intent(in)    :: history
        real(wp),             intent(in)    :: omega0
        integer  :: dummy_h
        real(wp) :: dummy_w
        dummy_h = history
        dummy_w = omega0
        self%n      = n
        self%factor = factor
        self%it     = 0
    end subroutine simple_init

    subroutine simple_mix(self, x_in, x_out, x_new)
        class(simpleMixer_t), intent(inout) :: self
        real(wp),             intent(in)    :: x_in(:)
        real(wp),             intent(in)    :: x_out(:)
        real(wp),             intent(out)   :: x_new(:)
        x_new = x_in + self%factor * (x_out - x_in)
        self%it = self%it + 1
    end subroutine simple_mix

    subroutine simple_free(self)
        class(simpleMixer_t), intent(inout) :: self
        self%n  = 0
        self%it = 0
    end subroutine simple_free
end module simple_mixer
