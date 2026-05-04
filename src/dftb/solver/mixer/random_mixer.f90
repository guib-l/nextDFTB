!> Mixer aléatoire : remplace x par un vecteur uniforme dans [-factor, factor].
module random_mixer
    use kinds, only: wp
    use mixer, only: mixer_t
    implicit none
    private

    type, extends(mixer_t), public :: randomMixer_t
    contains
        procedure :: init => random_init
        procedure :: mix  => random_mix
        procedure :: free => random_free
    end type randomMixer_t

contains

    subroutine random_init(self, n, factor, history, omega0)
        class(randomMixer_t), intent(inout) :: self
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
    end subroutine random_init

    subroutine random_mix(self, x_in, x_out, x_new)
        class(randomMixer_t), intent(inout) :: self
        real(wp),             intent(in)    :: x_in(:)
        real(wp),             intent(in)    :: x_out(:)
        real(wp),             intent(out)   :: x_new(:)
        real(wp) :: r(size(x_new))
        real(wp) :: dummy
        dummy = sum(x_in) + sum(x_out)
        call random_number(r)
        x_new = (2.0_wp * r - 1.0_wp) * self%factor
        self%it = self%it + 1
    end subroutine random_mix

    subroutine random_free(self)
        class(randomMixer_t), intent(inout) :: self
        self%n  = 0
        self%it = 0
    end subroutine random_free
end module random_mixer
