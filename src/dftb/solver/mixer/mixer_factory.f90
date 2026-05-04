!> Fabrique polymorphe de mixers à partir d'une `property_mixer_t`.
module mixer_factory
    use kinds,         only: wp
    use property,      only: property_mixer_t
    use mixer,         only: mixer_t
    use simple_mixer,  only: simpleMixer_t
    use random_mixer,  only: randomMixer_t
    use broyden_mixer, only: broydenMixer_t
    use errors,        only: fatal
    implicit none
    private

    public :: make_mixer

contains

    subroutine make_mixer(spec, n, mx)
        type(property_mixer_t),       intent(in)  :: spec
        integer,                      intent(in)  :: n
        class(mixer_t), allocatable,  intent(out) :: mx
        character(len=len(spec%type)) :: kind

        kind = upper(trim(spec%type))

        select case (trim(kind))
        case ("SIMPLE")
            allocate(simpleMixer_t :: mx)
        case ("RANDOM")
            allocate(randomMixer_t :: mx)
        case ("BROYDEN")
            allocate(broydenMixer_t :: mx)
        case default
            call fatal("mixer_factory", "unknown MIXING.TYPE: "//trim(spec%type))
        end select

        call mx%init(n, spec%factor, spec%history, spec%omega0)
    end subroutine make_mixer

    pure function upper(s) result(o)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: o
        integer :: i, c
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= iachar('a') .and. c <= iachar('z')) then
                o(i:i) = achar(c - 32)
            else
                o(i:i) = s(i:i)
            end if
        end do
    end function upper
end module mixer_factory
