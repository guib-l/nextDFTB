!> CLI minimaliste : récupère le nom du fichier d'input.
module cli
    implicit none
    private

    public :: parse_cli

contains

    subroutine parse_cli(input_file)
        character(len=:), allocatable, intent(out) :: input_file
        integer :: nargs, l
        character(len=512) :: arg

        nargs = command_argument_count()
        if (nargs < 1) then
            input_file = "input.dat"
        else
            call get_command_argument(1, arg, length=l)
            input_file = trim(arg(1:l))
        end if
    end subroutine parse_cli
end module cli
