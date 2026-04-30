!> nextDFTB — programme principal.
!>
!> Pipeline :
!>   input.dat (+ geometry.dat)
!>     → cli → parse_input → parse_geometry
!>     → driver → calculateur (DFTB/SKF/DFT)
!>     → write_output
program nextdftb
    use cli,            only: parse_cli
    use parse_input,    only: input_t, read_input
    use parse_geometry, only: read_geometry
    use structure_mod,  only: structure_t
    use driver,         only: run_default
    implicit none

    character(len=:), allocatable :: input_file
    type(input_t)     :: inp
    type(structure_t) :: struct

    call parse_cli(input_file)
    call read_input(input_file, inp)
    call read_geometry(trim(inp%geom_file), struct)
    call run_default(struct, inp)
end program nextdftb
