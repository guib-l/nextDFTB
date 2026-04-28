!> nextDFTB — programme principal.
!>
!> Pipeline :
!>   main → cli → parse_input → driver → dftb → write_output
program nextdftb
    use globals,     only: input_t
    use cli,         only: parse_cli
    use parse_input, only: read_input
    use driver,      only: run_default
    implicit none

    character(len=:), allocatable :: input_file
    type(input_t) :: inp

    call parse_cli(input_file)
    call read_input(input_file, inp)
    call run_default(inp)
end program nextdftb
