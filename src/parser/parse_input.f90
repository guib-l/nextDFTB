!> Parser de l'input utilisateur (`input.dat`).
!>
!> Format générique : sections SECTION = { ... } contenant soit
!> des assignations KEY = VALUE, soit des sous-sections SECTION = { ... }.
!> Le module définit aussi les types de configuration (basis, calc,
!> output, input) qui sont remplis à la lecture.
module parse_input
    use kinds,    only: wp
    use constants, only: PATH_LEN, SYMBOL_LEN
    use defaults, only: DEFAULT_GEO, DEFAULT_EXT, DEFAULT_SEP, DEFAULT_SRC, &
                        DEFAULT_TYPE, DEFAULT_OUT, DEFAULT_LOG, DEFAULT_LOG_ON, &
                        DEFAULT_MAXSCC, DEFAULT_TOLSCC, DEFAULT_SCC
    use keywords
    use errors,   only: fatal, warn
    implicit none
    private

    type, public :: orbital_spec_t
        character(len=SYMBOL_LEN) :: symbol = ""
        character(len=64)         :: orbitals = ""    ! ex: "1s 2s"
    end type orbital_spec_t

    type, public :: basis_t
        character(len=PATH_LEN) :: src  = ""
        character(len=16)       :: ext  = ".skf"
        character(len=4)        :: sep  = "-"
        character(len=16)       :: type = "spd"
        type(orbital_spec_t), allocatable :: orbitals(:)
    end type basis_t

    type, public :: calc_t
        character(len=8) :: kind   = "DFTB"   ! DFTB | DFT | SKF
        logical          :: scc    = .false.
        integer          :: maxscc = 100
        real(wp)         :: tolscc = 1.0e-5_wp
    end type calc_t

    type, public :: output_t
        character(len=PATH_LEN) :: out    = "results.out"
        character(len=PATH_LEN) :: log    = "logs"
        logical                 :: log_on = .false.
    end type output_t

    type, public :: input_t
        character(len=PATH_LEN) :: geom_file = "geometry.dat"
        type(basis_t)           :: basis
        type(calc_t)            :: calc
        type(output_t)          :: out
        logical                 :: has_driver  = .false.
        character(len=16)       :: driver_type = "SINGLE"
    end type input_t

    public :: read_input

contains

    subroutine read_input(filename, inp)
        character(len=*), intent(in)  :: filename
        type(input_t),    intent(out) :: inp

        integer :: u, ios
        character(len=512) :: raw
        character(len=:),  allocatable :: line
        character(len=64), allocatable :: stack(:)
        integer :: depth

        ! Defaults
        inp%geom_file  = DEFAULT_GEO
        inp%basis%src  = DEFAULT_SRC
        inp%basis%ext  = DEFAULT_EXT
        inp%basis%sep  = DEFAULT_SEP
        inp%basis%type = DEFAULT_TYPE
        inp%out%out    = DEFAULT_OUT
        inp%out%log    = DEFAULT_LOG
        inp%out%log_on = DEFAULT_LOG_ON
        inp%calc%scc    = DEFAULT_SCC
        inp%calc%maxscc = DEFAULT_MAXSCC
        inp%calc%tolscc = DEFAULT_TOLSCC
        inp%calc%kind   = "DFTB"

        allocate(stack(8))
        depth = 0

        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) call fatal("parse_input", "cannot open: "//trim(filename))

        do
            read(u, '(a)', iostat=ios) raw
            if (ios /= 0) exit
            line = strip_comment(strip(raw))
            if (len_trim(line) == 0) cycle

            if (line == "}") then
                if (depth == 0) call fatal("parse_input", "unexpected '}'")
                depth = depth - 1
                cycle
            end if

            if (is_block_open(line)) then
                depth = depth + 1
                if (depth > size(stack)) call fatal("parse_input", "section nesting too deep")
                stack(depth) = upper(extract_block_name(line))
                cycle
            end if

            if (has_equal(line)) then
                call dispatch_assign(stack, depth, line, inp)
                cycle
            end if

            call warn("parse_input", "ignored line: "//trim(line))
        end do
        close(u)

        if (depth /= 0) call fatal("parse_input", "unclosed '{' at EOF")
        if (len_trim(inp%basis%src) == 0) call fatal("parse_input", "BASIS.SRC missing")

        deallocate(stack)
    end subroutine read_input


    subroutine dispatch_assign(stack, depth, line, inp)
        character(len=*), intent(in)    :: stack(:), line
        integer,          intent(in)    :: depth
        type(input_t),    intent(inout) :: inp
        character(len=:), allocatable :: key, val
        character(len=64) :: ukey

        call split_assign(line, key, val)
        ukey = upper(key)

        if (depth == 0) then
            call warn("parse_input", "top-level assignment ignored: "//trim(key))
            return
        end if

        select case (depth)
        case (1)
            select case (trim(stack(1)))
            case (KW_GEOMETRY)
                if (trim(ukey) == KW_GEO) then
                    inp%geom_file = unquote(val)
                else
                    call warn("parse_input", "unknown GEOMETRY key: "//trim(key))
                end if
            case (KW_BASIS)
                select case (trim(ukey))
                case (KW_SRC);  inp%basis%src  = unquote(val)
                case (KW_EXT);  inp%basis%ext  = unquote(val)
                case (KW_SEP);  inp%basis%sep  = unquote(val)
                case (KW_TYPE); inp%basis%type = unquote(val)
                case default
                    call warn("parse_input", "unknown BASIS key: "//trim(key))
                end select
            case (KW_OUTPUT)
                select case (trim(ukey))
                case (KW_OUT); inp%out%out = unquote(val)
                case (KW_LOG)
                    if (is_bool(val)) then
                        inp%out%log_on = bool_value(val)
                    else
                        inp%out%log    = unquote(val)
                        inp%out%log_on = .true.
                    end if
                case default
                    call warn("parse_input", "unknown OUTPUT key: "//trim(key))
                end select
            case (KW_DRIVER)
                if (trim(ukey) == KW_DRV_TYPE) then
                    inp%driver_type = upper(unquote(val))
                    inp%has_driver  = .true.
                else
                    call warn("parse_input", "unknown DRIVER key: "//trim(key))
                end if
            case default
                call warn("parse_input", "unknown section: "//trim(stack(1)))
            end select

        case (2)
            if (trim(stack(1)) == KW_CALC .and. trim(stack(2)) == KW_DFTB) then
                inp%calc%kind = "DFTB"
                select case (trim(ukey))
                case (KW_SCC);    inp%calc%scc    = bool_value(val)
                case (KW_MAXSCC); read(val, *) inp%calc%maxscc
                case (KW_TOLSCC); read(val, *) inp%calc%tolscc
                case default
                    call warn("parse_input", "unknown CALC.DFTB key: "//trim(key))
                end select
            else if (trim(stack(1)) == KW_BASIS .and. trim(stack(2)) == KW_ORBITALS) then
                call append_orbital(inp%basis, key, val)
            else
                call warn("parse_input", "ignored nested key in "//trim(stack(1))//"/"//trim(stack(2)))
            end if
        case default
            call warn("parse_input", "deep nesting not handled")
        end select
    end subroutine dispatch_assign


    subroutine append_orbital(basis, sym, orbs)
        type(basis_t),    intent(inout) :: basis
        character(len=*), intent(in)    :: sym, orbs
        type(orbital_spec_t), allocatable :: tmp(:)
        integer :: n
        n = 0
        if (allocated(basis%orbitals)) n = size(basis%orbitals)
        allocate(tmp(n + 1))
        if (n > 0) tmp(1:n) = basis%orbitals(1:n)
        tmp(n+1)%symbol   = sym
        tmp(n+1)%orbitals = orbs
        call move_alloc(tmp, basis%orbitals)
    end subroutine append_orbital


    !-- helpers ----------------------------------------------------------

    pure function strip(s) result(o)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: o
        integer :: i, j
        i = 1
        do while (i <= len(s) .and. (s(i:i) == ' ' .or. s(i:i) == char(9)))
            i = i + 1
        end do
        j = len(s)
        do while (j >= i .and. (s(j:j) == ' ' .or. s(j:j) == char(9) .or. s(j:j) == char(13)))
            j = j - 1
        end do
        if (j >= i) then
            o = s(i:j)
        else
            o = ""
        end if
    end function strip

    pure function strip_comment(s) result(o)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: o
        integer :: ih, ie
        ih = index(s, '#')
        ie = index(s, '!')
        if (ih == 0) ih = len(s) + 1
        if (ie == 0) ie = len(s) + 1
        o = strip(s(1:min(ih, ie) - 1))
    end function strip_comment

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

    pure function has_equal(s) result(b)
        character(len=*), intent(in) :: s
        logical :: b
        b = index(s, '=') > 0
    end function has_equal

    pure function is_block_open(s) result(b)
        character(len=*), intent(in) :: s
        logical :: b
        b = (index(s, '=') > 0) .and. (index(s, '{') > 0)
    end function is_block_open

    pure function extract_block_name(s) result(name)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: name
        integer :: ie
        ie = index(s, '=')
        name = strip(s(1:ie-1))
    end function extract_block_name

    subroutine split_assign(s, key, val)
        character(len=*), intent(in) :: s
        character(len=:), allocatable, intent(out) :: key, val
        integer :: ie
        ie = index(s, '=')
        key = strip(s(1:ie-1))
        val = strip(s(ie+1:))
    end subroutine split_assign

    pure function unquote(s) result(o)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: o
        integer :: n
        n = len_trim(s)
        if (n >= 2) then
            if ((s(1:1) == '"' .and. s(n:n) == '"') .or. &
                (s(1:1) == "'" .and. s(n:n) == "'")) then
                o = s(2:n-1)
                return
            end if
        end if
        o = trim(s)
    end function unquote

    pure function is_bool(s) result(b)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: u
        logical :: b
        u = upper(trim(s))
        b = (trim(u) == "TRUE" .or. trim(u) == "FALSE" .or. &
             trim(u) == "T"    .or. trim(u) == "F"     .or. &
             trim(u) == ".TRUE." .or. trim(u) == ".FALSE.")
    end function is_bool

    pure function bool_value(s) result(b)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: u
        logical :: b
        u = upper(trim(s))
        b = (trim(u) == "TRUE" .or. trim(u) == "T" .or. trim(u) == ".TRUE.")
    end function bool_value
end module parse_input
