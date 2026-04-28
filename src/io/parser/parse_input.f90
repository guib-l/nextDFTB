!> Parser d'input nextDFTB.
!>
!> Format générique : sections SECTION = { ... } contenant soit
!> des assignations KEY = VALUE, soit des sous-sections SECTION = { ... }.
!> GEOMETRY contient en plus des lignes atomiques :
!>     SYMBOL X Y Z [GROUP]
!>
!> Coordonnées attendues en angström (converties en bohr en interne).
module parse_input
    use kinds,    only: wp
    use defaults, only: DEFAULT_EXT, DEFAULT_TYPE, DEFAULT_OUT, DEFAULT_LOG, &
                        DEFAULT_MAXSCC, DEFAULT_TOLSCC, DEFAULT_SCC, DEFAULT_LOG_ON
    use globals,  only: input_t, geometry_t, SYMBOL_LEN
    use keywords
    use units,    only: ang_to_bohr
    use errors,   only: fatal, warn
    implicit none
    private

    public :: read_input

contains

    subroutine read_input(filename, inp)
        character(len=*), intent(in)  :: filename
        type(input_t),    intent(out) :: inp

        integer :: u, ios
        character(len=512) :: raw
        character(len=:), allocatable :: line
        character(len=64), allocatable :: stack(:)
        integer :: depth
        integer :: iatom

        ! Defaults
        inp%basis%ext  = DEFAULT_EXT
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
        iatom = 0

        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) call fatal("parse_input", "cannot open: "//trim(filename))

        do
            read(u, '(a)', iostat=ios) raw
            if (ios /= 0) exit
            line = strip(raw)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle    ! comment

            ! close brace
            if (line == "}") then
                if (depth == 0) call fatal("parse_input", "unexpected '}'")
                depth = depth - 1
                cycle
            end if

            ! open block: NAME = {
            if (is_block_open(line)) then
                depth = depth + 1
                if (depth > size(stack)) call fatal("parse_input", "section nesting too deep")
                stack(depth) = upper(extract_block_name(line))
                if (stack(depth) == KW_GEOMETRY) iatom = 0
                cycle
            end if

            ! inside GEOMETRY section : detect atom lines
            if (depth >= 1) then
                if (stack(depth) == KW_GEOMETRY .and. .not. has_equal(line)) then
                    call parse_atom_line(line, inp%geom, iatom)
                    cycle
                end if
            end if

            ! assignment : KEY = VALUE
            if (has_equal(line)) then
                call dispatch_assign(stack, depth, line, inp)
                cycle
            end if

            call warn("parse_input", "ignored line: "//trim(line))
        end do
        close(u)

        if (depth /= 0) call fatal("parse_input", "unclosed '{' at EOF")
        if (.not. allocated(inp%geom%symbols)) call fatal("parse_input", "GEOMETRY missing")
        if (inp%geom%natoms == 0)              call fatal("parse_input", "NATOMS not set")
        if (iatom /= inp%geom%natoms)          call fatal("parse_input", "atom count != NATOMS")
        if (len_trim(inp%basis%src) == 0)      call fatal("parse_input", "BASIS.SRC missing")

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

        select case (depth)
        case (1)
            select case (trim(stack(1)))
            case (KW_GEOMETRY)
                if (trim(ukey) == KW_NATOMS) then
                    read(val, *) inp%geom%natoms
                    allocate(inp%geom%symbols(inp%geom%natoms))
                    allocate(inp%geom%coords(3, inp%geom%natoms))
                    allocate(inp%geom%groups(inp%geom%natoms))
                    inp%geom%groups = 1
                else
                    call warn("parse_input", "unknown GEOMETRY key: "//trim(key))
                end if
            case (KW_BASIS)
                select case (trim(ukey))
                case (KW_SRC);  inp%basis%src  = unquote(val)
                case (KW_EXT);  inp%basis%ext  = unquote(val)
                case (KW_SEP);  inp%basis%sep  = unquote(val)
                case (KW_TYPE); inp%basis%type = upper(unquote(val))
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
                    call warn("parse_input", "unknown DFTB key: "//trim(key))
                end select
            else
                call warn("parse_input", "ignored nested key in "//trim(stack(1))//"/"//trim(stack(2)))
            end if
        case default
            call warn("parse_input", "deep nesting not handled")
        end select
    end subroutine dispatch_assign


    subroutine parse_atom_line(line, geom, iatom)
        character(len=*),  intent(in)    :: line
        type(geometry_t),  intent(inout) :: geom
        integer,           intent(inout) :: iatom

        character(len=SYMBOL_LEN) :: sym
        real(wp) :: x, y, z
        integer  :: grp, ios

        if (.not. allocated(geom%symbols)) &
            call fatal("parse_input", "atom line before NATOMS")

        iatom = iatom + 1
        if (iatom > geom%natoms) call fatal("parse_input", "too many atoms")

        grp = 1
        read(line, *, iostat=ios) sym, x, y, z, grp
        if (ios /= 0) then
            read(line, *, iostat=ios) sym, x, y, z
            if (ios /= 0) call fatal("parse_input", "bad atom line: "//trim(line))
            grp = 1
        end if

        geom%symbols(iatom)   = sym
        geom%coords(1, iatom) = ang_to_bohr(x)
        geom%coords(2, iatom) = ang_to_bohr(y)
        geom%coords(3, iatom) = ang_to_bohr(z)
        geom%groups(iatom)    = grp
    end subroutine parse_atom_line


    !-- helpers ---------------------------------------------------------

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
