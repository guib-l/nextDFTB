!> Driver simple-point : agnostique du calculateur concret.
!>
!> Sélectionne le calculateur concret en fonction de `inp%calc%kind`,
!> puis n'appelle que les méthodes garanties par l'interface abstraite
!> `method_calc_t` (init, build, execute, get_total_energy,
!> get_total_gradient, write_output).
module single_point
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use keywords,      only: input_kw_t
    use property,      only: property_method_t
    use method_calc,   only: method_calc_t
    use dftb,          only: dftb_calc_t
    use dft,           only: dft_calc_t
    use write_output,  only: section, subsection, line, &
                             kv_str, kv_int, kv_real, kv_real_es, kv_log
    use timer,         only: start_timer, stop_timer
    use units,         only: bohr_to_ang
    use errors,        only: fatal
    implicit none
    private

    public :: run_single_point

contains

    subroutine run_single_point(struct, inp)
        type(structure_t), intent(in) :: struct
        type(input_kw_t),  intent(in) :: inp

        class(method_calc_t),     allocatable :: calc
        class(property_method_t), allocatable :: prop
        character(len=128) :: buf
        integer :: i

        call section("Calculation properties")
        call kv_str("driver",     "single_point")
        call kv_str("calculator", trim(inp%calc%kind))

        call subsection("Geometry")
        call kv_int("natoms", struct%natoms)
        call line("    idx  symbol         x (Å)         y (Å)         z (Å)         q0")
        do i = 1, struct%natoms
            write(buf, '(a,i5,3x,a4,3(1x,f14.6),1x,f12.8)') "  ", i, &
                struct%atoms(i)%symbol, &
                bohr_to_ang(struct%atoms(i)%position(1)), &
                bohr_to_ang(struct%atoms(i)%position(2)), &
                bohr_to_ang(struct%atoms(i)%position(3)), &
                struct%atoms(i)%charge
            call line(buf)
        end do

        call subsection("Basis")
        call kv_str("src",  trim(inp%basis%src))
        call kv_str("ext",  trim(inp%basis%ext))
        call kv_str("sep",  trim(inp%basis%sep))
        call kv_str("type", trim(inp%basis%type))

        !-- Sélection du calculateur concret ---------------------------
        select case (trim(inp%calc%kind))
        case ("DFTB")
            allocate(dftb_calc_t :: calc)
            allocate(prop, source=inp%calc%dftb)
        case ("DFT")
            allocate(dft_calc_t :: calc)
            allocate(prop, source=inp%calc%dft)
        case default
            call fatal("single_point", "unsupported calc kind: "//trim(inp%calc%kind))
        end select

        !-- Pipeline générique : init → build → execute → write --------
        call start_timer("CALC_INIT_BUILD")
        call calc%init(struct, inp%basis, prop)
        call calc%build()
        call stop_timer("CALC_INIT_BUILD")

        call start_timer("CALC_EXECUTE")
        call calc%execute()
        call stop_timer("CALC_EXECUTE")

        call calc%write_output()
    end subroutine run_single_point
end module single_point
