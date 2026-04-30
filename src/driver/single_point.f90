!> Driver simple-point : appelle le calculateur DFTB en passant par
!> sa seule API publique. Écrit le rappel des propriétés du calcul,
!> laisse le calculateur produire sa sortie SCF, puis écrit le résultat
!> final via write_dftb. Les phases sont chronométrées.
module single_point
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use parse_input,   only: input_t
    use dftb,          only: dftb_init       => init,             &
                              dftb_exec       => execute,          &
                              get_total_energy, get_repulsive_energy, &
                              get_coulomb_energy, get_band_energy, &
                              get_charges
    use write_output,  only: section, subsection, line, &
                             kv_str, kv_int, kv_real, kv_real_es, kv_log
    use write_dftb,    only: write_dftb_final
    use timer,         only: timer_t, tic, toc, timer_record
    use units,         only: bohr_to_ang
    implicit none
    private

    public :: run_single_point

contains

    subroutine run_single_point(struct, inp)
        type(structure_t), intent(in) :: struct
        type(input_t),     intent(in) :: inp

        real(wp) :: e_tot, e_band, e_coul, e_rep
        real(wp), allocatable :: q(:)
        type(timer_t) :: t_init, t_exec
        character(len=128) :: buf
        integer :: i

        !-- Rappel des propriétés du calcul ----------------------------
        call section("Calculation properties")
        call kv_str("driver",    "single_point")
        call kv_str("calculator", trim(inp%calc%kind))
        call kv_log("scc",       inp%calc%scc)
        if (inp%calc%scc) then
            call kv_int("max_scc",   inp%calc%maxscc)
            call kv_real_es("tol_scc", inp%calc%tolscc)
        end if

        call subsection("Geometry")
        call kv_int("natoms", struct%natoms)
        call line("    idx  symbol         x (Å)         y (Å)         z (Å)     q0")
        do i = 1, struct%natoms
            write(buf, '(a,i5,3x,a4,3(1x,f14.6),1x,f8.3)') "  ", i, &
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

        !-- Initialisation du calculateur ------------------------------
        call tic(t_init)
        call dftb_init(struct, inp%basis, inp%calc)
        call toc(t_init)
        call timer_record("dftb_init", t_init)

        !-- Exécution (le calculateur écrit la boucle SCF) -------------
        call tic(t_exec)
        call dftb_exec(.false.)
        call toc(t_exec)
        call timer_record("dftb_execute", t_exec)

        e_tot  = get_total_energy()
        e_band = get_band_energy()
        e_coul = get_coulomb_energy()
        e_rep  = get_repulsive_energy()
        q      = get_charges()

        call write_dftb_final(e_tot, e_band, e_coul, e_rep, q)
    end subroutine run_single_point
end module single_point
