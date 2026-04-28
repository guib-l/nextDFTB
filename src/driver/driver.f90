!> Driver par défaut : effectue un calcul DFTB simple-point et écrit la sortie.
!>
!> Le driver n'utilise QUE les méthodes publiques de dftb.
module driver
    use kinds,         only: wp
    use globals,       only: input_t
    use dftb,          only: dftb_init   => init,    &
                              dftb_exec   => execute, &
                              get_total_energy, get_repulsive_energy, &
                              get_coulomb_energy, get_band_energy,    &
                              get_charges
    use write_output,  only: open_output, close_output, banner, section, line
    use dftb_results,  only: write_dftb_summary
    use logger,        only: log_open, log_close
    implicit none
    private

    public :: run_default

contains

    subroutine run_default(inp)
        type(input_t), intent(in) :: inp
        real(wp) :: e_tot, e_band, e_coul, e_rep
        real(wp), allocatable :: q(:)
        character(len=128) :: buf

        if (inp%out%log_on) call log_open(inp%out%log)
        call open_output(inp%out%out)
        call banner("nextDFTB — single point")

        call section("input summary")
        write(buf, '(a,i0)') "  natoms = ", inp%geom%natoms; call line(buf)
        write(buf, '(a,a)')  "  basis  = ", trim(inp%basis%src); call line(buf)
        write(buf, '(a,l1)') "  scc    = ", inp%calc%scc; call line(buf)

        call dftb_init(inp%geom, inp%basis, inp%calc)
        call dftb_exec(.false.)

        e_tot  = get_total_energy()
        e_band = get_band_energy()
        e_coul = get_coulomb_energy()
        e_rep  = get_repulsive_energy()
        q      = get_charges()

        call write_dftb_summary(e_tot, e_band, e_coul, e_rep, q)

        call close_output()
        call log_close()
    end subroutine run_default
end module driver
