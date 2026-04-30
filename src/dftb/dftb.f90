!> Façade publique du calculateur DFTB.
!>
!> Expose EXACTEMENT l'interface suivante :
!>   - subroutine init(geometry, basis)
!>   - subroutine execute(dograd)
!>   - function get_total_energy()
!>   - function get_repulsive_energy()
!>   - function get_coulomb_energy()
!>   - function get_band_energy()
!>   - function get_gradient()
!>   - function get_charges()
!>
!> Aucune autre méthode publique. Le module conserve en interne le
!> dernier état du calcul (`dftbstate_t`).
module dftb
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use parse_input,   only: basis_t, calc_t
    use dftbstate,     only: dftbstate_t
    use skf,           only: skf_init       => init,            &
                              skf_readslako  => readslako,       &
                              skf_build_rep  => build_repulsion, &
                              skf_build_elec => build_electronic
    use matel,         only: build_basis_system
    use scc,           only: solve_scc
    use repulsif,      only: repulsive_energy
    use coulomb,       only: coulomb_energy
    use dftb_energy,   only: total_energy
    use dftb_grad,     only: zero_gradient
    use errors,        only: fatal
    implicit none
    private

    type(structure_t),       save :: g_struct
    type(basis_t),           save :: g_basis
    type(calc_t),            save :: g_calc
    type(dftbstate_t),       save :: g_state
    logical, save :: g_init = .false.

    public :: init, execute
    public :: get_total_energy, get_repulsive_energy, get_coulomb_energy
    public :: get_band_energy, get_gradient, get_charges

contains

    subroutine init(geometry, basis_in, calcul)
        type(structure_t), intent(in) :: geometry
        type(basis_t),     intent(in) :: basis_in
        type(calc_t),      intent(in) :: calcul

        g_struct = geometry
        g_basis  = basis_in
        g_calc   = calcul

        call skf_init(geometry, basis_in)
        call skf_readslako()
        call skf_build_rep()
        call skf_build_elec()

        call build_basis_system(geometry, g_basis, g_state%bas)

        if (allocated(g_state%grad)) deallocate(g_state%grad)
        allocate(g_state%grad(3, geometry%natoms))
        g_state%grad = 0.0_wp

        g_init = .true.
    end subroutine init


    subroutine execute(dograd)
        logical, intent(in) :: dograd

        if (.not. g_init) call fatal("dftb", "execute called before init")

        call solve_scc(g_struct, g_calc%scc, g_calc%maxscc, g_calc%tolscc, g_state)

        g_state%e_rep = repulsive_energy(g_struct)
        if (g_calc%scc) then
            g_state%e_coul = coulomb_energy(g_state%gamma, g_state%dq)
        else
            g_state%e_coul = 0.0_wp
        end if
        g_state%e_total = total_energy(g_state%e_band, g_state%e_coul, g_state%e_rep)

        call zero_gradient(g_state%grad)
        if (dograd) then
            ! Gradient analytique non implémenté.
            continue
        end if
    end subroutine execute


    function get_total_energy() result(e)
        real(wp) :: e
        e = g_state%e_total
    end function get_total_energy

    function get_repulsive_energy() result(e)
        real(wp) :: e
        e = g_state%e_rep
    end function get_repulsive_energy

    function get_coulomb_energy() result(e)
        real(wp) :: e
        e = g_state%e_coul
    end function get_coulomb_energy

    function get_band_energy() result(e)
        real(wp) :: e
        e = g_state%e_band
    end function get_band_energy

    function get_gradient() result(g)
        real(wp), allocatable :: g(:,:)
        allocate(g(3, g_struct%natoms))
        g = g_state%grad
    end function get_gradient

    function get_charges() result(q)
        real(wp), allocatable :: q(:)
        allocate(q(g_struct%natoms))
        if (allocated(g_state%q)) then
            q = g_state%q
        else
            q = 0.0_wp
        end if
    end function get_charges
end module dftb
