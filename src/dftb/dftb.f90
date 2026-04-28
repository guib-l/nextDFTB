!> nextDFTB — façade publique du module de calcul DFTB.
!>
!> Interface publique (et SEULE) :
!>   - subroutine init(geom, basis, calc)
!>   - subroutine execute(dograd)
!>   - function get_total_energy()
!>   - function get_repulsive_energy()
!>   - function get_coulomb_energy()
!>   - function get_band_energy()
!>   - function get_gradient()
!>   - function get_charges()
module dftb
    use kinds,        only: wp
    use globals,      only: geometry_t, basis_t, calc_t, basis_system_t
    use readskf,      only: skf_store_t, load_skf_store
    use dftb_matrix,  only: build_basis_system
    use scc,          only: scc_result_t, solve_scc
    use repulsif,     only: repulsive_energy
    use coulomb,      only: coulomb_energy
    use dftb_grad,    only: zero_gradient
    use errors,       only: fatal
    implicit none
    private

    type(geometry_t),     save :: g_geom
    type(basis_t),        save :: g_basis_in
    type(calc_t),         save :: g_calc
    type(basis_system_t), save :: g_bas
    type(skf_store_t),    save :: g_store
    type(scc_result_t),   save :: g_res
    real(wp), allocatable, save :: g_grad(:,:)
    real(wp), save :: g_e_total = 0.0_wp
    real(wp), save :: g_e_band  = 0.0_wp
    real(wp), save :: g_e_coul  = 0.0_wp
    real(wp), save :: g_e_rep   = 0.0_wp
    logical, save  :: g_init    = .false.

    public :: init, execute
    public :: get_total_energy, get_repulsive_energy, get_coulomb_energy
    public :: get_band_energy, get_gradient, get_charges

contains

    subroutine init(geom, basis, calcul)
        type(geometry_t), intent(in) :: geom
        type(basis_t),    intent(in) :: basis
        type(calc_t),     intent(in) :: calcul

        character(len=4), allocatable :: sym_table(:)

        g_geom     = geom
        g_basis_in = basis
        g_calc     = calcul

        call unique_symbols(geom%symbols, sym_table)
        call load_skf_store(sym_table, basis%src, basis%ext, basis%sep, g_store)
        call build_basis_system(geom, g_store, sym_table, g_bas)

        if (allocated(g_grad)) deallocate(g_grad)
        allocate(g_grad(3, geom%natoms))
        g_grad = 0.0_wp

        g_init = .true.
    end subroutine init


    subroutine execute(dograd)
        logical, intent(in) :: dograd

        if (.not. g_init) call fatal("dftb", "execute() called before init()")

        call solve_scc(g_geom, g_bas, g_store, g_calc%scc, &
                       g_calc%maxscc, g_calc%tolscc, g_res)

        g_e_band = g_res%e_band
        g_e_rep  = repulsive_energy(g_geom, g_bas, g_store)
        if (g_calc%scc) then
            g_e_coul = coulomb_energy(g_res%gamma, g_res%dq)
            g_e_total = g_e_band - g_e_coul + g_e_rep
        else
            g_e_coul = 0.0_wp
            g_e_total = g_e_band + g_e_rep
        end if

        if (dograd) then
            ! Gradient analytique non implémenté — on retourne 0.
            call zero_gradient(g_grad)
        else
            call zero_gradient(g_grad)
        end if
    end subroutine execute


    function get_total_energy() result(e)
        real(wp) :: e
        e = g_e_total
    end function get_total_energy

    function get_repulsive_energy() result(e)
        real(wp) :: e
        e = g_e_rep
    end function get_repulsive_energy

    function get_coulomb_energy() result(e)
        real(wp) :: e
        e = g_e_coul
    end function get_coulomb_energy

    function get_band_energy() result(e)
        real(wp) :: e
        e = g_e_band
    end function get_band_energy

    function get_gradient() result(g)
        real(wp), allocatable :: g(:,:)
        allocate(g(3, g_geom%natoms))
        g = g_grad
    end function get_gradient

    function get_charges() result(q)
        real(wp), allocatable :: q(:)
        allocate(q(g_geom%natoms))
        if (allocated(g_res%q)) then
            q = g_res%q
        else
            q = 0.0_wp
        end if
    end function get_charges


    !-- helpers --------------------------------------------------------

    subroutine unique_symbols(syms, uniq)
        character(len=*),               intent(in)  :: syms(:)
        character(len=4), allocatable,  intent(out) :: uniq(:)
        integer :: i, j, n
        character(len=4) :: tmp(size(syms))
        logical :: seen
        n = 0
        do i = 1, size(syms)
            seen = .false.
            do j = 1, n
                if (trim(tmp(j)) == trim(syms(i))) then
                    seen = .true.; exit
                end if
            end do
            if (.not. seen) then
                n = n + 1
                tmp(n) = syms(i)
            end if
        end do
        allocate(uniq(n))
        uniq(1:n) = tmp(1:n)
    end subroutine unique_symbols
end module dftb
