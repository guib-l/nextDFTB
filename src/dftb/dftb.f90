!> Calculateur DFTB.
!>
!> Implémente l'interface abstraite `method_calc_t` :
!>   - init(struct, basis, method)  : stockage des entrées
!>   - build()                      : chargement SKF + base + allocation
!>   - execute()                    : SCC + énergies
!>   - get_total_energy()           : énergie totale
!>   - get_total_gradient()         : gradient (3, natoms)
!>   - write_output()               : écriture du résultat final
!>
!> Méthodes additionnelles spécifiques :
!>   get_band_energy, get_coulomb_energy, get_repulsive_energy,
!>   get_charges.
module dftb
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use property,      only: property_basis_t, property_dftb_t,        &
                              property_method_t,                        &
                              SCHEME_BASIC, SCHEME_NOSCC, SCHEME_SCC
    use method_calc,   only: method_calc_t
    use dftbstate,     only: dftbstate_t
    use skf,           only: skf_init       => init,            &
                              skf_readslako  => readslako,       &
                              skf_build_rep  => build_repulsion, &
                              skf_build_elec => build_electronic, &
                              skf_get_mass   => get_mass
    use matel,         only: build_basis_system
    use scc,           only: solve_scc
    use repulsif,      only: repulsive_energy
    use coulomb,       only: coulomb_energy
    use dftb_energy,   only: total_energy
    use dftb_grad,     only: zero_gradient
    use write_dftb,    only: write_dftb_final
    use errors,        only: fatal
    implicit none
    private

    type, extends(method_calc_t), public :: dftb_calc_t
        type(structure_t)      :: struct
        type(property_basis_t) :: basis
        type(property_dftb_t)  :: prop
        type(dftbstate_t)      :: state
        logical                :: initialized = .false.
        logical                :: built       = .false.
    contains
        procedure :: init               => dftb_init
        procedure :: build              => dftb_build
        procedure :: execute            => dftb_execute
        procedure :: get_total_energy   => dftb_get_total_energy
        procedure :: get_total_gradient => dftb_get_total_gradient
        procedure :: write_output       => dftb_write_output
        ! méthodes spécifiques DFTB
        procedure :: get_band_energy
        procedure :: get_coulomb_energy
        procedure :: get_repulsive_energy
        procedure :: get_charges
    end type dftb_calc_t

contains

    subroutine dftb_init(self, struct, basis, method)
        class(dftb_calc_t),       intent(inout) :: self
        type(structure_t),        intent(in)    :: struct
        type(property_basis_t),   intent(in)    :: basis
        class(property_method_t), intent(in)    :: method

        self%struct = struct
        self%basis  = basis
        select type (method)
        type is (property_dftb_t)
            self%prop = method
        class default
            call fatal("dftb", "init: expected property_dftb_t")
        end select
        call self%struct%initialize()
        self%initialized = .true.
        self%built       = .false.
    end subroutine dftb_init


    subroutine dftb_build(self)
        class(dftb_calc_t), intent(inout) :: self
        integer :: i

        if (.not. self%initialized) call fatal("dftb", "build before init")

        call skf_init(self%struct, self%basis)
        call skf_readslako()
        call skf_build_rep()
        call skf_build_elec()

        ! Masses atomiques lues depuis les fichiers SKF homonucléaires.
        do i = 1, self%struct%natoms
            self%struct%atoms(i)%mass = skf_get_mass(self%struct%atoms(i)%symbol)
        end do

        call build_basis_system(self%struct, self%basis, self%state%bas)

        if (allocated(self%state%grad)) deallocate(self%state%grad)
        allocate(self%state%grad(3, self%struct%natoms))
        self%state%grad = 0.0_wp

        self%built = .true.
    end subroutine dftb_build


    subroutine dftb_execute(self)
        class(dftb_calc_t), intent(inout) :: self
        if (.not. self%built) call fatal("dftb", "execute before build")

        call solve_scc(self%struct, self%prop%scheme, self%prop%gamma_kind, &
                       self%prop%maxscc, self%prop%tolscc,                  &
                       self%prop%write_matrix, self%prop%mixing, self%state)

        self%state%e_rep = repulsive_energy(self%struct)
        select case (self%prop%scheme)
        case (SCHEME_BASIC)
            self%state%e_coul = 0.0_wp
        case (SCHEME_NOSCC, SCHEME_SCC)
            self%state%e_coul = coulomb_energy(self%state%gamma, self%state%dq)
        end select
        self%state%e_total = total_energy(self%state%e_band, &
                                          self%state%e_coul, self%state%e_rep)

        call zero_gradient(self%state%grad)
    end subroutine dftb_execute


    function dftb_get_total_energy(self) result(e)
        class(dftb_calc_t), intent(in) :: self
        real(wp) :: e
        e = self%state%e_total
    end function dftb_get_total_energy


    function dftb_get_total_gradient(self) result(g)
        class(dftb_calc_t), intent(in) :: self
        real(wp), allocatable :: g(:,:)
        allocate(g(3, self%struct%natoms))
        g = self%state%grad
    end function dftb_get_total_gradient


    subroutine dftb_write_output(self)
        class(dftb_calc_t), intent(inout) :: self
        real(wp), allocatable :: q(:)
        q = self%get_charges()
        call write_dftb_final(self%state%e_total, self%state%e_band, &
                              self%state%e_coul,  self%state%e_rep, q)
    end subroutine dftb_write_output


    !-- méthodes spécifiques DFTB --------------------------------------

    function get_band_energy(self) result(e)
        class(dftb_calc_t), intent(in) :: self
        real(wp) :: e
        e = self%state%e_band
    end function get_band_energy

    function get_coulomb_energy(self) result(e)
        class(dftb_calc_t), intent(in) :: self
        real(wp) :: e
        e = self%state%e_coul
    end function get_coulomb_energy

    function get_repulsive_energy(self) result(e)
        class(dftb_calc_t), intent(in) :: self
        real(wp) :: e
        e = self%state%e_rep
    end function get_repulsive_energy

    function get_charges(self) result(q)
        class(dftb_calc_t), intent(in) :: self
        real(wp), allocatable :: q(:)
        allocate(q(self%struct%natoms))
        if (allocated(self%state%q)) then
            q = self%state%q
        else
            q = 0.0_wp
        end if
    end function get_charges
end module dftb
