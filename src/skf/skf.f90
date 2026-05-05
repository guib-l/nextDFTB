!> Façade publique du calculateur SKF.
!>
!> Expose EXACTEMENT l'interface suivante :
!>   - subroutine init(struct, basis)
!>   - subroutine readslako()
!>   - subroutine build_repulsion()
!>   - subroutine build_electronic()
!>   - function get_repulsive(atom_A, atom_B, r)             result(e)
!>   - function get_overlaps  (atom_A, atom_B, r [,binding]) result(s)
!>   - function get_hamiltonian(atom_A, atom_B, r [,binding]) result(h)
!>   - function get_hubbard(atom_A) result(u)
!>   - function get_eps    (atom_A) result(eps)
!>
!> Aucune autre méthode publique. L'état est interne au module.
module skf
    use kinds,         only: wp
    use constants,     only: SYMBOL_LEN
    use structure_mod, only: structure_t
    use property,      only: property_basis_t
    use method_basis,  only: method_basis_t
    use slakos,        only: skf_store_t
    use readskf,       only: load_skf_store
    use electronic,    only: hs_pair_integrals, dhs_pair_integrals, binding_idx
    use repulsive_skf, only: pair_repulsive, pair_repulsive_spline, &
                              pair_repulsive_poly, pair_drepulsive
    use write_skf,     only: write_skf_summary
    use errors,        only: fatal
    implicit none
    private

    type(skf_store_t),         save :: store
    character(len=SYMBOL_LEN), allocatable, save :: sym_table(:)
    character(len=:),          allocatable, save :: src, ext, sep
    logical, save :: is_init    = .false.
    logical, save :: is_loaded  = .false.

    public :: init, readslako
    public :: build_repulsion, build_electronic
    public :: get_repulsive, get_overlaps, get_hamiltonian
    public :: get_repulsion_spline, get_repulsion_poly
    public :: get_dhamiltonian, get_doverlaps, get_drepulsive
    public :: get_hubbard, get_eps, get_occupations
    public :: get_mass
    public :: nelements, element_symbol, element_index
    public :: write_summary

    !> Calculateur SKF — extension de l'interface abstraite method_basis.
    !> Délègue intégralement aux procédures de module pour préserver
    !> l'API existante (utilisée par matel et dftb).
    type, extends(method_basis_t), public :: skf_calc_t
        logical :: initialized = .false.
        logical :: built       = .false.
    contains
        procedure :: init         => skf_calc_init
        procedure :: build        => skf_calc_build
        procedure :: execute      => skf_calc_execute
        procedure :: write_output => skf_calc_write_output
    end type skf_calc_t

contains

    subroutine init(struct, basis)
        type(structure_t),      intent(in) :: struct
        type(property_basis_t), intent(in) :: basis
        integer :: i

        call unique_symbols(struct, sym_table)

        src = trim(basis%src)
        ext = trim(basis%ext)
        sep = trim(basis%sep)

        is_init   = .true.
        is_loaded = .false.

        ! Décharge éventuelle avant rechargement.
        if (allocated(store%pair)) deallocate(store%pair)
        store%nelem = 0
        ! (le contenu détaillé est rempli par readslako)
        i = 0
    end subroutine init


    subroutine readslako()
        if (.not. is_init) call fatal("skf", "readslako called before init")
        call load_skf_store(sym_table, src, ext, sep, store)
        is_loaded = .true.
    end subroutine readslako


    subroutine build_repulsion()
        ! La répulsion est tabulée à la lecture (poly + spline). Aucun
        ! pré-calcul supplémentaire n'est nécessaire pour l'instant.
        if (.not. is_loaded) call fatal("skf", "build_repulsion before readslako")
    end subroutine build_repulsion


    subroutine build_electronic()
        ! La table H/S est interpolée à la volée par get_overlaps /
        ! get_hamiltonian. Hook réservé pour de futures pré-tabulations.
        if (.not. is_loaded) call fatal("skf", "build_electronic before readslako")
    end subroutine build_electronic


    function get_repulsive(atom_A, atom_B, r) result(e)
        character(len=*), intent(in) :: atom_A, atom_B
        real(wp),         intent(in) :: r
        real(wp) :: e
        e = get_repulsion_spline(atom_A, atom_B, r)
    end function get_repulsive


    function get_repulsion_spline(atom_A, atom_B, r) result(e)
        character(len=*), intent(in) :: atom_A, atom_B
        real(wp),         intent(in) :: r
        real(wp) :: e
        integer  :: ia, ib
        ia = element_index(atom_A)
        ib = element_index(atom_B)
        e = pair_repulsive_spline(store, ia, ib, r)
    end function get_repulsion_spline


    function get_repulsion_poly(atom_A, atom_B, r) result(e)
        character(len=*), intent(in) :: atom_A, atom_B
        real(wp),         intent(in) :: r
        real(wp) :: e
        integer  :: ia, ib
        ia = element_index(atom_A)
        ib = element_index(atom_B)
        e = pair_repulsive_poly(store, ia, ib, r)
    end function get_repulsion_poly


    function get_drepulsive(atom_A, atom_B, r) result(de)
        character(len=*), intent(in) :: atom_A, atom_B
        real(wp),         intent(in) :: r
        real(wp) :: de
        integer  :: ia, ib
        ia = element_index(atom_A)
        ib = element_index(atom_B)
        de = pair_drepulsive(store, ia, ib, r)
    end function get_drepulsive


    function get_overlaps(atom_A, atom_B, r, binding) result(s)
        character(len=*), intent(in)           :: atom_A, atom_B
        real(wp),         intent(in)           :: r
        character(len=*), intent(in), optional :: binding
        real(wp), allocatable :: s(:)
        real(wp) :: h_(10), s_(10)
        integer  :: ia, ib, k
        ia = element_index(atom_A)
        ib = element_index(atom_B)
        call hs_pair_integrals(store, ia, ib, r, h_, s_)
        if (present(binding)) then
            k = binding_idx(binding)
            allocate(s(1)); s(1) = s_(k)
        else
            allocate(s(10)); s(:) = s_
        end if
    end function get_overlaps


    function get_hamiltonian(atom_A, atom_B, r, binding) result(h)
        character(len=*), intent(in)           :: atom_A, atom_B
        real(wp),         intent(in)           :: r
        character(len=*), intent(in), optional :: binding
        real(wp), allocatable :: h(:)
        real(wp) :: h_(10), s_(10)
        integer  :: ia, ib, k
        ia = element_index(atom_A)
        ib = element_index(atom_B)
        call hs_pair_integrals(store, ia, ib, r, h_, s_)
        if (present(binding)) then
            k = binding_idx(binding)
            allocate(h(1)); h(1) = h_(k)
        else
            allocate(h(10)); h(:) = h_
        end if
    end function get_hamiltonian


    function get_doverlaps(atom_A, atom_B, r, binding) result(ds)
        character(len=*), intent(in)           :: atom_A, atom_B
        real(wp),         intent(in)           :: r
        character(len=*), intent(in), optional :: binding
        real(wp), allocatable :: ds(:)
        real(wp) :: dh_(10), ds_(10)
        integer  :: ia, ib, k
        ia = element_index(atom_A)
        ib = element_index(atom_B)
        call dhs_pair_integrals(store, ia, ib, r, dh_, ds_)
        if (present(binding)) then
            k = binding_idx(binding)
            allocate(ds(1)); ds(1) = ds_(k)
        else
            allocate(ds(10)); ds(:) = ds_
        end if
    end function get_doverlaps


    function get_dhamiltonian(atom_A, atom_B, r, binding) result(dh)
        character(len=*), intent(in)           :: atom_A, atom_B
        real(wp),         intent(in)           :: r
        character(len=*), intent(in), optional :: binding
        real(wp), allocatable :: dh(:)
        real(wp) :: dh_(10), ds_(10)
        integer  :: ia, ib, k
        ia = element_index(atom_A)
        ib = element_index(atom_B)
        call dhs_pair_integrals(store, ia, ib, r, dh_, ds_)
        if (present(binding)) then
            k = binding_idx(binding)
            allocate(dh(1)); dh(1) = dh_(k)
        else
            allocate(dh(10)); dh(:) = dh_
        end if
    end function get_dhamiltonian


    !> Hubbard du fichier homonucléaire (d, p, s).
    function get_hubbard(atom_A) result(u)
        character(len=*), intent(in) :: atom_A
        real(wp) :: u(3)
        integer  :: ia
        ia = element_index(atom_A)
        u = store%pair(ia, ia)%hubbard
    end function get_hubbard


    !> Énergies onsite du fichier homonucléaire (d, p, s).
    function get_eps(atom_A) result(eps)
        character(len=*), intent(in) :: atom_A
        real(wp) :: eps(3)
        integer  :: ia
        ia = element_index(atom_A)
        eps = store%pair(ia, ia)%e_onsite
    end function get_eps


    !> Occupations de valence neutre du fichier homonucléaire (d, p, s).
    function get_occupations(atom_A) result(occ)
        character(len=*), intent(in) :: atom_A
        real(wp) :: occ(3)
        integer  :: ia
        ia = element_index(atom_A)
        occ = store%pair(ia, ia)%occ
    end function get_occupations


    !> Masse atomique de l'élément (unités du fichier SKF homonucléaire, u.m.a.).
    function get_mass(atom_A) result(m)
        character(len=*), intent(in) :: atom_A
        real(wp) :: m
        integer  :: ia
        ia = element_index(atom_A)
        m = store%pair(ia, ia)%mass
    end function get_mass


    !> Résumé du magasin SKF (éléments chargés, paires).
    subroutine write_summary()
        if (.not. is_loaded) return
        call write_skf_summary(sym_table, store)
    end subroutine write_summary


    !-- Type-bound procedures de skf_calc_t ----------------------------

    subroutine skf_calc_init(self, struct, basis)
        class(skf_calc_t),      intent(inout) :: self
        type(structure_t),      intent(in)    :: struct
        type(property_basis_t), intent(in)    :: basis
        call init(struct, basis)
        self%initialized = .true.
        self%built       = .false.
    end subroutine skf_calc_init

    subroutine skf_calc_build(self)
        class(skf_calc_t), intent(inout) :: self
        if (.not. self%initialized) call fatal("skf", "build before init")
        call readslako()
        call build_repulsion()
        call build_electronic()
        self%built = .true.
    end subroutine skf_calc_build

    subroutine skf_calc_execute(self)
        class(skf_calc_t), intent(inout) :: self
        if (.not. self%built) call fatal("skf", "execute before build")
        ! Les intégrales sont interpolées à la demande ; rien à exécuter
        ! ici. Conservé pour la conformité à l'interface method_basis.
    end subroutine skf_calc_execute

    subroutine skf_calc_write_output(self)
        class(skf_calc_t), intent(inout) :: self
        call write_summary()
    end subroutine skf_calc_write_output


    !-- helpers internes (visibles aux autres modules pour itérer la table) -

    function nelements() result(n)
        integer :: n
        n = 0
        if (allocated(sym_table)) n = size(sym_table)
    end function nelements

    function element_symbol(i) result(sym)
        integer, intent(in) :: i
        character(len=SYMBOL_LEN) :: sym
        sym = sym_table(i)
    end function element_symbol

    function element_index(sym) result(idx)
        character(len=*), intent(in) :: sym
        integer :: idx, i
        idx = 0
        do i = 1, size(sym_table)
            if (trim(sym_table(i)) == trim(sym)) then
                idx = i; return
            end if
        end do
        call fatal("skf", "unknown element: "//trim(sym))
    end function element_index


    subroutine unique_symbols(struct, uniq)
        type(structure_t), intent(in) :: struct
        character(len=SYMBOL_LEN), allocatable, intent(out) :: uniq(:)
        character(len=SYMBOL_LEN) :: tmp(struct%natoms)
        integer :: i, j, n
        logical :: seen
        n = 0
        do i = 1, struct%natoms
            seen = .false.
            do j = 1, n
                if (trim(tmp(j)) == trim(struct%atoms(i)%symbol)) then
                    seen = .true.; exit
                end if
            end do
            if (.not. seen) then
                n = n + 1
                tmp(n) = struct%atoms(i)%symbol
            end if
        end do
        if (allocated(uniq)) deallocate(uniq)
        allocate(uniq(n))
        uniq(1:n) = tmp(1:n)
    end subroutine unique_symbols
end module skf
