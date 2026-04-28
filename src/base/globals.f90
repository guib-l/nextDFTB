!> Types globaux partagés entre modules (géométrie, base, calcul, sortie).
module globals
    use kinds, only: wp
    implicit none
    private

    integer, parameter, public :: SYMBOL_LEN = 4
    integer, parameter, public :: PATH_LEN   = 256

    type, public :: geometry_t
        integer :: natoms = 0
        character(len=SYMBOL_LEN), allocatable :: symbols(:)   ! (natoms)
        real(wp),                  allocatable :: coords(:,:)  ! (3, natoms), bohr
        integer,                   allocatable :: groups(:)    ! (natoms)
    end type geometry_t

    type, public :: basis_t
        character(len=PATH_LEN) :: src  = ""
        character(len=16)       :: ext  = ".skf"
        character(len=4)        :: sep  = "-"
        character(len=16)       :: type = "VALENCE"
    end type basis_t

    type, public :: calc_t
        character(len=8) :: kind   = "DFTB"   ! DFTB | DFT | SKF
        logical          :: scc    = .false.
        integer          :: maxscc = 100
        real(wp)         :: tolscc = 1.0e-5_wp
    end type calc_t

    type, public :: output_t
        character(len=PATH_LEN) :: out = "results.out"
        character(len=PATH_LEN) :: log = "logs"
        logical                 :: log_on = .false.
    end type output_t

    type, public :: input_t
        type(geometry_t) :: geom
        type(basis_t)    :: basis
        type(calc_t)     :: calc
        type(output_t)   :: out
        logical          :: has_driver = .false.
    end type input_t

    !> Description d'un élément chimique tel que reconstruit depuis son SKF
    !> homonucléaire. Limité aux orbitales de valence (s, p, d).
    type, public :: element_t
        character(len=SYMBOL_LEN) :: symbol = ""
        integer  :: l_max  = 0          ! 0=s, 1=sp, 2=spd
        integer  :: n_orb  = 0          ! 1, 4, ou 9
        real(wp) :: e_s    = 0.0_wp
        real(wp) :: e_p    = 0.0_wp
        real(wp) :: e_d    = 0.0_wp
        real(wp) :: U_s    = 0.0_wp     ! Hubbard "atomique" (m. utilisé)
        real(wp) :: occ_s  = 0.0_wp
        real(wp) :: occ_p  = 0.0_wp
        real(wp) :: occ_d  = 0.0_wp
        real(wp) :: q_neutral = 0.0_wp  ! charge de référence (somme occ)
    end type element_t

    !> Carte (atome → orbitales globales) et table des éléments.
    type, public :: basis_system_t
        integer :: nelem = 0
        type(element_t), allocatable :: elems(:)        ! (nelem)
        integer, allocatable :: atom_elem(:)            ! (natoms) → idx in elems
        integer, allocatable :: atom_orb_start(:)       ! (natoms) 1-based
        integer, allocatable :: atom_norb(:)            ! (natoms)
        integer :: norb_total = 0
        integer :: nelec = 0
    end type basis_system_t
end module globals
