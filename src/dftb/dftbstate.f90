!> Objet `dftbstate` : état du dernier calcul DFTB.
!>
!> Réunit le système orbitalaire (mapping atome → orbitales) et les
!> matrices résultantes (H, S, C, P, gamma) ainsi que les énergies et
!> les charges Mulliken.
module dftbstate
    use kinds,         only: wp
    use constants,     only: SYMBOL_LEN
    implicit none
    private

    type, public :: element_info_t
        character(len=SYMBOL_LEN) :: symbol = ""
        integer  :: l_max  = 0          ! 0=s, 1=sp, 2=spd
        integer  :: n_orb  = 0          ! 1, 4, ou 9
        real(wp) :: e_s    = 0.0_wp
        real(wp) :: e_p    = 0.0_wp
        real(wp) :: e_d    = 0.0_wp
        real(wp) :: U_s    = 0.0_wp     ! Hubbard valence (m. utilisé)
        real(wp) :: occ_s  = 0.0_wp
        real(wp) :: occ_p  = 0.0_wp
        real(wp) :: occ_d  = 0.0_wp
        real(wp) :: q_neutral = 0.0_wp  ! occupation neutre cumulée
    end type element_info_t

    type, public :: basis_system_t
        integer :: nelem = 0
        type(element_info_t), allocatable :: elems(:)
        integer, allocatable :: atom_elem(:)        ! (natoms) idx in elems
        integer, allocatable :: atom_orb_start(:)   ! (natoms) 1-based
        integer, allocatable :: atom_norb(:)        ! (natoms)
        integer :: norb_total = 0
        integer :: nelec      = 0
    end type basis_system_t

    type, public :: dftbstate_t
        type(basis_system_t)  :: bas
        real(wp), allocatable :: H(:,:)
        real(wp), allocatable :: S(:,:)
        real(wp), allocatable :: C(:,:)
        real(wp), allocatable :: eig(:)
        real(wp), allocatable :: P(:,:)
        real(wp), allocatable :: occ(:)
        real(wp), allocatable :: q(:)
        real(wp), allocatable :: dq(:)
        real(wp), allocatable :: gamma(:,:)
        real(wp), allocatable :: grad(:,:)
        real(wp) :: e_total = 0.0_wp
        real(wp) :: e_band  = 0.0_wp
        real(wp) :: e_coul  = 0.0_wp
        real(wp) :: e_rep   = 0.0_wp
        integer  :: niter     = 0
        logical  :: converged = .false.
    end type dftbstate_t
end module dftbstate
