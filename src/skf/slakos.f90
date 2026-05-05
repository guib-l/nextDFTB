!> Objet `slakos` : conteneur des paramètres Slater-Koster.
!>
!> Définit les types persistants pour stocker une paire SKF (`skf_t`),
!> ses segments spline répulsifs (`spline_seg_t`), et le magasin global
!> de toutes les paires (`skf_store_t`).
module slakos
    use kinds, only: wp
    implicit none
    private

    integer, parameter, public :: SKF_NCOL = 20   ! taille d'une ligne du fichier
    integer, parameter, public :: SKF_NHS  = 10   ! nb d'intégrales SK distinctes

    type, public :: spline_seg_t
        real(wp) :: r1 = 0.0_wp
        real(wp) :: r2 = 0.0_wp
        real(wp) :: c(0:5) = 0.0_wp
        integer  :: order  = 3
    end type spline_seg_t

    type, public :: skf_t
        logical  :: homonuclear = .false.
        real(wp) :: dr     = 0.0_wp
        integer  :: ngrid  = 0
        ! homonuclear seulement : onsite (d, p, s), Hubbard (d, p, s), occ (d, p, s)
        real(wp) :: e_onsite(3) = 0.0_wp
        real(wp) :: spe         = 0.0_wp
        real(wp) :: hubbard(3)  = 0.0_wp
        real(wp) :: occ(3)      = 0.0_wp
        ! polynôme répulsif (zone proche)
        real(wp) :: mass        = 0.0_wp
        real(wp) :: c_poly(2:9) = 0.0_wp
        real(wp) :: rcut        = 0.0_wp
        real(wp) :: d_poly(10)  = 0.0_wp
        ! Tables H et S séparées : (SKF_NHS, ngrid)
        ! Convention DFTB+ :
        !   1=ddσ 2=ddπ 3=ddδ 4=pdσ 5=pdπ 6=ppσ 7=ppπ 8=sdσ 9=spσ 10=ssσ
        real(wp), allocatable :: h(:,:)
        real(wp), allocatable :: s(:,:)
        ! Spline répulsif
        logical  :: has_spline    = .false.
        real(wp) :: spline_a1     = 0.0_wp
        real(wp) :: spline_a2     = 0.0_wp
        real(wp) :: spline_a3     = 0.0_wp
        real(wp) :: spline_cutoff = 0.0_wp
        type(spline_seg_t), allocatable :: segs(:)
    end type skf_t

    type, public :: skf_store_t
        integer :: nelem = 0
        type(skf_t), allocatable :: pair(:,:)
    end type skf_store_t
end module slakos
