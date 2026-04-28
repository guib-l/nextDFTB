!> Lecteur de fichiers Slater-Koster (.skf), format standard DFTB+.
!>
!> Format supporté :
!>   - homonucléaire et hétéronucléaire (PAS le format étendu `@`)
!>   - grille uniforme (dr, n_grid)
!>   - section optionnelle "Spline" pour la répulsion
!>
!> Référence : manuel DFTB+ (file format SKF).
module readskf
    use kinds,  only: wp
    use errors, only: fatal, warn
    implicit none
    private

    integer, parameter, public :: SKF_NCOL = 20

    !> Une intervalle cubique du spline : r in [r1, r2], y(r) = sum c_k (r-r1)^k
    type, public :: spline_seg_t
        real(wp) :: r1, r2
        real(wp) :: c(0:5) = 0.0_wp   ! ordre <=5 (dernier segment)
        integer  :: order = 3
    end type spline_seg_t

    type, public :: skf_t
        logical  :: homonuclear = .false.
        real(wp) :: dr     = 0.0_wp
        integer  :: ngrid  = 0
        ! homonuclear only: onsite energies (d, p, s), Hubbard (d, p, s), occupations (d, p, s)
        real(wp) :: e_onsite(3) = 0.0_wp     ! ordre (d, p, s)
        real(wp) :: spe         = 0.0_wp
        real(wp) :: hubbard(3)  = 0.0_wp
        real(wp) :: occ(3)      = 0.0_wp
        ! repulsive polynomial (zone proche)
        real(wp) :: mass = 0.0_wp
        real(wp) :: c_poly(2:9) = 0.0_wp
        real(wp) :: rcut = 0.0_wp
        real(wp) :: d_poly(10) = 0.0_wp
        ! H/S table : (SKF_NCOL, ngrid)
        real(wp), allocatable :: hs(:,:)
        ! Spline répulsif
        logical :: has_spline = .false.
        real(wp) :: spline_a1 = 0.0_wp
        real(wp) :: spline_a2 = 0.0_wp
        real(wp) :: spline_a3 = 0.0_wp
        real(wp) :: spline_cutoff = 0.0_wp
        type(spline_seg_t), allocatable :: segs(:)
    end type skf_t

    !> Conteneur de toutes les paires SKF nécessaires au calcul.
    !> Indice = (i_elem, j_elem). Dimension nelem × nelem.
    type, public :: skf_store_t
        integer :: nelem = 0
        type(skf_t), allocatable :: pair(:,:)
    end type skf_store_t

    public :: load_skf, load_skf_store

contains

    subroutine load_skf_store(symbols, src, ext, sep, store)
        character(len=*), intent(in)  :: symbols(:)   ! liste unique d'éléments
        character(len=*), intent(in)  :: src, ext, sep
        type(skf_store_t), intent(out) :: store
        integer :: n, i, j
        character(len=:), allocatable :: path

        n = size(symbols)
        store%nelem = n
        allocate(store%pair(n, n))

        do i = 1, n
            do j = 1, n
                path = trim(src) // trim(symbols(i)) // trim(sep) // &
                       trim(symbols(j)) // trim(ext)
                call load_skf(path, i == j, store%pair(i, j))
            end do
        end do
    end subroutine load_skf_store


    subroutine load_skf(filename, homo, skf)
        character(len=*), intent(in)  :: filename
        logical,          intent(in)  :: homo
        type(skf_t),      intent(out) :: skf

        integer :: u, ios, i
        character(len=512) :: line
        real(wp) :: row20(SKF_NCOL)
        real(wp) :: r2(20)

        skf%homonuclear = homo

        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) call fatal("readskf", "cannot open: "//trim(filename))

        ! ligne 1 : dr, n_grid (peut être préfixée par "@" pour format étendu - non supporté)
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) call fatal("readskf", "empty file: "//trim(filename))
        if (line(1:1) == '@') call fatal("readskf", "extended SKF format not supported")
        read(line, *, iostat=ios) skf%dr, skf%ngrid
        if (ios /= 0) call fatal("readskf", "bad header in "//trim(filename))

        if (homo) then
            ! ligne 2 : Ed Ep Es SPE Ud Up Us fd fp fs
            read(u, *, iostat=ios) skf%e_onsite(1), skf%e_onsite(2), skf%e_onsite(3), &
                                   skf%spe, &
                                   skf%hubbard(1), skf%hubbard(2), skf%hubbard(3), &
                                   skf%occ(1),     skf%occ(2),     skf%occ(3)
            if (ios /= 0) call fatal("readskf", "bad onsite line")
        end if

        ! ligne suivante : 20 valeurs (mass c2..c9 rcut d1..d10)
        read(u, *, iostat=ios) r2
        if (ios /= 0) call fatal("readskf", "bad repulsive polynomial line")
        skf%mass        = r2(1)
        skf%c_poly(2:9) = r2(2:9)
        skf%rcut        = r2(10)
        skf%d_poly(1:10)= r2(11:20)

        allocate(skf%hs(SKF_NCOL, skf%ngrid))
        do i = 1, skf%ngrid
            read(u, *, iostat=ios) row20
            if (ios /= 0) call fatal("readskf", "bad HS row")
            skf%hs(:, i) = row20
        end do

        ! cherche éventuelle section "Spline"
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, "Spline") > 0) then
                call read_spline(u, skf)
                exit
            end if
        end do
        close(u)
    end subroutine load_skf


    subroutine read_spline(u, skf)
        integer,     intent(in)    :: u
        type(skf_t), intent(inout) :: skf
        integer :: nint_, i, ios
        real(wp) :: r1, r2, c0, c1, c2, c3, c4, c5

        read(u, *, iostat=ios) nint_, skf%spline_cutoff
        if (ios /= 0) call fatal("readskf", "bad spline header")
        read(u, *, iostat=ios) skf%spline_a1, skf%spline_a2, skf%spline_a3
        if (ios /= 0) call fatal("readskf", "bad spline exp")

        allocate(skf%segs(nint_))
        do i = 1, nint_ - 1
            read(u, *, iostat=ios) r1, r2, c0, c1, c2, c3
            if (ios /= 0) call fatal("readskf", "bad spline segment")
            skf%segs(i)%r1 = r1; skf%segs(i)%r2 = r2
            skf%segs(i)%c(0:3) = (/ c0, c1, c2, c3 /)
            skf%segs(i)%order  = 3
        end do
        ! dernier segment : c0..c5
        read(u, *, iostat=ios) r1, r2, c0, c1, c2, c3, c4, c5
        if (ios /= 0) call fatal("readskf", "bad final spline segment")
        skf%segs(nint_)%r1 = r1; skf%segs(nint_)%r2 = r2
        skf%segs(nint_)%c(0:5) = (/ c0, c1, c2, c3, c4, c5 /)
        skf%segs(nint_)%order  = 5

        skf%has_spline = .true.
    end subroutine read_spline
end module readskf
