!> Lecteur de fichiers Slater-Koster (.skf), format standard DFTB+.
!>
!> Référence : `slakoformat.pdf` du dépôt. Format supporté :
!>   - homonucléaire et hétéronucléaire (PAS le format étendu `@`)
!>   - grille uniforme (dr, n_grid)
!>   - section optionnelle "Spline" pour la répulsion
module readskf
    use kinds,  only: wp
    use slakos, only: skf_t, skf_store_t, spline_seg_t, SKF_NCOL, SKF_NHS
    use errors, only: fatal
    implicit none
    private

    public :: load_skf, load_skf_store

contains

    subroutine load_skf_store(symbols, src, ext, sep, store)
        character(len=*), intent(in)   :: symbols(:)
        character(len=*), intent(in)   :: src, ext, sep
        type(skf_store_t), intent(out) :: store
        integer :: n, i, j
        character(len=:), allocatable :: path, root

        n = size(symbols)
        store%nelem = n
        allocate(store%pair(n, n))

        root = trim(src)
        if (len(root) > 0) then
            if (root(len(root):len(root)) /= "/") root = root // "/"
        end if

        do i = 1, n
            do j = 1, n
                path = root // trim(symbols(i)) // trim(sep) // &
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

        read(u, '(a)', iostat=ios) line
        if (ios /= 0) call fatal("readskf", "empty file: "//trim(filename))
        if (line(1:1) == '@') call fatal("readskf", "extended SKF format not supported")
        read(line, *, iostat=ios) skf%dr, skf%ngrid
        if (ios /= 0) call fatal("readskf", "bad header in "//trim(filename))

        if (homo) then
            read(u, *, iostat=ios) skf%e_onsite(1), skf%e_onsite(2), skf%e_onsite(3), &
                                   skf%spe, &
                                   skf%hubbard(1), skf%hubbard(2), skf%hubbard(3), &
                                   skf%occ(1),     skf%occ(2),     skf%occ(3)
            if (ios /= 0) call fatal("readskf", "bad onsite line in "//trim(filename))
        end if

        read(u, *, iostat=ios) r2
        if (ios /= 0) call fatal("readskf", "bad repulsive polynomial line")
        skf%mass         = r2(1)
        skf%c_poly(2:9)  = r2(2:9)
        skf%rcut         = r2(10)
        skf%d_poly(1:10) = r2(11:20)

        allocate(skf%h(SKF_NHS, skf%ngrid))
        allocate(skf%s(SKF_NHS, skf%ngrid))

        do i = 1, skf%ngrid
            read(u, *, iostat=ios) row20
            if (ios /= 0) call fatal("readskf", "bad HS row")
            skf%h(1:SKF_NHS, i) = row20(1:SKF_NHS)
            skf%s(1:SKF_NHS, i) = row20(SKF_NHS+1:SKF_NCOL)
        end do

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
        if (ios /= 0) call fatal("readskf", "bad spline exponential")

        allocate(skf%segs(nint_))
        do i = 1, nint_ - 1
            read(u, *, iostat=ios) r1, r2, c0, c1, c2, c3
            if (ios /= 0) call fatal("readskf", "bad spline segment")
            skf%segs(i)%r1 = r1; skf%segs(i)%r2 = r2
            skf%segs(i)%c(0:3) = [ c0, c1, c2, c3 ]
            skf%segs(i)%order  = 3
        end do
        read(u, *, iostat=ios) r1, r2, c0, c1, c2, c3, c4, c5
        if (ios /= 0) call fatal("readskf", "bad final spline segment")
        skf%segs(nint_)%r1 = r1; skf%segs(nint_)%r2 = r2
        skf%segs(nint_)%c(0:5) = [ c0, c1, c2, c3, c4, c5 ]
        skf%segs(nint_)%order  = 5

        skf%has_spline = .true.
    end subroutine read_spline
end module readskf
