!> Interpolation des paramètres SKF (table H/S et énergie répulsive).
!>
!> Utilise polint (Neville, ordre fixé) sur une fenêtre centrée autour
!> de r. Les dérivées analytiques sont fournies pour la répulsion ;
!> celles de H/S passent par différence finie sur polint.
module interp
    use kinds,      only: wp
    use slakos,     only: skf_t, SKF_NHS
    use polint_mod, only: polint, dpolint, n_order_polint
    implicit none
    private

    public :: hs_at_r, dhs_at_r
    public :: vrep_at_r, dvrep_at_r

contains

    !> Renvoie h(1:10) et s(1:10) interpolés par polint à r (en bohr).
    !> Convention SKF (DFTB+) :
    !>   1=ddσ 2=ddπ 3=ddδ 4=pdσ 5=pdπ 6=ppσ 7=ppπ 8=sdσ 9=spσ 10=ssσ
    !> Retourne 0 si r est hors-grille.
    subroutine hs_at_r(skf, r, h, s)
        type(skf_t), intent(in)  :: skf
        real(wp),    intent(in)  :: r
        real(wp),    intent(out) :: h(SKF_NHS), s(SKF_NHS)

        real(wp) :: xa(n_order_polint), ya(n_order_polint)
        real(wp) :: y, dy
        integer  :: i0, k, c

        h = 0.0_wp; s = 0.0_wp
        if (r <= 0.0_wp .or. skf%ngrid < n_order_polint) return
        if (r >= real(skf%ngrid - 1, wp) * skf%dr) return

        call window_indices(skf, r, i0)

        do k = 1, n_order_polint
            xa(k) = real(i0 + k - 2, wp) * skf%dr
        end do

        do c = 1, SKF_NHS
            do k = 1, n_order_polint
                ya(k) = skf%h(c, i0 + k - 1)
            end do
            call polint(xa, ya, r, y, dy)
            h(c) = y
            do k = 1, n_order_polint
                ya(k) = skf%s(c, i0 + k - 1)
            end do
            call polint(xa, ya, r, y, dy)
            s(c) = y
        end do
    end subroutine hs_at_r


    !> Dérivée première de h et s par rapport à r (différence finie via
    !> dpolint sur la même fenêtre).
    subroutine dhs_at_r(skf, r, dh, ds)
        type(skf_t), intent(in)  :: skf
        real(wp),    intent(in)  :: r
        real(wp),    intent(out) :: dh(SKF_NHS), ds(SKF_NHS)

        real(wp) :: xa(n_order_polint), ya(n_order_polint)
        integer  :: i0, k, c

        dh = 0.0_wp; ds = 0.0_wp
        if (r <= 0.0_wp .or. skf%ngrid < n_order_polint) return
        if (r >= real(skf%ngrid - 1, wp) * skf%dr) return

        call window_indices(skf, r, i0)

        do k = 1, n_order_polint
            xa(k) = real(i0 + k - 2, wp) * skf%dr
        end do

        do c = 1, SKF_NHS
            do k = 1, n_order_polint
                ya(k) = skf%h(c, i0 + k - 1)
            end do
            call dpolint(xa, ya, r, dh(c))
            do k = 1, n_order_polint
                ya(k) = skf%s(c, i0 + k - 1)
            end do
            call dpolint(xa, ya, r, ds(c))
        end do
    end subroutine dhs_at_r


    !> Énergie répulsive paramétrée par le SKF.
    !> Si spline disponible, utilise (exp + segments). Sinon : polynôme
    !> c2(r-rcut)^2 + ... + c9(r-rcut)^9 pour r < rcut.
    function vrep_at_r(skf, r) result(v)
        type(skf_t), intent(in) :: skf
        real(wp),    intent(in) :: r
        real(wp) :: v
        integer  :: i
        real(wp) :: dr_

        v = 0.0_wp
        if (skf%has_spline) then
            if (r >= skf%spline_cutoff) return
            if (allocated(skf%segs)) then
                if (r < skf%segs(1)%r1) then
                    v = exp(-skf%spline_a1 * r + skf%spline_a2) + skf%spline_a3
                    return
                end if
                do i = 1, size(skf%segs)
                    if (r >= skf%segs(i)%r1 .and. r < skf%segs(i)%r2) then
                        dr_ = r - skf%segs(i)%r1
                        v = poly_eval(skf%segs(i)%c, skf%segs(i)%order, dr_)
                        return
                    end if
                end do
            end if
        else
            if (r >= skf%rcut) return
            dr_ = r - skf%rcut
            v =   skf%c_poly(2) * dr_**2 &
                + skf%c_poly(3) * dr_**3 &
                + skf%c_poly(4) * dr_**4 &
                + skf%c_poly(5) * dr_**5 &
                + skf%c_poly(6) * dr_**6 &
                + skf%c_poly(7) * dr_**7 &
                + skf%c_poly(8) * dr_**8 &
                + skf%c_poly(9) * dr_**9
        end if
    end function vrep_at_r


    !> Dérivée première analytique de vrep par rapport à r.
    function dvrep_at_r(skf, r) result(dv)
        type(skf_t), intent(in) :: skf
        real(wp),    intent(in) :: r
        real(wp) :: dv
        integer  :: i, k
        real(wp) :: dr_

        dv = 0.0_wp
        if (skf%has_spline) then
            if (r >= skf%spline_cutoff) return
            if (allocated(skf%segs)) then
                if (r < skf%segs(1)%r1) then
                    dv = -skf%spline_a1 * exp(-skf%spline_a1 * r + skf%spline_a2)
                    return
                end if
                do i = 1, size(skf%segs)
                    if (r >= skf%segs(i)%r1 .and. r < skf%segs(i)%r2) then
                        dr_ = r - skf%segs(i)%r1
                        dv = poly_deriv(skf%segs(i)%c, skf%segs(i)%order, dr_)
                        return
                    end if
                end do
            end if
        else
            if (r >= skf%rcut) return
            dr_ = r - skf%rcut
            do k = 2, 9
                dv = dv + real(k, wp) * skf%c_poly(k) * dr_**(k - 1)
            end do
        end if
    end function dvrep_at_r


    !> Sélectionne i0 tel que la fenêtre [i0 .. i0+n-1] (1-based) reste
    !> dans la grille et soit centrée au mieux autour de r.
    pure subroutine window_indices(skf, r, i0)
        type(skf_t), intent(in)  :: skf
        real(wp),    intent(in)  :: r
        integer,     intent(out) :: i0
        integer :: ic
        ic = int(r / skf%dr) + 1
        i0 = ic - n_order_polint / 2 + 1
        if (i0 < 1) i0 = 1
        if (i0 > skf%ngrid - n_order_polint + 1) &
            i0 = skf%ngrid - n_order_polint + 1
    end subroutine window_indices


    pure function poly_eval(c, ord, x) result(y)
        real(wp), intent(in) :: c(0:5), x
        integer,  intent(in) :: ord
        real(wp) :: y
        integer  :: k
        y = 0.0_wp
        do k = 0, ord
            y = y + c(k) * x**k
        end do
    end function poly_eval


    pure function poly_deriv(c, ord, x) result(y)
        real(wp), intent(in) :: c(0:5), x
        integer,  intent(in) :: ord
        real(wp) :: y
        integer  :: k
        y = 0.0_wp
        do k = 1, ord
            y = y + real(k, wp) * c(k) * x**(k - 1)
        end do
    end function poly_deriv
end module interp
