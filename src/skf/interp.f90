!> Interpolation des paramètres SKF (table H/S et énergie répulsive).
module interp
    use kinds,  only: wp
    use slakos, only: skf_t
    implicit none
    private

    public :: hs_at_r, vrep_at_r

contains

    !> Renvoie h(1:10) et s(1:10) interpolés linéairement à r (en bohr).
    !> Convention SKF (DFTB+) :
    !>   1=ddσ 2=ddπ 3=ddδ 4=pdσ 5=pdπ 6=ppσ 7=ppπ 8=sdσ 9=spσ 10=ssσ
    !> Retourne 0 si r est hors-grille.
    subroutine hs_at_r(skf, r, h, s)
        type(skf_t), intent(in)  :: skf
        real(wp),    intent(in)  :: r
        real(wp),    intent(out) :: h(10), s(10)
        real(wp) :: t
        integer  :: i

        h = 0.0_wp; s = 0.0_wp
        if (r <= 0.0_wp .or. skf%ngrid < 2) return

        t = r / skf%dr
        i = int(t) + 1
        if (i < 1 .or. i >= skf%ngrid) return

        t = t - real(i - 1, wp)
        h(1:10) = (1.0_wp - t) * skf%hs(1:10,  i) + t * skf%hs(1:10,  i+1)
        s(1:10) = (1.0_wp - t) * skf%hs(11:20, i) + t * skf%hs(11:20, i+1)
    end subroutine hs_at_r


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
end module interp
