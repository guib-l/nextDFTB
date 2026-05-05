!> Énergie répulsive paramétrée Slater-Koster (côté SKF).
!>
!> Deux variantes exposées :
!>   - pair_repulsive_spline : spline si disponible, polynôme sinon
!>   - pair_repulsive_poly   : polynôme c_poly seul (zone proche)
!> `pair_repulsive` est conservé comme alias par défaut sur la spline.
module repulsive_skf
    use kinds,  only: wp
    use slakos, only: skf_t, skf_store_t
    use interp, only: vrep_at_r, dvrep_at_r
    implicit none
    private

    public :: pair_repulsive
    public :: pair_repulsive_spline, pair_repulsive_poly
    public :: pair_drepulsive

contains

    function pair_repulsive(store, ea, eb, r) result(e)
        type(skf_store_t), intent(in) :: store
        integer,           intent(in) :: ea, eb
        real(wp),          intent(in) :: r
        real(wp) :: e
        e = pair_repulsive_spline(store, ea, eb, r)
    end function pair_repulsive


    !> Variante spline : utilise les segments si présents, sinon retombe
    !> sur le polynôme via vrep_at_r (comportement existant).
    function pair_repulsive_spline(store, ea, eb, r) result(e)
        type(skf_store_t), intent(in) :: store
        integer,           intent(in) :: ea, eb
        real(wp),          intent(in) :: r
        real(wp) :: e
        e = vrep_at_r(store%pair(ea, eb), r)
    end function pair_repulsive_spline


    !> Variante polynôme seul (3e ligne du SKF) : ignore une éventuelle
    !> spline et évalue c_poly sur (r - rcut) pour r < rcut.
    function pair_repulsive_poly(store, ea, eb, r) result(e)
        type(skf_store_t), intent(in) :: store
        integer,           intent(in) :: ea, eb
        real(wp),          intent(in) :: r
        real(wp) :: e
        type(skf_t) :: tmp
        e = 0.0_wp
        if (r >= store%pair(ea, eb)%rcut) return
        tmp = store%pair(ea, eb)
        tmp%has_spline = .false.
        if (allocated(tmp%segs)) deallocate(tmp%segs)
        e = vrep_at_r(tmp, r)
    end function pair_repulsive_poly


    !> Dérivée première de la répulsion (variante spline par défaut).
    function pair_drepulsive(store, ea, eb, r) result(de)
        type(skf_store_t), intent(in) :: store
        integer,           intent(in) :: ea, eb
        real(wp),          intent(in) :: r
        real(wp) :: de
        de = dvrep_at_r(store%pair(ea, eb), r)
    end function pair_drepulsive
end module repulsive_skf
