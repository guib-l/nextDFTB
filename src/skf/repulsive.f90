!> Énergie répulsive paramétrée Slater-Koster (côté SKF).
module repulsive_skf
    use kinds,  only: wp
    use slakos, only: skf_store_t
    use interp, only: vrep_at_r
    implicit none
    private

    public :: pair_repulsive

contains

    function pair_repulsive(store, ea, eb, r) result(e)
        type(skf_store_t), intent(in) :: store
        integer,           intent(in) :: ea, eb
        real(wp),          intent(in) :: r
        real(wp) :: e
        e = vrep_at_r(store%pair(ea, eb), r)
    end function pair_repulsive
end module repulsive_skf
