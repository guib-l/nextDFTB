!> Module electronic : interpolation des intégrales Slater-Koster
!> et accès aux propriétés électroniques par élément.
!>
!> Le module définit l'index `binding_idx` qui transforme un libellé de
!> liaison ("sss", "ppσ", ...) en indice (1..10) de la table SKF.
module electronic
    use kinds,  only: wp
    use slakos, only: skf_t, skf_store_t
    use interp, only: hs_at_r, dhs_at_r
    use errors, only: fatal
    implicit none
    private

    public :: hs_pair_integrals, dhs_pair_integrals
    public :: binding_idx

contains

    !> Renvoie les 10 intégrales SK (h, s) pour la paire (ea, eb) à r.
    subroutine hs_pair_integrals(store, ea, eb, r, h, s)
        type(skf_store_t), intent(in)  :: store
        integer,           intent(in)  :: ea, eb
        real(wp),          intent(in)  :: r
        real(wp),          intent(out) :: h(10), s(10)
        call hs_at_r(store%pair(ea, eb), r, h, s)
    end subroutine hs_pair_integrals


    !> Dérivée première par rapport à r des 10 intégrales SK (h, s).
    subroutine dhs_pair_integrals(store, ea, eb, r, dh, ds)
        type(skf_store_t), intent(in)  :: store
        integer,           intent(in)  :: ea, eb
        real(wp),          intent(in)  :: r
        real(wp),          intent(out) :: dh(10), ds(10)
        call dhs_at_r(store%pair(ea, eb), r, dh, ds)
    end subroutine dhs_pair_integrals


    !> Index 1..10 de la table SKF correspondant au libellé `binding`.
    !> Convention :
    !>   1=ddσ 2=ddπ 3=ddδ 4=pdσ 5=pdπ 6=ppσ 7=ppπ 8=sdσ 9=spσ 10=ssσ
    function binding_idx(binding) result(k)
        character(len=*), intent(in) :: binding
        integer :: k
        character(len=8) :: u
        integer :: i, c
        u = " "
        do i = 1, min(len(binding), len(u))
            c = iachar(binding(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) c = c + 32
            u(i:i) = achar(c)
        end do
        select case (trim(u))
        case ("sss", "sssigma", "ssσ");          k = 10
        case ("sps", "spsigma", "spσ");          k = 9
        case ("sds", "sdsigma", "sdσ");          k = 8
        case ("pps", "ppsigma", "ppσ");          k = 6
        case ("ppp", "pppi",    "ppπ");          k = 7
        case ("pds", "pdsigma", "pdσ");          k = 4
        case ("pdp", "pdpi",    "pdπ");          k = 5
        case ("dds", "ddsigma", "ddσ");          k = 1
        case ("ddp", "ddpi",    "ddπ");          k = 2
        case ("ddd", "dddelta", "ddδ");          k = 3
        case default
            call fatal("electronic", "unknown binding: "//trim(binding))
            k = 0
        end select
    end function binding_idx
end module electronic
