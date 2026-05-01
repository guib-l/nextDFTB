!> Transformations Slater-Koster : rotation des intégrales tabulées
!> (ssσ, spσ, ppσ, ppπ, ...) vers les blocs angulaires (s, px, py, pz)
!> dans le repère cartésien.
!>
!> La table d'intégrales suit la convention DFTB+ :
!>   1=ddσ 2=ddπ 3=ddδ 4=pdσ 5=pdπ 6=ppσ 7=ppπ 8=sdσ 9=spσ 10=ssσ
!>
!> Ordre des orbitales par atome : 1=s, 2=px, 3=py, 4=pz.
!>
!> Référence : Slater & Koster, Phys. Rev. 94, 1498 (1954), Table I.
module skrot
    use kinds, only: wp
    implicit none
    private

    public :: sk_block

contains

    !> Construit les blocs (na × nb) de H et S pour une paire d'atomes.
    !>
    !>   h_ab, s_ab : intégrales SK lues de A→B (10 valeurs chacune)
    !>   h_ba, s_ba : intégrales SK lues de B→A (utilisées pour les
    !>                composantes asymétriques s↔p)
    !>   dir(3)     : vecteur unitaire A → B
    !>   na, nb     : nombre d'orbitales (1 ou 4 ; le cas 9 (d) n'est
    !>                pas pris en charge dans cette version minimale)
    subroutine sk_block(h_ab, s_ab, h_ba, s_ba, dir, na, nb, hb, sb)
        real(wp), intent(in)  :: h_ab(10), s_ab(10)
        real(wp), intent(in)  :: h_ba(10), s_ba(10)
        real(wp), intent(in)  :: dir(3)
        integer,  intent(in)  :: na, nb
        real(wp), intent(out) :: hb(:,:), sb(:,:)

        real(wp) :: l_, m_, n_

        hb = 0.0_wp; sb = 0.0_wp
        l_ = dir(1); m_ = dir(2); n_ = dir(3)

        call rotate_one(h_ab, h_ba, l_, m_, n_, na, nb, hb)
        call rotate_one(s_ab, s_ba, l_, m_, n_, na, nb, sb)
    end subroutine sk_block


    subroutine rotate_one(v_ab, v_ba, l_, m_, n_, na, nb, blk)
        real(wp), intent(in)  :: v_ab(10), v_ba(10)
        real(wp), intent(in)  :: l_, m_, n_
        integer,  intent(in)  :: na, nb
        real(wp), intent(out) :: blk(:,:)
        real(wp) :: vss_s, vsp_s, vps_s, vpp_s, vpp_p

        vss_s = v_ab(10)
        vsp_s = v_ab(9)            ! s sur A, p sur B
        vps_s = -v_ba(9)           ! p sur A, s sur B (signe (-1)^(l_A+l_B) = -1)
        vpp_s = v_ab(6)
        vpp_p = v_ab(7)

        blk(1,1) = vss_s

        if (nb >= 4) then
            blk(1, 2) = l_ * vsp_s
            blk(1, 3) = m_ * vsp_s
            blk(1, 4) = n_ * vsp_s
        end if

        if (na >= 4) then
            blk(2, 1) = -l_ * vps_s
            blk(3, 1) = -m_ * vps_s
            blk(4, 1) = -n_ * vps_s
        end if

        if (na >= 4 .and. nb >= 4) then
            blk(2, 2) = l_*l_ * vpp_s + (1.0_wp - l_*l_) * vpp_p
            blk(2, 3) = l_*m_ * (vpp_s - vpp_p)
            blk(2, 4) = l_*n_ * (vpp_s - vpp_p)
            blk(3, 2) = blk(2, 3)
            blk(3, 3) = m_*m_ * vpp_s + (1.0_wp - m_*m_) * vpp_p
            blk(3, 4) = m_*n_ * (vpp_s - vpp_p)
            blk(4, 2) = blk(2, 4)
            blk(4, 3) = blk(3, 4)
            blk(4, 4) = n_*n_ * vpp_s + (1.0_wp - n_*n_) * vpp_p
        end if
    end subroutine rotate_one
end module skrot
