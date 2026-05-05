!> Transformations Slater-Koster : rotation des intégrales tabulées vers
!> les blocs angulaires cartésiens.
!>
!> Convention DFTB+ pour la table d'intégrales :
!>   1=ddσ 2=ddπ 3=ddδ 4=pdσ 5=pdπ 6=ppσ 7=ppπ 8=sdσ 9=spσ 10=ssσ
!>
!> Ordre des orbitales par atome :
!>   1=s, 2=px, 3=py, 4=pz,
!>   5=dxy, 6=dyz, 7=dzx, 8=dx2−y2, 9=d3z2−r2
!>
!> Les directions sont exprimées en cosinus directeurs (l, m, n) avec
!> α = l²+m², β = l²−m². Référence : Slater & Koster Phys. Rev. 94,
!> 1498 (1954), Table I.
module skrot
    use kinds, only: wp
    implicit none
    private

    public :: sk_block

    real(wp), parameter :: SQRT3 = 1.7320508075688772_wp

contains

    !> Construit les blocs (na × nb) de H et S pour une paire d'atomes.
    !>
    !>   h_ab, s_ab : intégrales SK lues de A→B (10 valeurs chacune)
    !>   h_ba, s_ba : intégrales SK lues de B→A (utilisées pour les
    !>                composantes asymétriques où l_A ≠ l_B)
    !>   dir(3)     : vecteur unitaire A → B
    !>   na, nb     : nombre d'orbitales (1, 4 ou 9)
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


    !> Remplit blk(1:na,1:nb). Les éléments sont issus de v_ab pour les
    !> sous-blocs symétriques (ss, pp, dd) et combinés avec v_ba pour
    !> les anti-diagonaux (sp, sd, pd) avec le signe (-1)^(l_A+l_B).
    subroutine rotate_one(v_ab, v_ba, l_, m_, n_, na, nb, blk)
        real(wp), intent(in)  :: v_ab(10), v_ba(10)
        real(wp), intent(in)  :: l_, m_, n_
        integer,  intent(in)  :: na, nb
        real(wp), intent(out) :: blk(:,:)

        ! ss
        blk(1, 1) = v_ab(10)

        ! s (A) ↔ p, d (B)
        if (nb >= 4) call sp_row(v_ab, l_, m_, n_, blk, 1)
        if (nb >= 9) call sd_row(v_ab, l_, m_, n_, blk, 1)

        ! p, d (A) ↔ s (B)  : signe (-1)^(l_A+l_B) sur v_ba
        if (na >= 4) call ps_col(v_ba, l_, m_, n_, blk, 1)
        if (na >= 9) call ds_col(v_ba, l_, m_, n_, blk, 1)

        ! pp
        if (na >= 4 .and. nb >= 4) call pp_block(v_ab, l_, m_, n_, blk)

        ! p (A) ↔ d (B)
        if (na >= 4 .and. nb >= 9) call pd_block(v_ab, l_, m_, n_, blk)

        ! d (A) ↔ p (B)  : signe (-1)^(l_A+l_B) = -1
        if (na >= 9 .and. nb >= 4) call dp_block(v_ba, l_, m_, n_, blk)

        ! dd
        if (na >= 9 .and. nb >= 9) call dd_block(v_ab, l_, m_, n_, blk)
    end subroutine rotate_one


    !> s (ligne mu_s) ↔ p (colonnes 2..4) sur l'atome B.
    pure subroutine sp_row(v, l_, m_, n_, blk, mu_s)
        real(wp), intent(in)    :: v(10), l_, m_, n_
        real(wp), intent(inout) :: blk(:,:)
        integer,  intent(in)    :: mu_s
        real(wp) :: vsps
        vsps = v(9)
        blk(mu_s, 2) = l_ * vsps
        blk(mu_s, 3) = m_ * vsps
        blk(mu_s, 4) = n_ * vsps
    end subroutine sp_row


    !> s (ligne mu_s) ↔ d (colonnes 5..9) sur l'atome B.
    pure subroutine sd_row(v, l_, m_, n_, blk, mu_s)
        real(wp), intent(in)    :: v(10), l_, m_, n_
        real(wp), intent(inout) :: blk(:,:)
        integer,  intent(in)    :: mu_s
        real(wp) :: vsds, alpha, beta, z2
        vsds  = v(8)
        alpha = l_*l_ + m_*m_
        beta  = l_*l_ - m_*m_
        z2    = n_*n_
        blk(mu_s, 5) =  SQRT3 * l_ * m_                    * vsds
        blk(mu_s, 6) =  SQRT3 * m_ * n_                    * vsds
        blk(mu_s, 7) =  SQRT3 * n_ * l_                    * vsds
        blk(mu_s, 8) =  0.5_wp * SQRT3 * beta              * vsds
        blk(mu_s, 9) = (z2 - 0.5_wp * alpha)               * vsds
    end subroutine sd_row


    !> p (lignes 2..4) ↔ s (colonne nu_s). v est v_ba ; signe (-1).
    pure subroutine ps_col(v_ba, l_, m_, n_, blk, nu_s)
        real(wp), intent(in)    :: v_ba(10), l_, m_, n_
        real(wp), intent(inout) :: blk(:,:)
        integer,  intent(in)    :: nu_s
        real(wp) :: vsps
        vsps = v_ba(9)
        blk(2, nu_s) = -l_ * vsps
        blk(3, nu_s) = -m_ * vsps
        blk(4, nu_s) = -n_ * vsps
    end subroutine ps_col


    !> d (lignes 5..9) ↔ s (colonne nu_s). v est v_ba ; signe (+1) car
    !> (-1)^(l_A+l_B) = (-1)^(2+0) = +1.
    pure subroutine ds_col(v_ba, l_, m_, n_, blk, nu_s)
        real(wp), intent(in)    :: v_ba(10), l_, m_, n_
        real(wp), intent(inout) :: blk(:,:)
        integer,  intent(in)    :: nu_s
        real(wp) :: vsds, alpha, beta, z2
        vsds  = v_ba(8)
        alpha = l_*l_ + m_*m_
        beta  = l_*l_ - m_*m_
        z2    = n_*n_
        blk(5, nu_s) =  SQRT3 * l_ * m_       * vsds
        blk(6, nu_s) =  SQRT3 * m_ * n_       * vsds
        blk(7, nu_s) =  SQRT3 * n_ * l_       * vsds
        blk(8, nu_s) =  0.5_wp * SQRT3 * beta * vsds
        blk(9, nu_s) = (z2 - 0.5_wp * alpha)  * vsds
    end subroutine ds_col


    !> Bloc pp (px, py, pz) × (px, py, pz).
    pure subroutine pp_block(v, l_, m_, n_, blk)
        real(wp), intent(in)    :: v(10), l_, m_, n_
        real(wp), intent(inout) :: blk(:,:)
        real(wp) :: vpps, vppp
        vpps = v(6); vppp = v(7)
        blk(2,2) = l_*l_ * vpps + (1.0_wp - l_*l_) * vppp
        blk(2,3) = l_*m_ * (vpps - vppp)
        blk(2,4) = l_*n_ * (vpps - vppp)
        blk(3,2) = blk(2,3)
        blk(3,3) = m_*m_ * vpps + (1.0_wp - m_*m_) * vppp
        blk(3,4) = m_*n_ * (vpps - vppp)
        blk(4,2) = blk(2,4)
        blk(4,3) = blk(3,4)
        blk(4,4) = n_*n_ * vpps + (1.0_wp - n_*n_) * vppp
    end subroutine pp_block


    !> Bloc p (A) ↔ d (B). Lignes 2..4, colonnes 5..9. Utilise v_ab.
    pure subroutine pd_block(v, l_, m_, n_, blk)
        real(wp), intent(in)    :: v(10), l_, m_, n_
        real(wp), intent(inout) :: blk(:,:)
        real(wp) :: vpds, vpdp, alpha, beta, z2

        vpds  = v(4); vpdp = v(5)
        alpha = l_*l_ + m_*m_
        beta  = l_*l_ - m_*m_
        z2    = n_*n_

        ! px ↔ d
        blk(2,5) = SQRT3 * l_*l_ * m_              * vpds + m_ * (1.0_wp - 2.0_wp*l_*l_) * vpdp
        blk(2,6) = SQRT3 * l_ * m_ * n_            * vpds - 2.0_wp * l_ * m_ * n_        * vpdp
        blk(2,7) = SQRT3 * l_*l_ * n_              * vpds + n_ * (1.0_wp - 2.0_wp*l_*l_) * vpdp
        blk(2,8) = 0.5_wp * SQRT3 * l_ * beta      * vpds + l_ * (1.0_wp - beta)         * vpdp
        blk(2,9) = l_ * (z2 - 0.5_wp*alpha)        * vpds - SQRT3 * l_ * z2              * vpdp

        ! py ↔ d
        blk(3,5) = SQRT3 * m_*m_ * l_              * vpds + l_ * (1.0_wp - 2.0_wp*m_*m_) * vpdp
        blk(3,6) = SQRT3 * m_*m_ * n_              * vpds + n_ * (1.0_wp - 2.0_wp*m_*m_) * vpdp
        blk(3,7) = SQRT3 * l_ * m_ * n_            * vpds - 2.0_wp * l_ * m_ * n_        * vpdp
        blk(3,8) = 0.5_wp * SQRT3 * m_ * beta      * vpds - m_ * (1.0_wp + beta)         * vpdp
        blk(3,9) = m_ * (z2 - 0.5_wp*alpha)        * vpds - SQRT3 * m_ * z2              * vpdp

        ! pz ↔ d
        blk(4,5) = SQRT3 * l_ * m_ * n_            * vpds - 2.0_wp * l_ * m_ * n_        * vpdp
        blk(4,6) = SQRT3 * n_*n_ * m_              * vpds + m_ * (1.0_wp - 2.0_wp*n_*n_) * vpdp
        blk(4,7) = SQRT3 * n_*n_ * l_              * vpds + l_ * (1.0_wp - 2.0_wp*n_*n_) * vpdp
        blk(4,8) = 0.5_wp * SQRT3 * n_ * beta      * vpds - n_ * beta                    * vpdp
        blk(4,9) = n_ * (z2 - 0.5_wp*alpha)        * vpds + SQRT3 * n_ * alpha           * vpdp
    end subroutine pd_block


    !> Bloc d (A) ↔ p (B). v est v_ba ; signe (-1)^(l_A+l_B) = -1.
    pure subroutine dp_block(v_ba, l_, m_, n_, blk)
        real(wp), intent(in)    :: v_ba(10), l_, m_, n_
        real(wp), intent(inout) :: blk(:,:)
        real(wp) :: tblk(9, 9)
        integer :: i, j

        tblk = 0.0_wp
        ! Calcule pd avec v_ba en réutilisant la formule, puis transpose
        ! et applique le signe (-1)^(l_A+l_B) = -1.
        call pd_block(v_ba, l_, m_, n_, tblk)
        do i = 5, 9
            do j = 2, 4
                blk(i, j) = -tblk(j, i)
            end do
        end do
    end subroutine dp_block


    !> Bloc dd (5..9) × (5..9).
    pure subroutine dd_block(v, l_, m_, n_, blk)
        real(wp), intent(in)    :: v(10), l_, m_, n_
        real(wp), intent(inout) :: blk(:,:)
        real(wp) :: vdds, vddp, vddd
        real(wp) :: alpha, beta, z2, l2, m2, n2, lm, mn, nl

        vdds = v(1); vddp = v(2); vddd = v(3)
        l2 = l_*l_; m2 = m_*m_; n2 = n_*n_
        lm = l_*m_; mn = m_*n_; nl = n_*l_
        alpha = l2 + m2
        beta  = l2 - m2
        z2    = n2

        ! dxy / dxy
        blk(5,5) = 3.0_wp*l2*m2 * vdds + (alpha - 4.0_wp*l2*m2) * vddp &
                 + (z2 + l2*m2)        * vddd
        ! dxy / dyz
        blk(5,6) = 3.0_wp*l_*m2*n_ * vdds + l_*n_*(1.0_wp - 4.0_wp*m2) * vddp &
                 + l_*n_*(m2 - 1.0_wp)    * vddd
        ! dxy / dzx
        blk(5,7) = 3.0_wp*l2*m_*n_ * vdds + m_*n_*(1.0_wp - 4.0_wp*l2) * vddp &
                 + m_*n_*(l2 - 1.0_wp)    * vddd
        ! dxy / dx2-y2
        blk(5,8) = 1.5_wp*lm*beta  * vdds - 2.0_wp*lm*beta * vddp &
                 + 0.5_wp*lm*beta  * vddd
        ! dxy / d3z2-r2
        blk(5,9) = SQRT3*lm*(z2 - 0.5_wp*alpha) * vdds - 2.0_wp*SQRT3*lm*z2 * vddp &
                 + 0.5_wp*SQRT3*lm*(1.0_wp + z2) * vddd

        ! dyz / dyz
        blk(6,6) = 3.0_wp*m2*n2 * vdds + (m2 + n2 - 4.0_wp*m2*n2) * vddp &
                 + (l2 + m2*n2) * vddd
        ! dyz / dzx
        blk(6,7) = 3.0_wp*l_*m_*n2 * vdds + l_*m_*(1.0_wp - 4.0_wp*n2) * vddp &
                 + l_*m_*(n2 - 1.0_wp)    * vddd
        ! dyz / dx2-y2
        blk(6,8) = 1.5_wp*mn*beta * vdds - mn*(1.0_wp + 2.0_wp*beta) * vddp &
                 + mn*(1.0_wp + 0.5_wp*beta) * vddd
        ! dyz / d3z2-r2
        blk(6,9) = SQRT3*mn*(z2 - 0.5_wp*alpha) * vdds + SQRT3*mn*(alpha - z2) * vddp &
                 - 0.5_wp*SQRT3*mn*alpha        * vddd

        ! dzx / dzx
        blk(7,7) = 3.0_wp*n2*l2 * vdds + (n2 + l2 - 4.0_wp*n2*l2) * vddp &
                 + (m2 + n2*l2) * vddd
        ! dzx / dx2-y2
        blk(7,8) = 1.5_wp*nl*beta * vdds + nl*(1.0_wp - 2.0_wp*beta) * vddp &
                 - nl*(1.0_wp - 0.5_wp*beta) * vddd
        ! dzx / d3z2-r2
        blk(7,9) = SQRT3*nl*(z2 - 0.5_wp*alpha) * vdds + SQRT3*nl*(alpha - z2) * vddp &
                 - 0.5_wp*SQRT3*nl*alpha        * vddd

        ! dx2-y2 / dx2-y2
        blk(8,8) = 0.75_wp*beta*beta * vdds + (alpha - beta*beta) * vddp &
                 + (z2 + 0.25_wp*beta*beta) * vddd
        ! dx2-y2 / d3z2-r2
        blk(8,9) = 0.5_wp*SQRT3*beta*(z2 - 0.5_wp*alpha) * vdds &
                 - SQRT3*z2*beta * vddp &
                 + 0.25_wp*SQRT3*(1.0_wp + z2)*beta * vddd

        ! d3z2-r2 / d3z2-r2
        blk(9,9) = (z2 - 0.5_wp*alpha)**2 * vdds + 3.0_wp*z2*alpha * vddp &
                 + 0.75_wp*alpha*alpha    * vddd

        ! Symétrise le sous-bloc dd
        blk(6,5) = blk(5,6)
        blk(7,5) = blk(5,7); blk(7,6) = blk(6,7)
        blk(8,5) = blk(5,8); blk(8,6) = blk(6,8); blk(8,7) = blk(7,8)
        blk(9,5) = blk(5,9); blk(9,6) = blk(6,9); blk(9,7) = blk(7,9)
        blk(9,8) = blk(8,9)
    end subroutine dd_block
end module skrot
