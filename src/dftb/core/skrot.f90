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
!> Les directions sont exprimées en cosinus directeurs (x, y, z) avec
!> α = x²+y², β = x²−y². Référence : Slater & Koster Phys. Rev. 94,
!> 1498 (1954), Table I.
module skrot
    use kinds, only: wp
    implicit none
    private

    public :: transform_sk, transform_sk_opt

    real(wp), parameter :: SQRT3 = 1.7320508075688772_wp

contains

    !> Retourne l'élément (mu, nu) du bloc rotationné Slater-Koster pour
    !> une paire d'atomes A (en 0) et B (en R), où :
    !>   x, y, z : cosinus directeurs du vecteur unitaire A → B
    !>   mu, nu  : indices d'orbitales (1..9) sur A et B
    !>   sk_ab   : 10 intégrales SK lues A → B
    !>   sk_ba   : 10 intégrales SK lues B → A (utilisées pour les
    !>             composantes asymétriques où l_A ≠ l_B, avec
    !>             le signe (-1)^(l_A+l_B))
    pure function transform_sk(x, y, z, mu, nu, sk_ab, sk_ba) result(v)
        real(wp), intent(in) :: x, y, z
        integer,  intent(in) :: mu, nu
        real(wp), intent(in) :: sk_ab(10), sk_ba(10)
        real(wp) :: v
        real(wp) :: alpha, beta, z2

        alpha = x*x + y*y
        beta  = x*x - y*y
        z2    = z*z
        v     = 0.0_wp

        select case (mu)

        ! --------------------------------------------------------------
        case (1)  ! s
            select case (nu)
            case (1)  ! s ↔ s
                v = sk_ab(10)
            case (2)  ! s ↔ px
                v = x * sk_ab(9)
            case (3)  ! s ↔ py
                v = y * sk_ab(9)
            case (4)  ! s ↔ pz
                v = z * sk_ab(9)
            case (5)  ! s ↔ dxy
                v = SQRT3 * x * y * sk_ab(8)
            case (6)  ! s ↔ dyz
                v = SQRT3 * y * z * sk_ab(8)
            case (7)  ! s ↔ dzx
                v = SQRT3 * z * x * sk_ab(8)
            case (8)  ! s ↔ dx2-y2
                v = 0.5_wp * SQRT3 * beta * sk_ab(8)
            case (9)  ! s ↔ d3z2-r2
                v = (z2 - 0.5_wp * alpha) * sk_ab(8)
            end select

        ! --------------------------------------------------------------
        case (2)  ! px
            select case (nu)
            case (1)  ! px ↔ s : signe (-1)^(1+0) = -1, sk_ba
                v = -x * sk_ba(9)
            case (2)  ! px ↔ px
                v = x*x * sk_ab(6) + (1.0_wp - x*x) * sk_ab(7)
            case (3)  ! px ↔ py
                v = x*y * (sk_ab(6) - sk_ab(7))
            case (4)  ! px ↔ pz
                v = x*z * (sk_ab(6) - sk_ab(7))
            case (5)  ! px ↔ dxy
                v = SQRT3 * x*x * y * sk_ab(4) &
                    + y * (1.0_wp - 2.0_wp*x*x) * sk_ab(5)
            case (6)  ! px ↔ dyz
                v = SQRT3 * x * y * z * sk_ab(4) &
                    - 2.0_wp * x * y * z * sk_ab(5)
            case (7)  ! px ↔ dzx
                v = SQRT3 * x*x * z * sk_ab(4) &
                    + z * (1.0_wp - 2.0_wp*x*x) * sk_ab(5)
            case (8)  ! px ↔ dx2-y2
                v = 0.5_wp * SQRT3 * x * beta * sk_ab(4) &
                    + x * (1.0_wp - beta) * sk_ab(5)
            case (9)  ! px ↔ d3z2-r2
                v = x * (z2 - 0.5_wp*alpha) * sk_ab(4) &
                    - SQRT3 * x * z2 * sk_ab(5)
            end select

        ! --------------------------------------------------------------
        case (3)  ! py
            select case (nu)
            case (1)  ! py ↔ s : signe -1, sk_ba
                v = -y * sk_ba(9)
            case (2)  ! py ↔ px
                v = x*y * (sk_ab(6) - sk_ab(7))
            case (3)  ! py ↔ py
                v = y*y * sk_ab(6) + (1.0_wp - y*y) * sk_ab(7)
            case (4)  ! py ↔ pz
                v = y*z * (sk_ab(6) - sk_ab(7))
            case (5)  ! py ↔ dxy
                v = SQRT3 * y*y * x * sk_ab(4) &
                    + x * (1.0_wp - 2.0_wp*y*y) * sk_ab(5)
            case (6)  ! py ↔ dyz
                v = SQRT3 * y*y * z * sk_ab(4) &
                    + z * (1.0_wp - 2.0_wp*y*y) * sk_ab(5)
            case (7)  ! py ↔ dzx
                v = SQRT3 * x * y * z * sk_ab(4) &
                    - 2.0_wp * x * y * z * sk_ab(5)
            case (8)  ! py ↔ dx2-y2
                v = 0.5_wp * SQRT3 * y * beta * sk_ab(4) &
                    - y * (1.0_wp + beta) * sk_ab(5)
            case (9)  ! py ↔ d3z2-r2
                v = y * (z2 - 0.5_wp*alpha) * sk_ab(4) &
                    - SQRT3 * y * z2 * sk_ab(5)
            end select

        ! --------------------------------------------------------------
        case (4)  ! pz
            select case (nu)
            case (1)  ! pz ↔ s : signe -1, sk_ba
                v = -z * sk_ba(9)
            case (2)  ! pz ↔ px
                v = x*z * (sk_ab(6) - sk_ab(7))
            case (3)  ! pz ↔ py
                v = y*z * (sk_ab(6) - sk_ab(7))
            case (4)  ! pz ↔ pz
                v = z*z * sk_ab(6) + (1.0_wp - z*z) * sk_ab(7)
            case (5)  ! pz ↔ dxy
                v = SQRT3 * x * y * z * sk_ab(4) &
                    - 2.0_wp * x * y * z * sk_ab(5)
            case (6)  ! pz ↔ dyz
                v = SQRT3 * z*z * y * sk_ab(4) &
                    + y * (1.0_wp - 2.0_wp*z*z) * sk_ab(5)
            case (7)  ! pz ↔ dzx
                v = SQRT3 * z*z * x * sk_ab(4) &
                    + x * (1.0_wp - 2.0_wp*z*z) * sk_ab(5)
            case (8)  ! pz ↔ dx2-y2
                v = 0.5_wp * SQRT3 * z * beta * sk_ab(4) &
                    - z * beta * sk_ab(5)
            case (9)  ! pz ↔ d3z2-r2
                v = z * (z2 - 0.5_wp*alpha) * sk_ab(4) &
                    + SQRT3 * z * alpha * sk_ab(5)
            end select

        ! --------------------------------------------------------------
        case (5)  ! dxy
            select case (nu)
            case (1)  ! dxy ↔ s : signe (-1)^(2+0) = +1, sk_ba
                v = SQRT3 * x * y * sk_ba(8)
            case (2)  ! dxy ↔ px : signe (-1)^(2+1) = -1, sk_ba
                v = -(SQRT3 * x*x * y * sk_ba(4) &
                      + y * (1.0_wp - 2.0_wp*x*x) * sk_ba(5))
            case (3)  ! dxy ↔ py
                v = -(SQRT3 * y*y * x * sk_ba(4) &
                      + x * (1.0_wp - 2.0_wp*y*y) * sk_ba(5))
            case (4)  ! dxy ↔ pz
                v = -(SQRT3 * x * y * z * sk_ba(4) &
                      - 2.0_wp * x * y * z * sk_ba(5))
            case (5)  ! dxy ↔ dxy
                v = 3.0_wp*x*x*y*y * sk_ab(1) &
                    + (alpha - 4.0_wp*x*x*y*y) * sk_ab(2) &
                    + (z2 + x*x*y*y) * sk_ab(3)
            case (6)  ! dxy ↔ dyz
                v = 3.0_wp*x*y*y*z * sk_ab(1) &
                    + x*z*(1.0_wp - 4.0_wp*y*y) * sk_ab(2) &
                    + x*z*(y*y - 1.0_wp) * sk_ab(3)
            case (7)  ! dxy ↔ dzx
                v = 3.0_wp*x*x*y*z * sk_ab(1) &
                    + y*z*(1.0_wp - 4.0_wp*x*x) * sk_ab(2) &
                    + y*z*(x*x - 1.0_wp) * sk_ab(3)
            case (8)  ! dxy ↔ dx2-y2
                v = 1.5_wp * x*y * beta * sk_ab(1) &
                    - 2.0_wp * x*y * beta * sk_ab(2) &
                    + 0.5_wp * x*y * beta * sk_ab(3)
            case (9)  ! dxy ↔ d3z2-r2
                v = SQRT3 * x*y * (z2 - 0.5_wp*alpha) * sk_ab(1) &
                    - 2.0_wp * SQRT3 * x*y * z2 * sk_ab(2) &
                    + 0.5_wp * SQRT3 * x*y * (1.0_wp + z2) * sk_ab(3)
            end select

        ! --------------------------------------------------------------
        case (6)  ! dyz
            select case (nu)
            case (1)  ! dyz ↔ s : signe +1, sk_ba
                v = SQRT3 * y * z * sk_ba(8)
            case (2)  ! dyz ↔ px : signe -1, sk_ba
                v = -(SQRT3 * x * y * z * sk_ba(4) &
                      - 2.0_wp * x * y * z * sk_ba(5))
            case (3)  ! dyz ↔ py
                v = -(SQRT3 * y*y * z * sk_ba(4) &
                      + z * (1.0_wp - 2.0_wp*y*y) * sk_ba(5))
            case (4)  ! dyz ↔ pz
                v = -(SQRT3 * z*z * y * sk_ba(4) &
                      + y * (1.0_wp - 2.0_wp*z*z) * sk_ba(5))
            case (5)  ! dyz ↔ dxy : symétrie => même formule que dxy/dyz
                v = 3.0_wp*x*y*y*z * sk_ab(1) &
                    + x*z*(1.0_wp - 4.0_wp*y*y) * sk_ab(2) &
                    + x*z*(y*y - 1.0_wp) * sk_ab(3)
            case (6)  ! dyz ↔ dyz
                v = 3.0_wp*y*y*z*z * sk_ab(1) &
                    + (y*y + z*z - 4.0_wp*y*y*z*z) * sk_ab(2) &
                    + (x*x + y*y*z*z) * sk_ab(3)
            case (7)  ! dyz ↔ dzx
                v = 3.0_wp*x*y*z*z * sk_ab(1) &
                    + x*y*(1.0_wp - 4.0_wp*z*z) * sk_ab(2) &
                    + x*y*(z*z - 1.0_wp) * sk_ab(3)
            case (8)  ! dyz ↔ dx2-y2
                v = 1.5_wp * y*z * beta * sk_ab(1) &
                    - y*z*(1.0_wp + 2.0_wp*beta) * sk_ab(2) &
                    + y*z*(1.0_wp + 0.5_wp*beta) * sk_ab(3)
            case (9)  ! dyz ↔ d3z2-r2
                v = SQRT3 * y*z * (z2 - 0.5_wp*alpha) * sk_ab(1) &
                    + SQRT3 * y*z * (alpha - z2) * sk_ab(2) &
                    - 0.5_wp * SQRT3 * y*z * alpha * sk_ab(3)
            end select

        ! --------------------------------------------------------------
        case (7)  ! dzx
            select case (nu)
            case (1)  ! dzx ↔ s
                v = SQRT3 * z * x * sk_ba(8)
            case (2)  ! dzx ↔ px
                v = -(SQRT3 * x*x * z * sk_ba(4) &
                      + z * (1.0_wp - 2.0_wp*x*x) * sk_ba(5))
            case (3)  ! dzx ↔ py
                v = -(SQRT3 * x * y * z * sk_ba(4) &
                      - 2.0_wp * x * y * z * sk_ba(5))
            case (4)  ! dzx ↔ pz
                v = -(SQRT3 * z*z * x * sk_ba(4) &
                      + x * (1.0_wp - 2.0_wp*z*z) * sk_ba(5))
            case (5)  ! dzx ↔ dxy
                v = 3.0_wp*x*x*y*z * sk_ab(1) &
                    + y*z*(1.0_wp - 4.0_wp*x*x) * sk_ab(2) &
                    + y*z*(x*x - 1.0_wp) * sk_ab(3)
            case (6)  ! dzx ↔ dyz
                v = 3.0_wp*x*y*z*z * sk_ab(1) &
                    + x*y*(1.0_wp - 4.0_wp*z*z) * sk_ab(2) &
                    + x*y*(z*z - 1.0_wp) * sk_ab(3)
            case (7)  ! dzx ↔ dzx
                v = 3.0_wp*z*z*x*x * sk_ab(1) &
                    + (z*z + x*x - 4.0_wp*z*z*x*x) * sk_ab(2) &
                    + (y*y + z*z*x*x) * sk_ab(3)
            case (8)  ! dzx ↔ dx2-y2
                v = 1.5_wp * z*x * beta * sk_ab(1) &
                    + z*x * (1.0_wp - 2.0_wp*beta) * sk_ab(2) &
                    - z*x * (1.0_wp - 0.5_wp*beta) * sk_ab(3)
            case (9)  ! dzx ↔ d3z2-r2
                v = SQRT3 * z*x * (z2 - 0.5_wp*alpha) * sk_ab(1) &
                    + SQRT3 * z*x * (alpha - z2) * sk_ab(2) &
                    - 0.5_wp * SQRT3 * z*x * alpha * sk_ab(3)
            end select

        ! --------------------------------------------------------------
        case (8)  ! dx2-y2
            select case (nu)
            case (1)  ! dx2-y2 ↔ s
                v = 0.5_wp * SQRT3 * beta * sk_ba(8)
            case (2)  ! dx2-y2 ↔ px
                v = -(0.5_wp * SQRT3 * x * beta * sk_ba(4) &
                      + x * (1.0_wp - beta) * sk_ba(5))
            case (3)  ! dx2-y2 ↔ py
                v = -(0.5_wp * SQRT3 * y * beta * sk_ba(4) &
                      - y * (1.0_wp + beta) * sk_ba(5))
            case (4)  ! dx2-y2 ↔ pz
                v = -(0.5_wp * SQRT3 * z * beta * sk_ba(4) &
                      - z * beta * sk_ba(5))
            case (5)  ! dx2-y2 ↔ dxy
                v = 1.5_wp * x*y * beta * sk_ab(1) &
                    - 2.0_wp * x*y * beta * sk_ab(2) &
                    + 0.5_wp * x*y * beta * sk_ab(3)
            case (6)  ! dx2-y2 ↔ dyz
                v = 1.5_wp * y*z * beta * sk_ab(1) &
                    - y*z*(1.0_wp + 2.0_wp*beta) * sk_ab(2) &
                    + y*z*(1.0_wp + 0.5_wp*beta) * sk_ab(3)
            case (7)  ! dx2-y2 ↔ dzx
                v = 1.5_wp * z*x * beta * sk_ab(1) &
                    + z*x * (1.0_wp - 2.0_wp*beta) * sk_ab(2) &
                    - z*x * (1.0_wp - 0.5_wp*beta) * sk_ab(3)
            case (8)  ! dx2-y2 ↔ dx2-y2
                v = 0.75_wp * beta*beta * sk_ab(1) &
                    + (alpha - beta*beta) * sk_ab(2) &
                    + (z2 + 0.25_wp * beta*beta) * sk_ab(3)
            case (9)  ! dx2-y2 ↔ d3z2-r2
                v = 0.5_wp * SQRT3 * beta * (z2 - 0.5_wp*alpha) * sk_ab(1) &
                    - SQRT3 * z2 * beta * sk_ab(2) &
                    + 0.25_wp * SQRT3 * (1.0_wp + z2) * beta * sk_ab(3)
            end select

        ! --------------------------------------------------------------
        case (9)  ! d3z2-r2
            select case (nu)
            case (1)  ! d3z2-r2 ↔ s
                v = (z2 - 0.5_wp * alpha) * sk_ba(8)
            case (2)  ! d3z2-r2 ↔ px
                v = -(x * (z2 - 0.5_wp*alpha) * sk_ba(4) &
                      - SQRT3 * x * z2 * sk_ba(5))
            case (3)  ! d3z2-r2 ↔ py
                v = -(y * (z2 - 0.5_wp*alpha) * sk_ba(4) &
                      - SQRT3 * y * z2 * sk_ba(5))
            case (4)  ! d3z2-r2 ↔ pz
                v = -(z * (z2 - 0.5_wp*alpha) * sk_ba(4) &
                      + SQRT3 * z * alpha * sk_ba(5))
            case (5)  ! d3z2-r2 ↔ dxy
                v = SQRT3 * x*y * (z2 - 0.5_wp*alpha) * sk_ab(1) &
                    - 2.0_wp * SQRT3 * x*y * z2 * sk_ab(2) &
                    + 0.5_wp * SQRT3 * x*y * (1.0_wp + z2) * sk_ab(3)
            case (6)  ! d3z2-r2 ↔ dyz
                v = SQRT3 * y*z * (z2 - 0.5_wp*alpha) * sk_ab(1) &
                    + SQRT3 * y*z * (alpha - z2) * sk_ab(2) &
                    - 0.5_wp * SQRT3 * y*z * alpha * sk_ab(3)
            case (7)  ! d3z2-r2 ↔ dzx
                v = SQRT3 * z*x * (z2 - 0.5_wp*alpha) * sk_ab(1) &
                    + SQRT3 * z*x * (alpha - z2) * sk_ab(2) &
                    - 0.5_wp * SQRT3 * z*x * alpha * sk_ab(3)
            case (8)  ! d3z2-r2 ↔ dx2-y2
                v = 0.5_wp * SQRT3 * beta * (z2 - 0.5_wp*alpha) * sk_ab(1) &
                    - SQRT3 * z2 * beta * sk_ab(2) &
                    + 0.25_wp * SQRT3 * (1.0_wp + z2) * beta * sk_ab(3)
            case (9)  ! d3z2-r2 ↔ d3z2-r2
                v = (z2 - 0.5_wp*alpha)**2 * sk_ab(1) &
                    + 3.0_wp * z2 * alpha * sk_ab(2) &
                    + 0.75_wp * alpha*alpha * sk_ab(3)
            end select

        end select
    end function transform_sk


    !> Variante optimisée de `transform_sk` : précalcule les monômes
    !> géométriques et factorise les coefficients (σ/π/δ) avant
    !> multiplication par sk_ab / sk_ba. Le résultat numérique est
    !> identique à `transform_sk` (à la précision flottante près).
    pure function transform_sk_opt(x, y, z, mu, nu, sk_ab, sk_ba) result(v)
        real(wp), intent(in) :: x, y, z
        integer,  intent(in) :: mu, nu
        real(wp), intent(in) :: sk_ab(10), sk_ba(10)
        real(wp) :: v
        real(wp) :: x2, y2, z2, xy, yz, zx
        real(wp) :: alpha, beta, zh, x2y2, xyz
        real(wp) :: cs, cp, cd
        integer  :: key

        x2 = x*x; y2 = y*y; z2 = z*z
        xy = x*y; yz = y*z; zx = z*x
        alpha = x2 + y2
        beta  = x2 - y2
        zh    = z2 - 0.5_wp * alpha
        x2y2  = x2 * y2
        xyz   = x * y * z

        v   = 0.0_wp
        key = 9 * (mu - 1) + nu

        select case (key)

        ! ---- mu = 1 (s) ----
        case (1)            ! s, s
            v = sk_ab(10)
        case (2)            ! s, px
            v = x * sk_ab(9)
        case (3)            ! s, py
            v = y * sk_ab(9)
        case (4)            ! s, pz
            v = z * sk_ab(9)
        case (5)            ! s, dxy
            v = SQRT3 * xy * sk_ab(8)
        case (6)            ! s, dyz
            v = SQRT3 * yz * sk_ab(8)
        case (7)            ! s, dzx
            v = SQRT3 * zx * sk_ab(8)
        case (8)            ! s, dx2-y2
            v = 0.5_wp * SQRT3 * beta * sk_ab(8)
        case (9)            ! s, d3z2-r2
            v = zh * sk_ab(8)

        ! ---- mu = 2 (px) ----
        case (10)           ! px, s
            v = -x * sk_ba(9)
        case (11)           ! px, px
            v = x2 * sk_ab(6) + (1.0_wp - x2) * sk_ab(7)
        case (12)           ! px, py
            v = xy * (sk_ab(6) - sk_ab(7))
        case (13)           ! px, pz
            v = zx * (sk_ab(6) - sk_ab(7))
        case (14)           ! px, dxy
            cs = SQRT3 * x2 * y
            cp = y * (1.0_wp - 2.0_wp*x2)
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (15)           ! px, dyz
            cs = SQRT3 * xyz
            cp = -2.0_wp * xyz
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (16)           ! px, dzx
            cs = SQRT3 * x2 * z
            cp = z * (1.0_wp - 2.0_wp*x2)
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (17)           ! px, dx2-y2
            cs = 0.5_wp * SQRT3 * x * beta
            cp = x * (1.0_wp - beta)
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (18)           ! px, d3z2-r2
            cs = x * zh
            cp = -SQRT3 * x * z2
            v  = cs * sk_ab(4) + cp * sk_ab(5)

        ! ---- mu = 3 (py) ----
        case (19)           ! py, s
            v = -y * sk_ba(9)
        case (20)           ! py, px
            v = xy * (sk_ab(6) - sk_ab(7))
        case (21)           ! py, py
            v = y2 * sk_ab(6) + (1.0_wp - y2) * sk_ab(7)
        case (22)           ! py, pz
            v = yz * (sk_ab(6) - sk_ab(7))
        case (23)           ! py, dxy
            cs = SQRT3 * y2 * x
            cp = x * (1.0_wp - 2.0_wp*y2)
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (24)           ! py, dyz
            cs = SQRT3 * y2 * z
            cp = z * (1.0_wp - 2.0_wp*y2)
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (25)           ! py, dzx
            cs = SQRT3 * xyz
            cp = -2.0_wp * xyz
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (26)           ! py, dx2-y2
            cs = 0.5_wp * SQRT3 * y * beta
            cp = -y * (1.0_wp + beta)
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (27)           ! py, d3z2-r2
            cs = y * zh
            cp = -SQRT3 * y * z2
            v  = cs * sk_ab(4) + cp * sk_ab(5)

        ! ---- mu = 4 (pz) ----
        case (28)           ! pz, s
            v = -z * sk_ba(9)
        case (29)           ! pz, px
            v = zx * (sk_ab(6) - sk_ab(7))
        case (30)           ! pz, py
            v = yz * (sk_ab(6) - sk_ab(7))
        case (31)           ! pz, pz
            v = z2 * sk_ab(6) + (1.0_wp - z2) * sk_ab(7)
        case (32)           ! pz, dxy
            cs = SQRT3 * xyz
            cp = -2.0_wp * xyz
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (33)           ! pz, dyz
            cs = SQRT3 * z2 * y
            cp = y * (1.0_wp - 2.0_wp*z2)
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (34)           ! pz, dzx
            cs = SQRT3 * z2 * x
            cp = x * (1.0_wp - 2.0_wp*z2)
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (35)           ! pz, dx2-y2
            cs = 0.5_wp * SQRT3 * z * beta
            cp = -z * beta
            v  = cs * sk_ab(4) + cp * sk_ab(5)
        case (36)           ! pz, d3z2-r2
            cs = z * zh
            cp = SQRT3 * z * alpha
            v  = cs * sk_ab(4) + cp * sk_ab(5)

        ! ---- mu = 5 (dxy) ----
        case (37)           ! dxy, s
            v = SQRT3 * xy * sk_ba(8)
        case (38)           ! dxy, px
            cs = SQRT3 * x2 * y
            cp = y * (1.0_wp - 2.0_wp*x2)
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (39)           ! dxy, py
            cs = SQRT3 * y2 * x
            cp = x * (1.0_wp - 2.0_wp*y2)
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (40)           ! dxy, pz
            cs = SQRT3 * xyz
            cp = -2.0_wp * xyz
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (41)           ! dxy, dxy
            cs = 3.0_wp * x2y2
            cp = alpha - 4.0_wp * x2y2
            cd = z2 + x2y2
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (42)           ! dxy, dyz
            cs = 3.0_wp * xy * yz
            cp = zx * (1.0_wp - 4.0_wp*y2)
            cd = zx * (y2 - 1.0_wp)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (43)           ! dxy, dzx
            cs = 3.0_wp * x2 * yz
            cp = yz * (1.0_wp - 4.0_wp*x2)
            cd = yz * (x2 - 1.0_wp)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (44)           ! dxy, dx2-y2
            cs = 1.5_wp * xy * beta
            cp = -2.0_wp * xy * beta
            cd = 0.5_wp * xy * beta
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (45)           ! dxy, d3z2-r2
            cs = SQRT3 * xy * zh
            cp = -2.0_wp * SQRT3 * xy * z2
            cd = 0.5_wp * SQRT3 * xy * (1.0_wp + z2)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)

        ! ---- mu = 6 (dyz) ----
        case (46)           ! dyz, s
            v = SQRT3 * yz * sk_ba(8)
        case (47)           ! dyz, px
            cs = SQRT3 * xyz
            cp = -2.0_wp * xyz
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (48)           ! dyz, py
            cs = SQRT3 * y2 * z
            cp = z * (1.0_wp - 2.0_wp*y2)
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (49)           ! dyz, pz
            cs = SQRT3 * z2 * y
            cp = y * (1.0_wp - 2.0_wp*z2)
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (50)           ! dyz, dxy
            cs = 3.0_wp * xy * yz
            cp = zx * (1.0_wp - 4.0_wp*y2)
            cd = zx * (y2 - 1.0_wp)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (51)           ! dyz, dyz
            cs = 3.0_wp * y2 * z2
            cp = y2 + z2 - 4.0_wp * y2 * z2
            cd = x2 + y2 * z2
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (52)           ! dyz, dzx
            cs = 3.0_wp * xy * z2
            cp = xy * (1.0_wp - 4.0_wp*z2)
            cd = xy * (z2 - 1.0_wp)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (53)           ! dyz, dx2-y2
            cs = 1.5_wp * yz * beta
            cp = -yz * (1.0_wp + 2.0_wp*beta)
            cd = yz * (1.0_wp + 0.5_wp*beta)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (54)           ! dyz, d3z2-r2
            cs = SQRT3 * yz * zh
            cp = SQRT3 * yz * (alpha - z2)
            cd = -0.5_wp * SQRT3 * yz * alpha
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)

        ! ---- mu = 7 (dzx) ----
        case (55)           ! dzx, s
            v = SQRT3 * zx * sk_ba(8)
        case (56)           ! dzx, px
            cs = SQRT3 * x2 * z
            cp = z * (1.0_wp - 2.0_wp*x2)
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (57)           ! dzx, py
            cs = SQRT3 * xyz
            cp = -2.0_wp * xyz
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (58)           ! dzx, pz
            cs = SQRT3 * z2 * x
            cp = x * (1.0_wp - 2.0_wp*z2)
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (59)           ! dzx, dxy
            cs = 3.0_wp * x2 * yz
            cp = yz * (1.0_wp - 4.0_wp*x2)
            cd = yz * (x2 - 1.0_wp)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (60)           ! dzx, dyz
            cs = 3.0_wp * xy * z2
            cp = xy * (1.0_wp - 4.0_wp*z2)
            cd = xy * (z2 - 1.0_wp)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (61)           ! dzx, dzx
            cs = 3.0_wp * z2 * x2
            cp = z2 + x2 - 4.0_wp * z2 * x2
            cd = y2 + z2 * x2
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (62)           ! dzx, dx2-y2
            cs = 1.5_wp * zx * beta
            cp = zx * (1.0_wp - 2.0_wp*beta)
            cd = -zx * (1.0_wp - 0.5_wp*beta)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (63)           ! dzx, d3z2-r2
            cs = SQRT3 * zx * zh
            cp = SQRT3 * zx * (alpha - z2)
            cd = -0.5_wp * SQRT3 * zx * alpha
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)

        ! ---- mu = 8 (dx2-y2) ----
        case (64)           ! dx2-y2, s
            v = 0.5_wp * SQRT3 * beta * sk_ba(8)
        case (65)           ! dx2-y2, px
            cs = 0.5_wp * SQRT3 * x * beta
            cp = x * (1.0_wp - beta)
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (66)           ! dx2-y2, py
            cs = 0.5_wp * SQRT3 * y * beta
            cp = -y * (1.0_wp + beta)
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (67)           ! dx2-y2, pz
            cs = 0.5_wp * SQRT3 * z * beta
            cp = -z * beta
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (68)           ! dx2-y2, dxy
            cs = 1.5_wp * xy * beta
            cp = -2.0_wp * xy * beta
            cd = 0.5_wp * xy * beta
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (69)           ! dx2-y2, dyz
            cs = 1.5_wp * yz * beta
            cp = -yz * (1.0_wp + 2.0_wp*beta)
            cd = yz * (1.0_wp + 0.5_wp*beta)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (70)           ! dx2-y2, dzx
            cs = 1.5_wp * zx * beta
            cp = zx * (1.0_wp - 2.0_wp*beta)
            cd = -zx * (1.0_wp - 0.5_wp*beta)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (71)           ! dx2-y2, dx2-y2
            cs = 0.75_wp * beta * beta
            cp = alpha - beta * beta
            cd = z2 + 0.25_wp * beta * beta
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (72)           ! dx2-y2, d3z2-r2
            cs = 0.5_wp * SQRT3 * beta * zh
            cp = -SQRT3 * z2 * beta
            cd = 0.25_wp * SQRT3 * (1.0_wp + z2) * beta
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)

        ! ---- mu = 9 (d3z2-r2) ----
        case (73)           ! d3z2-r2, s
            v = zh * sk_ba(8)
        case (74)           ! d3z2-r2, px
            cs = x * zh
            cp = -SQRT3 * x * z2
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (75)           ! d3z2-r2, py
            cs = y * zh
            cp = -SQRT3 * y * z2
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (76)           ! d3z2-r2, pz
            cs = z * zh
            cp = SQRT3 * z * alpha
            v  = -(cs * sk_ba(4) + cp * sk_ba(5))
        case (77)           ! d3z2-r2, dxy
            cs = SQRT3 * xy * zh
            cp = -2.0_wp * SQRT3 * xy * z2
            cd = 0.5_wp * SQRT3 * xy * (1.0_wp + z2)
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (78)           ! d3z2-r2, dyz
            cs = SQRT3 * yz * zh
            cp = SQRT3 * yz * (alpha - z2)
            cd = -0.5_wp * SQRT3 * yz * alpha
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (79)           ! d3z2-r2, dzx
            cs = SQRT3 * zx * zh
            cp = SQRT3 * zx * (alpha - z2)
            cd = -0.5_wp * SQRT3 * zx * alpha
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (80)           ! d3z2-r2, dx2-y2
            cs = 0.5_wp * SQRT3 * beta * zh
            cp = -SQRT3 * z2 * beta
            cd = 0.25_wp * SQRT3 * (1.0_wp + z2) * beta
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)
        case (81)           ! d3z2-r2, d3z2-r2
            cs = zh * zh
            cp = 3.0_wp * z2 * alpha
            cd = 0.75_wp * alpha * alpha
            v  = cs * sk_ab(1) + cp * sk_ab(2) + cd * sk_ab(3)

        end select
    end function transform_sk_opt

end module skrot
