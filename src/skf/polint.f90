!> Interpolation polynomiale par schéma de Neville (Numerical Recipes).
!>
!> Fournit `polint` (valeur + estimation d'erreur) et `dpolint` (dérivée
!> première par différence finie centrée). L'ordre de la fenêtre est
!> figé à `n_order_polint = 8`.
module polint_mod
    use kinds, only: wp
    implicit none
    private

    integer, parameter, public :: n_order_polint = 8

    public :: polint, dpolint

contains

    !> Interpole y = P(x) par schéma de Neville sur (xa, ya), tableaux de
    !> taille n_order_polint. Retourne `dy` : estimation d'erreur (et non
    !> dérivée), conformément à la recette originale.
    pure subroutine polint(xa, ya, x, y, dy)
        real(wp), intent(in)  :: xa(n_order_polint), ya(n_order_polint)
        real(wp), intent(in)  :: x
        real(wp), intent(out) :: y, dy

        real(wp) :: c(n_order_polint), d(n_order_polint)
        real(wp) :: dif, dift, ho, hp, w, den
        integer  :: i, m, ns, n

        n   = n_order_polint
        ns  = 1
        dif = abs(x - xa(1))
        do i = 1, n
            dift = abs(x - xa(i))
            if (dift < dif) then
                ns  = i
                dif = dift
            end if
            c(i) = ya(i)
            d(i) = ya(i)
        end do

        y  = ya(ns)
        dy = 0.0_wp
        ns = ns - 1
        do m = 1, n - 1
            do i = 1, n - m
                ho  = xa(i)     - x
                hp  = xa(i + m) - x
                w   = c(i + 1) - d(i)
                den = ho - hp
                if (den == 0.0_wp) then
                    y  = 0.0_wp
                    dy = 1.0e9_wp
                    return
                end if
                den   = w / den
                d(i)  = hp * den
                c(i)  = ho * den
            end do
            if (2 * ns < n - m) then
                dy = c(ns + 1)
            else
                dy = d(ns)
                ns = ns - 1
            end if
            y = y + dy
        end do
    end subroutine polint


    !> Dérivée première de l'interpolant polint par différence finie
    !> centrée. h fixé à 1e-4 : compromis générique pour des grilles SKF
    !> de pas typique 0.02 bohr.
    pure subroutine dpolint(xa, ya, x, dydx)
        real(wp), intent(in)  :: xa(n_order_polint), ya(n_order_polint)
        real(wp), intent(in)  :: x
        real(wp), intent(out) :: dydx
        real(wp), parameter   :: H_STEP = 1.0e-4_wp
        real(wp) :: yp, ym, dyp, dym

        call polint(xa, ya, x + H_STEP, yp, dyp)
        call polint(xa, ya, x - H_STEP, ym, dym)
        dydx = (yp - ym) / (2.0_wp * H_STEP)
    end subroutine dpolint
end module polint_mod
