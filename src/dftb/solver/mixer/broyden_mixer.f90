!> Mixer Broyden II modifié (Johnson 1988 / Vanderbilt-Louie).
!>
!> Met à jour itérativement l'inverse du jacobien à partir de
!> l'historique des résidus F_k = x_out_k - x_in_k. À chaque
!> itération k > 1 :
!>
!>   dF_k = F_k - F_{k-1},  dx_k = x_in_k - x_in_{k-1}
!>   u_k  = factor * dF_k + dx_k        (toutes les paires normalisées
!>                                       par |dF_k|)
!>   on résout (w0^2 I + a) c = w * <dF | F_k>
!>   x_new = x_in + factor * F_k - sum_k w_k c_k u_k
!>
!> où a_{ij} = w_i w_j <dF_i | dF_j>. Les poids w_k sont pris constants
!> à 1 (choix usuel quand aucune information de fiabilité par itération
!> n'est disponible). w0 = OMEGA0.
module broyden_mixer
    use kinds,  only: wp
    use mixer,  only: mixer_t
    use errors, only: fatal
    implicit none
    private

    type, extends(mixer_t), public :: broydenMixer_t
        integer  :: history = 6
        real(wp) :: omega0  = 0.01_wp
        integer  :: m       = 0          ! nombre de paires stockées
        real(wp), allocatable :: dF(:,:) ! (n, history) résidus normalisés
        real(wp), allocatable :: u(:,:)  ! (n, history) vecteurs u_k
        real(wp), allocatable :: F_prev(:)
        real(wp), allocatable :: x_prev(:)
    contains
        procedure :: init => broyden_init
        procedure :: mix  => broyden_mix
        procedure :: free => broyden_free
    end type broydenMixer_t

contains

    subroutine broyden_init(self, n, factor, history, omega0)
        class(broydenMixer_t), intent(inout) :: self
        integer,               intent(in)    :: n
        real(wp),              intent(in)    :: factor
        integer,               intent(in)    :: history
        real(wp),              intent(in)    :: omega0
        self%n       = n
        self%factor  = factor
        self%history = max(1, history)
        self%omega0  = omega0
        self%it      = 0
        self%m       = 0
        allocate(self%dF(n, self%history))
        allocate(self%u(n,  self%history))
        allocate(self%F_prev(n))
        allocate(self%x_prev(n))
        self%dF     = 0.0_wp
        self%u      = 0.0_wp
        self%F_prev = 0.0_wp
        self%x_prev = 0.0_wp
    end subroutine broyden_init

    subroutine broyden_mix(self, x_in, x_out, x_new)
        class(broydenMixer_t), intent(inout) :: self
        real(wp),              intent(in)    :: x_in(:)
        real(wp),              intent(in)    :: x_out(:)
        real(wp),              intent(out)   :: x_new(:)

        real(wp) :: F(size(x_in))
        real(wp) :: dF_new(size(x_in)), dx_new(size(x_in))
        real(wp) :: nrm
        real(wp), allocatable :: a(:,:), rhs(:), c(:)
        integer,  allocatable :: ipiv(:)
        integer :: i, j, k, m, info

        F = x_out - x_in

        if (self%it == 0) then
            x_new = x_in + self%factor * F
            self%F_prev = F
            self%x_prev = x_in
            self%it = 1
            return
        end if

        dF_new = F - self%F_prev
        dx_new = x_in - self%x_prev
        nrm    = sqrt(sum(dF_new * dF_new))
        if (nrm <= 0.0_wp) then
            x_new = x_in + self%factor * F
            self%F_prev = F
            self%x_prev = x_in
            self%it = self%it + 1
            return
        end if
        dF_new = dF_new / nrm
        dx_new = dx_new / nrm

        if (self%m < self%history) then
            self%m = self%m + 1
            self%dF(:, self%m) = dF_new
            self%u (:, self%m) = self%factor * dF_new + dx_new
        else
            do k = 1, self%history - 1
                self%dF(:, k) = self%dF(:, k + 1)
                self%u (:, k) = self%u (:, k + 1)
            end do
            self%dF(:, self%history) = dF_new
            self%u (:, self%history) = self%factor * dF_new + dx_new
        end if

        m = self%m
        allocate(a(m, m), rhs(m), c(m), ipiv(m))
        do i = 1, m
            rhs(i) = dot_product(self%dF(:, i), F)
            do j = 1, m
                a(i, j) = dot_product(self%dF(:, i), self%dF(:, j))
            end do
            a(i, i) = a(i, i) + self%omega0 * self%omega0
        end do

        c = rhs
        call dgesv(m, 1, a, m, ipiv, c, m, info)
        if (info /= 0) call fatal("broyden_mixer", "dgesv failed")

        x_new = x_in + self%factor * F
        do k = 1, m
            x_new = x_new - c(k) * self%u(:, k)
        end do

        deallocate(a, rhs, c, ipiv)

        self%F_prev = F
        self%x_prev = x_in
        self%it = self%it + 1
    end subroutine broyden_mix

    subroutine broyden_free(self)
        class(broydenMixer_t), intent(inout) :: self
        if (allocated(self%dF))     deallocate(self%dF)
        if (allocated(self%u))      deallocate(self%u)
        if (allocated(self%F_prev)) deallocate(self%F_prev)
        if (allocated(self%x_prev)) deallocate(self%x_prev)
        self%m  = 0
        self%it = 0
        self%n  = 0
    end subroutine broyden_free
end module broyden_mixer
