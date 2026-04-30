!> Mixer de charges : simple (linéaire) et historique de type Anderson.
module mixer
    use kinds,    only: wp
    use defaults, only: DEFAULT_MIX_ALPHA, DEFAULT_BROYDEN_HISTORY
    implicit none
    private

    type, public :: mixer_t
        integer  :: kind     = 0           ! 0 = simple, 1 = anderson/broyden
        real(wp) :: alpha    = 0.2_wp
        integer  :: hist_max = 6
        integer  :: it       = 0
        real(wp), allocatable :: dF(:,:)
        real(wp), allocatable :: dx(:,:)
        real(wp), allocatable :: F_prev(:)
        real(wp), allocatable :: x_prev(:)
    end type mixer_t

    public :: init_mixer, mix_charges, free_mixer

contains

    subroutine init_mixer(mx, n, kind)
        type(mixer_t), intent(out) :: mx
        integer,       intent(in)  :: n, kind
        mx%kind     = kind
        mx%alpha    = DEFAULT_MIX_ALPHA
        mx%hist_max = DEFAULT_BROYDEN_HISTORY
        mx%it = 0
        if (kind == 1) then
            allocate(mx%dF(n, mx%hist_max), mx%dx(n, mx%hist_max))
            allocate(mx%F_prev(n), mx%x_prev(n))
            mx%dF = 0.0_wp; mx%dx = 0.0_wp
            mx%F_prev = 0.0_wp; mx%x_prev = 0.0_wp
        end if
    end subroutine init_mixer

    subroutine free_mixer(mx)
        type(mixer_t), intent(inout) :: mx
        if (allocated(mx%dF))     deallocate(mx%dF)
        if (allocated(mx%dx))     deallocate(mx%dx)
        if (allocated(mx%F_prev)) deallocate(mx%F_prev)
        if (allocated(mx%x_prev)) deallocate(mx%x_prev)
    end subroutine free_mixer

    !> Mise à jour de x_in vers x_new à partir du résidu F = x_out - x_in.
    !>   simple   : x_new = x_in + alpha * F
    !>   anderson : mixing linéaire à un pas avec historique stocké.
    subroutine mix_charges(mx, x_in, x_out, x_new)
        type(mixer_t), intent(inout) :: mx
        real(wp),      intent(in)    :: x_in(:), x_out(:)
        real(wp),      intent(out)   :: x_new(:)

        real(wp) :: F(size(x_in))
        integer  :: i

        F = x_out - x_in

        if (mx%kind == 0) then
            x_new = x_in + mx%alpha * F
            mx%it = mx%it + 1
            return
        end if

        if (mx%it == 0) then
            x_new = x_in + mx%alpha * F
            mx%F_prev = F
            mx%x_prev = x_in
            mx%it = 1
            return
        end if

        do i = mx%hist_max, 2, -1
            mx%dF(:, i) = mx%dF(:, i - 1)
            mx%dx(:, i) = mx%dx(:, i - 1)
        end do
        mx%dF(:, 1) = F    - mx%F_prev
        mx%dx(:, 1) = x_in - mx%x_prev

        x_new = x_in + mx%alpha * F

        mx%F_prev = F
        mx%x_prev = x_in
        mx%it = mx%it + 1
    end subroutine mix_charges
end module mixer
