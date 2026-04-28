!> Outils d'affichage des matrices H, S, C, density, ...
module write_matrix
    use kinds, only: wp
    implicit none
    private

    public :: print_matrix

contains

    subroutine print_matrix(unit, label, M)
        integer,          intent(in) :: unit
        character(len=*), intent(in) :: label
        real(wp),         intent(in) :: M(:,:)
        integer :: i, j, n, m_

        n  = size(M, 1)
        m_ = size(M, 2)
        write(unit, '(a)') trim(label)
        do i = 1, n
            do j = 1, m_
                write(unit, '(es14.6)', advance='no') M(i, j)
            end do
            write(unit, '(a)') ""
        end do
    end subroutine print_matrix
end module write_matrix
