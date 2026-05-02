!> Énergie répulsive DFTB :
!>
!>   E_rep = Σ_{A<B}  V_rep(elem_A, elem_B; r_AB)
module repulsif
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use skf,           only: skf_get_repulsive => get_repulsive
    implicit none
    private

    public :: repulsive_energy

contains

    function repulsive_energy(struct) result(e)
        type(structure_t), intent(in) :: struct
        real(wp) :: e
        integer  :: i, j

        e = 0.0_wp
        do i = 1, struct%natoms
            do j = i + 1, struct%natoms
                e = e + skf_get_repulsive(trim(struct%atoms(i)%symbol), &
                                          trim(struct%atoms(j)%symbol), &
                                          struct%dist(i, j))
            end do
        end do
    end function repulsive_energy
end module repulsif
