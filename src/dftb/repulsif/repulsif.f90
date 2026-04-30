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
        real(wp) :: dx, dy, dz, r

        e = 0.0_wp
        do i = 1, struct%natoms
            do j = i + 1, struct%natoms
                dx = struct%atoms(j)%position(1) - struct%atoms(i)%position(1)
                dy = struct%atoms(j)%position(2) - struct%atoms(i)%position(2)
                dz = struct%atoms(j)%position(3) - struct%atoms(i)%position(3)
                r  = sqrt(dx*dx + dy*dy + dz*dz)
                e = e + skf_get_repulsive(trim(struct%atoms(i)%symbol), &
                                          trim(struct%atoms(j)%symbol), r)
            end do
        end do
    end function repulsive_energy
end module repulsif
