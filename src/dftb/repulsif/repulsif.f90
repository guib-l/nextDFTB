!> Énergie répulsive DFTB :
!>
!>     E_rep = Σ_{A<B}  V_rep(elem_A, elem_B; r_AB)
module repulsif
    use kinds,   only: wp
    use globals, only: geometry_t, basis_system_t
    use readskf, only: skf_store_t
    use interp,  only: vrep_at_r
    implicit none
    private

    public :: repulsive_energy

contains

    function repulsive_energy(geom, bas, store) result(e)
        type(geometry_t),     intent(in) :: geom
        type(basis_system_t), intent(in) :: bas
        type(skf_store_t),    intent(in) :: store
        real(wp) :: e
        integer  :: i, j, ei, ej
        real(wp) :: dx, dy, dz, r

        e = 0.0_wp
        do i = 1, geom%natoms
            do j = i + 1, geom%natoms
                ei = bas%atom_elem(i); ej = bas%atom_elem(j)
                dx = geom%coords(1, j) - geom%coords(1, i)
                dy = geom%coords(2, j) - geom%coords(2, i)
                dz = geom%coords(3, j) - geom%coords(3, i)
                r  = sqrt(dx*dx + dy*dy + dz*dz)
                e = e + vrep_at_r(store%pair(ei, ej), r)
            end do
        end do
    end function repulsive_energy
end module repulsif
