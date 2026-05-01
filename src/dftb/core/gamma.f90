!> Matrice gamma DFTB (interaction Coulomb-like entre charges Mulliken).
!>
!> Forme Klopman-Ohno (approximation simple, suffisante pour SCC standard) :
!>   γ_AB(r) = 1 / sqrt(r² + (1/(2·U_avg))²),    U_avg = (U_A+U_B)/2
!>   γ_AA    = U_A
module gamma_mod
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use dftbstate,     only: basis_system_t
    implicit none
    private

    public :: build_gamma

contains

    subroutine build_gamma(struct, bas, gamma)
        type(structure_t),    intent(in)  :: struct
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: gamma(:,:)

        integer  :: i, j, ei, ej
        real(wp) :: r, dx, dy, dz, U_avg, a2

        do i = 1, struct%natoms
            ei = bas%atom_elem(i)
            gamma(i, i) = bas%elems(ei)%U_s
            do j = i + 1, struct%natoms
                ej = bas%atom_elem(j)
                dx = struct%atoms(j)%position(1) - struct%atoms(i)%position(1)
                dy = struct%atoms(j)%position(2) - struct%atoms(i)%position(2)
                dz = struct%atoms(j)%position(3) - struct%atoms(i)%position(3)
                r  = sqrt(dx*dx + dy*dy + dz*dz)
                U_avg = 0.5_wp * (bas%elems(ei)%U_s + bas%elems(ej)%U_s)
                a2 = (1.0_wp / (2.0_wp * U_avg))**2
                gamma(i, j) = 1.0_wp / sqrt(r*r + a2)
                gamma(j, i) = gamma(i, j)
            end do
        end do
    end subroutine build_gamma
end module gamma_mod
