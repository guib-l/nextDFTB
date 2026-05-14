!> Énergie répulsive DFTB :
!>
!>   E_rep = Σ_{A<B}  V_rep(elem_A, elem_B; r_AB)
module repulsif
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use skf,           only: skf_get_repulsive => get_repulsive
    implicit none
    private

    real(wp), parameter :: EPS_FD = 1.0e-4_wp

    public :: repulsive_energy, repulsive_grad

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


    !> Composante `dof` (1=x, 2=y, 3=z) de la dérivée de l'énergie
    !> répulsive par rapport à la position de l'atome `katom`,
    !> calculée par différence finie centrée. Seules les paires (k, j)
    !> contribuent : on reconstruit donc uniquement les termes affectés.
    function repulsive_grad(struct, katom, dof) result(de)
        type(structure_t), intent(in) :: struct
        integer,           intent(in) :: katom, dof
        real(wp) :: de

        integer  :: j
        real(wp) :: pos_k(3), pos_kp(3), pos_km(3), rp, rm, ep, em

        pos_k  = struct%atoms(katom)%position
        pos_kp = pos_k
        pos_km = pos_k
        pos_kp(dof) = pos_k(dof) + EPS_FD
        pos_km(dof) = pos_k(dof) - EPS_FD

        ep = 0.0_wp
        em = 0.0_wp
        do j = 1, struct%natoms
            if (j == katom) cycle
            rp = norm2(pos_kp - struct%atoms(j)%position)
            rm = norm2(pos_km - struct%atoms(j)%position)
            ep = ep + skf_get_repulsive(trim(struct%atoms(katom)%symbol), &
                                        trim(struct%atoms(j)%symbol), rp)
            em = em + skf_get_repulsive(trim(struct%atoms(katom)%symbol), &
                                        trim(struct%atoms(j)%symbol), rm)
        end do

        de = (ep - em) / (2.0_wp * EPS_FD)
    end function repulsive_grad

end module repulsif
