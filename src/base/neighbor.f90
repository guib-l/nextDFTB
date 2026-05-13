!> Liste de voisins du système : paires (i < j) telles que la distance
!> r_ij est inférieure ou égale à une limite donnée.
!>
!> Ce module n'importe pas `structure_mod` : il prend en argument les
!> données brutes (natoms, positions) afin d'éviter une dépendance
!> circulaire avec `structure_mod` qui contient un champ de type
!> `neighbor_list_t`.
module neighbor_mod
    use kinds, only: wp
    implicit none
    private

    !> Distance limite par défaut (en bohr) : suffisamment large pour
    !> être traitée comme « infini » dans l'usage actuel sans imposer
    !> de seuil arbitraire.
    real(wp), parameter, public :: DEFAULT_NEIGH_CUTOFF = 100.0_wp

    type, public :: neighbor_list_t
        integer               :: npairs = 0
        integer,  allocatable :: pair_i(:)
        integer,  allocatable :: pair_j(:)
        real(wp), allocatable :: dist(:)
        real(wp), allocatable :: vector(:,:)   ! (3, npairs) : r_j - r_i
    end type neighbor_list_t

    public :: build_neigh_list

contains

    !> Construit la liste des paires (i<j) telles que d(i,j) <= limite_r.
    !> Si `limite_r` est absent, utilise DEFAULT_NEIGH_CUTOFF.
    !> Si `dist` est fourni, l'utilise pour éviter de recalculer les
    !> normes ; sinon les distances sont calculées à la volée.
    subroutine build_neigh_list(natoms, positions, nlist, limite_r, dist)
        integer,               intent(in)  :: natoms
        real(wp),              intent(in)  :: positions(3, natoms)
        type(neighbor_list_t), intent(out) :: nlist
        real(wp), optional,    intent(in)  :: limite_r
        real(wp), optional,    intent(in)  :: dist(natoms, natoms)

        real(wp) :: cutoff, d
        real(wp) :: vec(3)
        integer  :: i, j, k, npairs

        if (present(limite_r)) then
            cutoff = limite_r
        else
            cutoff = DEFAULT_NEIGH_CUTOFF
        end if

        npairs = 0
        do i = 1, natoms
            do j = i + 1, natoms
                if (present(dist)) then
                    d = dist(i, j)
                else
                    d = norm2(positions(:, j) - positions(:, i))
                end if
                if (d <= cutoff) npairs = npairs + 1
            end do
        end do

        nlist%npairs = npairs
        allocate(nlist%pair_i(npairs))
        allocate(nlist%pair_j(npairs))
        allocate(nlist%dist(npairs))
        allocate(nlist%vector(3, npairs))

        k = 0
        do i = 1, natoms
            do j = i + 1, natoms
                vec = positions(:, j) - positions(:, i)
                if (present(dist)) then
                    d = dist(i, j)
                else
                    d = norm2(vec)
                end if
                if (d <= cutoff) then
                    k = k + 1
                    nlist%pair_i(k)    = i
                    nlist%pair_j(k)    = j
                    nlist%dist(k)      = d
                    nlist%vector(:, k) = vec
                end if
            end do
        end do
    end subroutine build_neigh_list

end module neighbor_mod
