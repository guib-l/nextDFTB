!> Matrice gamma DFTB (interaction Coulomb-like entre charges Mulliken).
!>
!> Trois variantes de calcul exposées via un sélecteur :
!>   - GK_BASE : point-charges purs, γ_AB = 1/R_AB
!>   - GK_MEAN : forme Klopman-Ohno avec U_avg = (U_A+U_B)/2
!>                  γ_AB = 1 / sqrt(R² + (1/(2 U_avg))²)
!>   - GK_STDR : non implémenté pour le moment
!>
!> Le terme diagonal γ_AA = U_A est conservé pour les trois variantes
!> (convention DFTB standard pour le self-Coulomb atomique).
module gamma_mod
    use kinds,         only: wp
    use structure_mod, only: structure_t
    use dftbstate,     only: basis_system_t
    use property,      only: GK_BASE, GK_MEAN, GK_STDR
    use errors,        only: fatal
    implicit none
    private

    public :: build_gamma, select_gamma
    public :: gamma_base, gamma_mean, gamma_stdr

contains

    !> Construit la matrice gamma (natoms × natoms).
    !> `kind` est optionnel ; à défaut, le sélecteur GK_MEAN est utilisé.
    subroutine build_gamma(struct, bas, gamma, kind)
        type(structure_t),    intent(in)  :: struct
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: gamma(:,:)
        integer, optional,    intent(in)  :: kind

        integer  :: i, j, ei, ej, k
        real(wp) :: r

        k = GK_MEAN
        if (present(kind)) k = kind

        do i = 1, struct%natoms
            ei = bas%atom_elem(i)
            gamma(i, i) = bas%elems(ei)%U_s
            do j = i + 1, struct%natoms
                ej = bas%atom_elem(j)
                r  = struct%dist(i, j)
                gamma(i, j) = select_gamma(k, bas%elems(ei)%U_s, &
                                              bas%elems(ej)%U_s, r)
                gamma(j, i) = gamma(i, j)
            end do
        end do
    end subroutine build_gamma


    !> Sélecteur : dispatche le calcul de γ_AB selon la variante choisie.
    function select_gamma(kind, U_a, U_b, r) result(g)
        integer,  intent(in) :: kind
        real(wp), intent(in) :: U_a, U_b, r
        real(wp) :: g

        select case (kind)
        case (GK_BASE)
            g = gamma_base(r)
        case (GK_MEAN)
            g = gamma_mean(U_a, U_b, r)
        case (GK_STDR)
            g = gamma_stdr(U_a, U_b, r)
        case default
            g = gamma_mean(U_a, U_b, r)
        end select
    end function select_gamma


    !> γ_AB = 1 / R_AB (point-charges purs).
    pure function gamma_base(r) result(g)
        real(wp), intent(in) :: r
        real(wp) :: g
        g = 1.0_wp / r
    end function gamma_base


    !> γ_AB Klopman-Ohno : 1 / sqrt(R² + (1/(2 U_avg))²),
    !> avec U_avg = (U_A + U_B) / 2.
    pure function gamma_mean(U_a, U_b, r) result(g)
        real(wp), intent(in) :: U_a, U_b, r
        real(wp) :: g, U_avg, a2
        U_avg = 0.5_wp * (U_a + U_b)
        a2    = (1.0_wp / (2.0_wp * U_avg))**2
        g     = 1.0_wp / sqrt(r*r + a2)
    end function gamma_mean


    !> γ_AB selon le schéma standard DFTB+ (analytique pour fonctions
    !> de Slater de Hubbard différents). Non implémenté pour le moment.
    function gamma_stdr(U_a, U_b, r) result(g)
        real(wp), intent(in) :: U_a, U_b, r
        real(wp) :: g
        g = 0.0_wp * (U_a + U_b + r)
        call fatal("gamma_mod", "gamma_stdr not implemented yet")
    end function gamma_stdr
end module gamma_mod
