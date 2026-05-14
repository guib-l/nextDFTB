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
    implicit none
    private

    real(wp), parameter :: EPS_FD = 1.0e-4_wp

    public :: build_gamma, select_gamma
    public :: gamma_base, gamma_mean, gamma_stdr
    public :: build_dgamma, select_dgamma
    public :: dgamma_base, dgamma_mean, dgamma_stdr

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
                gamma(i, j) = 1/r - select_gamma(k, bas%elems(ei)%U_s, &
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
    !> de Slater de Hubbard différents).
    function gamma_stdr(U_a, U_b, r) result(g)
        real(wp), intent(in) :: U_a, U_b, r
        real(wp) :: g
        real(wp) :: tauA, tauB, tauM
        real(wp), parameter :: TOL_R = 1.0e-6_wp
        real(wp), parameter :: TOL_U = 1.0e-6_wp

        tauA = 3.2_wp * U_a
        tauB = 3.2_wp * U_b

        if (r < TOL_R) then
            if (abs(U_a - U_b) < TOL_U) then
                g = -0.5_wp * (U_a + U_b)
            else
                g = -0.5_wp * ( (tauA * tauB) / (tauA + tauB) &
                              + (tauA * tauB)**2 * (tauA + tauB)**(-3) )
            end if
        else if (abs(U_a - U_b) < TOL_U) then
            tauM = 0.5_wp * (tauA + tauB)
            g = exp(-tauM * r) * ( 1.0_wp / r &
                + 0.6875_wp * tauM &
                + 0.1875_wp * r * tauM**2 &
                + 0.020833333333_wp * r**2 * tauM**3 )
        else
            g = gamma_E(r, tauA, tauB) + gamma_E(r, tauB, tauA)
        end if
    end function gamma_stdr


    !> Terme exponentiel auxiliaire pour gamma_stdr (Hubbard distincts).
    pure function gamma_E(r, a, b) result(v)
        real(wp), intent(in) :: r, a, b
        real(wp) :: v, d
        d = a*a - b*b
        v = exp(-a * r) * ( (0.5_wp * b**4 * a) / (d**2) - &
                          ( b**6 - 3.0_wp * b**4 * a**2 ) * (r * d**3)**(-1) )
    end function gamma_E


    !-- Dérivées numériques par rapport à r ----------------------------

    !> Dérivée par différence finie de γ_AB = 1/r par rapport à r.
    function dgamma_base(r) result(dg)
        real(wp), intent(in) :: r
        real(wp) :: dg
        dg = (gamma_base(r + EPS_FD) - gamma_base(r - EPS_FD)) &
             / (2.0_wp * EPS_FD)
    end function dgamma_base


    !> Dérivée par différence finie de γ_AB (Klopman-Ohno) par rapport à r.
    function dgamma_mean(U_a, U_b, r) result(dg)
        real(wp), intent(in) :: U_a, U_b, r
        real(wp) :: dg
        dg = (gamma_mean(U_a, U_b, r + EPS_FD) &
              - gamma_mean(U_a, U_b, r - EPS_FD)) / (2.0_wp * EPS_FD)
    end function dgamma_mean


    !> Dérivée par différence finie de γ_AB (Slater) par rapport à r.
    function dgamma_stdr(U_a, U_b, r) result(dg)
        real(wp), intent(in) :: U_a, U_b, r
        real(wp) :: dg
        dg = (gamma_stdr(U_a, U_b, r + EPS_FD) &
              - gamma_stdr(U_a, U_b, r - EPS_FD)) / (2.0_wp * EPS_FD)
    end function dgamma_stdr


    !> Sélecteur de la dérivée des éléments. Retourne 0 par défaut.
    function select_dgamma(kind, U_a, U_b, r) result(dg)
        integer,  intent(in) :: kind
        real(wp), intent(in) :: U_a, U_b, r
        real(wp) :: dg

        select case (kind)
        case (GK_BASE)
            dg = dgamma_base(r)
        case (GK_MEAN)
            dg = dgamma_mean(U_a, U_b, r)
        case (GK_STDR)
            dg = dgamma_stdr(U_a, U_b, r)
        case default
            dg = 0.0_wp
        end select
    end function select_dgamma


    !> Construit la dérivée de la matrice γ (natoms × natoms) par
    !> rapport à la coordonnée `dof` (1=x, 2=y, 3=z) de l'atome `katom`.
    !> Différence finie centrée sur les positions atomiques.
    subroutine build_dgamma(struct, bas, katom, dof, dgamma, kind)
        type(structure_t),    intent(in)  :: struct
        type(basis_system_t), intent(in)  :: bas
        integer,              intent(in)  :: katom, dof
        real(wp),             intent(out) :: dgamma(:,:)
        integer, optional,    intent(in)  :: kind

        integer  :: i, j, ei, ej, k, natoms
        real(wp) :: pos_k(3), pos_kp(3), pos_km(3)
        real(wp) :: rp, rm, gp, gm
        real(wp) :: Ua, Ub

        k = GK_MEAN
        if (present(kind)) k = kind

        natoms = struct%natoms
        dgamma = 0.0_wp

        pos_k  = struct%atoms(katom)%position
        pos_kp = pos_k
        pos_km = pos_k
        pos_kp(dof) = pos_k(dof) + EPS_FD
        pos_km(dof) = pos_k(dof) - EPS_FD

        ! Seules les paires (katom, j) avec j /= katom dépendent de la
        ! position de katom. Le terme diagonal γ_AA = U_A est constant.
        ei = bas%atom_elem(katom)
        Ua = bas%elems(ei)%U_s
        do j = 1, natoms
            if (j == katom) cycle
            ej = bas%atom_elem(j)
            Ub = bas%elems(ej)%U_s
            rp = norm2(pos_kp - struct%atoms(j)%position)
            rm = norm2(pos_km - struct%atoms(j)%position)
            gp = 1.0_wp / rp - select_gamma(k, Ua, Ub, rp)
            gm = 1.0_wp / rm - select_gamma(k, Ua, Ub, rm)
            dgamma(katom, j) = (gp - gm) / (2.0_wp * EPS_FD)
            dgamma(j, katom) = dgamma(katom, j)
        end do

        ! Silence unused: i (boucle non utilisée ici)
        i = 0
    end subroutine build_dgamma

end module gamma_mod
