!> Outils orbitales DFTB : construction du système de base et
!> utilitaires de parsing de la spécification ORBITALS.
module orbitals_mod
    use kinds,         only: wp
    use constants,     only: SYMBOL_LEN
    use structure_mod, only: structure_t
    use property,      only: property_basis_t
    use dftbstate,     only: basis_system_t
    use skf,           only: skf_get_eps         => get_eps,         &
                              skf_get_hubbard     => get_hubbard,     &
                              skf_get_occupations => get_occupations, &
                              skf_nelements       => nelements,       &
                              skf_element_symbol  => element_symbol
    use errors,        only: fatal
    implicit none
    private

    public :: build_basis_system, count_l_shells, lookup_orbital_spec, &
              match_symbol

contains

    !> Construit la table des éléments uniques et le mapping
    !> atome → orbitales en interrogeant l'API publique du module skf.
    !> Le nombre d'orbitales est déterminé par le mot-clé ORBITALS de
    !> l'input (s'il est fourni), sinon par le pattern non-nul des
    !> occupations lues dans les fichiers SKF.
    !> Propage également `size_orb` vers chaque atome de `struct`.
    subroutine build_basis_system(struct, basis_in, bas)
        type(structure_t),    intent(inout) :: struct
        type(property_basis_t),        intent(in)  :: basis_in
        type(basis_system_t), intent(out) :: bas
        integer :: i, k, n_elem, ia, n_s, n_p, n_d, n_f, nls
        real(wp) :: eps(3), hub(3), occ(3)
        real(wp) :: nelec_real
        character(len=SYMBOL_LEN) :: sym
        character(len=64) :: orb_str
        logical :: has_spec

        n_elem = skf_nelements()
        bas%nelem = n_elem
        allocate(bas%elems(n_elem))

        do i = 1, n_elem
            sym = skf_element_symbol(i)
            eps = skf_get_eps(trim(sym))
            hub = skf_get_hubbard(trim(sym))
            occ = skf_get_occupations(trim(sym))

            bas%elems(i)%symbol = sym
            bas%elems(i)%e_d    = eps(1)
            bas%elems(i)%e_p    = eps(2)
            bas%elems(i)%e_s    = eps(3)
            bas%elems(i)%U_s    = hub(3)
            bas%elems(i)%occ_d  = occ(1)
            bas%elems(i)%occ_p  = occ(2)
            bas%elems(i)%occ_s  = occ(3)
            bas%elems(i)%q_neutral = occ(1) + occ(2) + occ(3)

            call lookup_orbital_spec(basis_in, trim(sym), orb_str, has_spec)

            if (has_spec) then
                call count_l_shells(orb_str, n_s, n_p, n_d, n_f)
                if (n_s > 1 .or. n_p > 1 .or. n_d > 1 .or. n_f > 1) &
                    call fatal("orbitals_mod", &
                        "Not implemented yet: multiple shells of same angular momentum in ORBITALS for "//trim(sym))
                if (n_f > 0) &
                    call fatal("orbitals_mod", "f-orbitals not supported yet for "//trim(sym))
                if (n_s + n_p + n_d == 0) &
                    call fatal("orbitals_mod", "no orbitals listed in ORBITALS for "//trim(sym))
                if (n_d == 1) then
                    bas%elems(i)%l_max = 2
                else if (n_p == 1) then
                    bas%elems(i)%l_max = 1
                else
                    bas%elems(i)%l_max = 0
                end if
                bas%elems(i)%n_orb = n_s + 3*n_p + 5*n_d
            else
                ! fallback : déduit du pattern d'occupations SKF
                if (occ(1) > 0.0_wp) then
                    bas%elems(i)%l_max = 2
                else if (occ(2) > 0.0_wp) then
                    bas%elems(i)%l_max = 1
                else
                    bas%elems(i)%l_max = 0
                end if
                select case (bas%elems(i)%l_max)
                case (0); bas%elems(i)%n_orb = 1
                case (1); bas%elems(i)%n_orb = 4
                case (2); bas%elems(i)%n_orb = 9
                end select
            end if

        end do

        allocate(bas%atom_elem(struct%natoms))
        allocate(bas%atom_norb(struct%natoms))
        allocate(bas%atom_orb_start(struct%natoms))
        allocate(bas%atom_nlshell(struct%natoms))
        allocate(bas%atom_lshell_start(struct%natoms))

        bas%norb_total  = 0
        bas%lshell_orbs = 0
        do ia = 1, struct%natoms
            k = match_symbol(bas, struct%atoms(ia)%symbol)
            if (k == 0) call fatal("orbitals_mod", "unknown element: "//trim(struct%atoms(ia)%symbol))
            nls = bas%elems(k)%l_max + 1
            bas%atom_elem(ia)         = k
            bas%atom_norb(ia)         = bas%elems(k)%n_orb
            bas%atom_orb_start(ia)    = bas%norb_total + 1
            bas%norb_total            = bas%norb_total + bas%elems(k)%n_orb
            bas%atom_nlshell(ia)      = nls
            bas%atom_lshell_start(ia) = bas%lshell_orbs + 1
            bas%lshell_orbs           = bas%lshell_orbs + nls
            struct%atoms(ia)%size_orb = bas%elems(k)%n_orb
        end do

        ! Somme en flottant puis arrondi unique : évite que `nint`
        ! par atome ne crée un déséquilibre quand les charges
        ! formelles fractionnaires des atomes se compensent globalement.
        nelec_real = 0.0_wp
        do ia = 1, struct%natoms
            k = bas%atom_elem(ia)
            nelec_real = nelec_real + bas%elems(k)%q_neutral - struct%atoms(ia)%charge
        end do
        bas%nelec = nint(nelec_real)
    end subroutine build_basis_system


    !> Recherche la chaîne ORBITALS associée au symbole `sym` dans
    !> l'input parsé. Retourne `found=.false.` si non spécifiée.
    subroutine lookup_orbital_spec(basis_in, sym, orb_str, found)
        type(property_basis_t),    intent(in)  :: basis_in
        character(len=*), intent(in)  :: sym
        character(len=*), intent(out) :: orb_str
        logical,          intent(out) :: found
        integer :: j
        orb_str = ""
        found   = .false.
        if (.not. allocated(basis_in%orbitals)) return
        do j = 1, size(basis_in%orbitals)
            if (trim(basis_in%orbitals(j)%symbol) == trim(sym)) then
                orb_str = basis_in%orbitals(j)%orbitals
                found   = .true.
                return
            end if
        end do
    end subroutine lookup_orbital_spec


    !> Compte le nombre de shells par moment angulaire dans une chaîne
    !> du type "1s 2s 2p" (insensible à la casse). Les chiffres
    !> (nombre quantique principal) sont ignorés.
    subroutine count_l_shells(orb_str, n_s, n_p, n_d, n_f)
        character(len=*), intent(in)  :: orb_str
        integer,          intent(out) :: n_s, n_p, n_d, n_f
        integer :: i
        character(len=1) :: c
        n_s = 0; n_p = 0; n_d = 0; n_f = 0
        do i = 1, len_trim(orb_str)
            c = orb_str(i:i)
            select case (c)
            case ('s', 'S'); n_s = n_s + 1
            case ('p', 'P'); n_p = n_p + 1
            case ('d', 'D'); n_d = n_d + 1
            case ('f', 'F'); n_f = n_f + 1
            end select
        end do
    end subroutine count_l_shells


    pure function match_symbol(bas, sym) result(idx)
        type(basis_system_t), intent(in) :: bas
        character(len=*),     intent(in) :: sym
        integer :: idx, i
        idx = 0
        do i = 1, bas%nelem
            if (trim(bas%elems(i)%symbol) == trim(sym)) then
                idx = i; return
            end if
        end do
    end function match_symbol

end module orbitals_mod
