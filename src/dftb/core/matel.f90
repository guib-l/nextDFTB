!> Construction de l'hamiltonien H et de la matrice d'overlap S DFTB.
!>
!> - blocs diagonaux (i,i) : H_ii = onsite, S_ii = identité
!> - blocs off-diagonaux  : intégrales SK (via API publique skf)
!>                          + rotation cartésienne (skrot)
!>
!> Ordre des orbitales par atome : 1=s, 2=px, 3=py, 4=pz.
!> La version actuelle prend en charge les bases s et s+p (les d-orbitales
!> ne sont pas encore couvertes par la rotation).
module matel
    use kinds,         only: wp
    use constants,     only: SYMBOL_LEN
    use structure_mod, only: structure_t
    use property,      only: property_basis_t
    use dftbstate,     only: basis_system_t, element_info_t
    use skf,           only: skf_get_overlaps    => get_overlaps,    &
                              skf_get_hamiltonian => get_hamiltonian, &
                              skf_get_eps         => get_eps,         &
                              skf_get_hubbard     => get_hubbard,     &
                              skf_get_occupations => get_occupations, &
                              skf_nelements       => nelements,       &
                              skf_element_symbol  => element_symbol
    use errors,        only: fatal
    use skrot,         only: sk_block
    implicit none
    private

    public :: build_basis_system, build_hs, build_hs_diag, build_hs_offdiag

contains

    !> Construit la table des éléments uniques et le mapping
    !> atome → orbitales en interrogeant l'API publique du module skf.
    !> Le nombre d'orbitales est déterminé par le mot-clé ORBITALS de
    !> l'input (s'il est fourni), sinon par le pattern non-nul des
    !> occupations lues dans les fichiers SKF.
    subroutine build_basis_system(struct, basis_in, bas)
        type(structure_t),    intent(in)  :: struct
        type(property_basis_t),        intent(in)  :: basis_in
        type(basis_system_t), intent(out) :: bas
        integer :: i, k, n_elem, ia, n_s, n_p, n_d, n_f
        real(wp) :: eps(3), hub(3), occ(3)
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
                    call fatal("matel", &
                        "Not implemented yet: multiple shells of same angular momentum in ORBITALS for "//trim(sym))
                if (n_f > 0) &
                    call fatal("matel", "f-orbitals not supported yet for "//trim(sym))
                if (n_s + n_p + n_d == 0) &
                    call fatal("matel", "no orbitals listed in ORBITALS for "//trim(sym))
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

            if (bas%elems(i)%l_max == 2) &
                call fatal("matel", "d-orbitals not supported yet for "//trim(sym))
        end do

        allocate(bas%atom_elem(struct%natoms))
        allocate(bas%atom_norb(struct%natoms))
        allocate(bas%atom_orb_start(struct%natoms))

        bas%norb_total = 0
        do ia = 1, struct%natoms
            k = match_symbol(bas, struct%atoms(ia)%symbol)
            if (k == 0) call fatal("matel", "unknown element: "//trim(struct%atoms(ia)%symbol))
            bas%atom_elem(ia)      = k
            bas%atom_norb(ia)      = bas%elems(k)%n_orb
            bas%atom_orb_start(ia) = bas%norb_total + 1
            bas%norb_total         = bas%norb_total + bas%elems(k)%n_orb
        end do

        bas%nelec = 0
        do ia = 1, struct%natoms
            k = bas%atom_elem(ia)
            bas%nelec = bas%nelec + nint(bas%elems(k)%q_neutral - struct%atoms(ia)%charge)
        end do
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


    !> Construit H et S (norb × norb) pour la structure donnée.
    !> Wrapper : remet H et S à zéro puis remplit la diagonale et la
    !> partie off-diagonale via les deux sous-routines dédiées.
    subroutine build_hs(struct, bas, H, S)
        type(structure_t),    intent(in)  :: struct
        type(basis_system_t), intent(in)  :: bas
        real(wp),             intent(out) :: H(:,:), S(:,:)

        H = 0.0_wp
        S = 0.0_wp
        call build_hs_diag(struct, bas, H, S)
        call build_hs_offdiag(struct, bas, H, S)
    end subroutine build_hs


    !> Remplit les blocs diagonaux onsite (i,i) de H et S.
    !> H et S doivent être déjà alloués ; les éléments hors diagonale
    !> ne sont pas modifiés.
    subroutine build_hs_diag(struct, bas, H, S)
        type(structure_t),    intent(in)    :: struct
        type(basis_system_t), intent(in)    :: bas
        real(wp),             intent(inout) :: H(:,:), S(:,:)

        integer :: i, ei, oi, ni

        do i = 1, struct%natoms
            ei = bas%atom_elem(i)
            oi = bas%atom_orb_start(i)
            ni = bas%atom_norb(i)

            S(oi, oi) = 1.0_wp
            H(oi, oi) = bas%elems(ei)%e_s
            if (ni >= 4) then
                S(oi+1, oi+1) = 1.0_wp; H(oi+1, oi+1) = bas%elems(ei)%e_p
                S(oi+2, oi+2) = 1.0_wp; H(oi+2, oi+2) = bas%elems(ei)%e_p
                S(oi+3, oi+3) = 1.0_wp; H(oi+3, oi+3) = bas%elems(ei)%e_p
            end if
        end do
    end subroutine build_hs_diag


    !> Remplit les blocs off-diagonaux (i ≠ j) de H et S.
    !> Exploite la symétrie : seul le triangle supérieur (j > i) est
    !> calculé, le triangle inférieur est obtenu par miroir.
    subroutine build_hs_offdiag(struct, bas, H, S)
        type(structure_t),    intent(in)    :: struct
        type(basis_system_t), intent(in)    :: bas
        real(wp),             intent(inout) :: H(:,:), S(:,:)

        integer :: i, j, oi, oj, ni, nj, mu, nu
        real(wp) :: r, dir(3)
        real(wp) :: hblk(4,4), sblk(4,4)
        real(wp), allocatable :: h_ab(:), s_ab(:), h_ba(:), s_ba(:)
        character(len=SYMBOL_LEN) :: sym_i, sym_j

        do i = 1, struct%natoms
            sym_i = struct%atoms(i)%symbol
            do j = i + 1, struct%natoms
                sym_j = struct%atoms(j)%symbol
                oi = bas%atom_orb_start(i); oj = bas%atom_orb_start(j)
                ni = bas%atom_norb(i);      nj = bas%atom_norb(j)

                r = struct%dist(i, j)
                if (r < 1.0e-8_wp) cycle
                dir = (struct%atoms(j)%position - struct%atoms(i)%position) / r

                h_ab = skf_get_hamiltonian(trim(sym_i), trim(sym_j), r)
                s_ab = skf_get_overlaps   (trim(sym_i), trim(sym_j), r)
                h_ba = skf_get_hamiltonian(trim(sym_j), trim(sym_i), r)
                s_ba = skf_get_overlaps   (trim(sym_j), trim(sym_i), r)

                call sk_block(h_ab, s_ab, h_ba, s_ba, dir, ni, nj, hblk, sblk)

                do mu = 1, ni
                    do nu = 1, nj
                        H(oi+mu-1, oj+nu-1) = hblk(mu, nu)
                        S(oi+mu-1, oj+nu-1) = sblk(mu, nu)
                        H(oj+nu-1, oi+mu-1) = hblk(mu, nu)
                        S(oj+nu-1, oi+mu-1) = sblk(mu, nu)
                    end do
                end do
            end do
        end do

        if (allocated(h_ab)) deallocate(h_ab)
        if (allocated(s_ab)) deallocate(s_ab)
        if (allocated(h_ba)) deallocate(h_ba)
        if (allocated(s_ba)) deallocate(s_ba)
    end subroutine build_hs_offdiag
end module matel
