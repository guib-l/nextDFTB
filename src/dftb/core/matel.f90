!> Construction de l'hamiltonien H et de la matrice d'overlap S DFTB.
!>
!> - blocs diagonaux (i,i) : H_ii = onsite, S_ii = identité
!> - blocs off-diagonaux  : intégrales SK (via API publique skf)
!>                          + rotation cartésienne (skrot)
!>
!> Ordre des orbitales par atome :
!>   1=s, 2=px, 3=py, 4=pz,
!>   5=dxy, 6=dyz, 7=dzx, 8=dx2−y2, 9=d3z2−r2
module matel
    use kinds,         only: wp
    use constants,     only: SYMBOL_LEN
    use structure_mod, only: structure_t
    use dftbstate,     only: basis_system_t
    use skf,           only: skf_get_overlaps    => get_overlaps,    &
                              skf_get_hamiltonian => get_hamiltonian
    use skrot,         only: transform_sk
    implicit none
    private

    real(wp), parameter :: EPS_FD = 1.0e-4_wp

    public :: build_hs, build_hs_diag, build_hs_offdiag
    public :: build_off_block, build_diag_block
    public :: build_off_dblock, build_diag_dblock

contains

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
    subroutine build_hs_diag(struct, bas, H, S)
        type(structure_t),    intent(in)    :: struct
        type(basis_system_t), intent(in)    :: bas
        real(wp),             intent(inout) :: H(:,:), S(:,:)

        integer :: i, oi, ni, mu
        real(wp), allocatable :: hblk(:,:), sblk(:,:)

        do i = 1, struct%natoms
            ni = bas%atom_norb(i)
            oi = bas%atom_orb_start(i)
            allocate(hblk(ni, ni), sblk(ni, ni))
            call build_diag_block(bas, i, hblk, sblk)
            do mu = 1, ni
                H(oi+mu-1, oi+mu-1) = hblk(mu, mu)
                S(oi+mu-1, oi+mu-1) = sblk(mu, mu)
            end do
            deallocate(hblk, sblk)
        end do
    end subroutine build_hs_diag


    !> Remplit les blocs off-diagonaux (i ≠ j) de H et S.
    subroutine build_hs_offdiag(struct, bas, H, S)
        type(structure_t),    intent(in)    :: struct
        type(basis_system_t), intent(in)    :: bas
        real(wp),             intent(inout) :: H(:,:), S(:,:)

        integer :: i, j, oi, oj, ni, nj, mu, nu
        real(wp), allocatable :: hblk(:,:), sblk(:,:)

        do i = 1, struct%natoms
            do j = 1 + i, struct%natoms
                ni = bas%atom_norb(i); nj = bas%atom_norb(j)
                oi = bas%atom_orb_start(i); oj = bas%atom_orb_start(j)
                allocate(hblk(ni, nj), sblk(ni, nj))
                call build_off_block(struct, bas, i, j, hblk, sblk)
                do mu = 1, ni
                    do nu = 1, nj
                        H(oi+mu-1, oj+nu-1) = hblk(mu, nu)
                        S(oi+mu-1, oj+nu-1) = sblk(mu, nu)
                        H(oj+nu-1, oi+mu-1) = hblk(mu, nu)
                        S(oj+nu-1, oi+mu-1) = sblk(mu, nu)
                    end do
                end do
                deallocate(hblk, sblk)
            end do
        end do
    end subroutine build_hs_offdiag


    !> Calcule le bloc diagonal (ni × ni) pour l'atome i :
    !> S = identité, H = diag(eps_s, eps_p, eps_p, eps_p, eps_d×5).
    subroutine build_diag_block(bas, i, hblk, sblk)
        type(basis_system_t), intent(in)  :: bas
        integer,              intent(in)  :: i
        real(wp),             intent(out) :: hblk(:,:), sblk(:,:)
        integer :: ei, ni, k

        hblk = 0.0_wp; sblk = 0.0_wp
        ei = bas%atom_elem(i)
        ni = bas%atom_norb(i)

        sblk(1, 1) = 1.0_wp
        hblk(1, 1) = bas%elems(ei)%e_s
        if (ni >= 4) then
            do k = 2, 4
                sblk(k, k) = 1.0_wp
                hblk(k, k) = bas%elems(ei)%e_p
            end do
        end if
        if (ni >= 9) then
            do k = 5, 9
                sblk(k, k) = 1.0_wp
                hblk(k, k) = bas%elems(ei)%e_d
            end do
        end if
    end subroutine build_diag_block


    !> Calcule le bloc hors-diagonal (ni × nj) entre les atomes i et j
    !> via skf + skrot. Boucle élément par élément en appelant
    !> `transform_sk` pour chaque (mu, nu).
    subroutine build_off_block(struct, bas, i, j, hblk, sblk)
        type(structure_t),    intent(in)  :: struct
        type(basis_system_t), intent(in)  :: bas
        integer,              intent(in)  :: i, j
        real(wp),             intent(out) :: hblk(:,:), sblk(:,:)

        call off_block_at(struct%atoms(i)%position, struct%atoms(j)%position, &
                          struct%atoms(i)%symbol, struct%atoms(j)%symbol,     &
                          bas%atom_norb(i), bas%atom_norb(j), hblk, sblk)
    end subroutine build_off_block


    !> Dérivée numérique centrée du bloc hors-diagonal (ni × nj) entre
    !> les atomes i et j par rapport à la coordonnée `dof` (1=x, 2=y,
    !> 3=z) de l'atome `katom`. Si katom n'est ni i ni j la dérivée est
    !> nulle. Pas de différence finie : EPS_FD (paramètre de module).
    subroutine build_off_dblock(struct, bas, i, j, katom, dof, dhblk, dsblk)
        type(structure_t),    intent(in)  :: struct
        type(basis_system_t), intent(in)  :: bas
        integer,              intent(in)  :: i, j, katom, dof
        real(wp),             intent(out) :: dhblk(:,:), dsblk(:,:)

        real(wp) :: pos_i(3), pos_j(3), pos_ip(3), pos_im(3)
        real(wp) :: pos_jp(3), pos_jm(3)
        real(wp), allocatable :: h_p(:,:), s_p(:,:), h_m(:,:), s_m(:,:)
        integer  :: ni, nj

        dhblk = 0.0_wp
        dsblk = 0.0_wp
        if (katom /= i .and. katom /= j) return

        ni = bas%atom_norb(i)
        nj = bas%atom_norb(j)
        allocate(h_p(ni, nj), s_p(ni, nj), h_m(ni, nj), s_m(ni, nj))

        pos_i = struct%atoms(i)%position
        pos_j = struct%atoms(j)%position
        pos_ip = pos_i
        pos_im = pos_i
        pos_jp = pos_j
        pos_jm = pos_j

        if (katom == i) then
            pos_ip(dof) = pos_i(dof) + EPS_FD
            pos_im(dof) = pos_i(dof) - EPS_FD
        else
            pos_jp(dof) = pos_j(dof) + EPS_FD
            pos_jm(dof) = pos_j(dof) - EPS_FD
        end if

        call off_block_at(pos_ip, pos_jp, struct%atoms(i)%symbol,            &
                          struct%atoms(j)%symbol, ni, nj, h_p, s_p)
        call off_block_at(pos_im, pos_jm, struct%atoms(i)%symbol,            &
                          struct%atoms(j)%symbol, ni, nj, h_m, s_m)

        dhblk = (h_p - h_m) / (2.0_wp * EPS_FD)
        dsblk = (s_p - s_m) / (2.0_wp * EPS_FD)

        deallocate(h_p, s_p, h_m, s_m)
    end subroutine build_off_dblock


    !> Dérivée du bloc diagonal (ni × ni). Les onsites H_ii et S_ii ne
    !> dépendent pas des positions atomiques : la dérivée est nulle.
    subroutine build_diag_dblock(bas, i, katom, dof, dhblk, dsblk)
        type(basis_system_t), intent(in)  :: bas
        integer,              intent(in)  :: i, katom, dof
        real(wp),             intent(out) :: dhblk(:,:), dsblk(:,:)
        integer :: dummy
        dummy = bas%atom_norb(i) + katom + dof  ! silence unused warnings
        dhblk = 0.0_wp
        dsblk = 0.0_wp
    end subroutine build_diag_dblock


    !-- Helpers internes -----------------------------------------------

    !> Construit un bloc off-diagonal (ni × nj) entre deux atomes
    !> caractérisés par leurs positions et symboles. Implémente la
    !> logique partagée par `build_off_block` et la différence finie.
    subroutine off_block_at(pos_i, pos_j, sym_i, sym_j, ni, nj, hblk, sblk)
        real(wp),         intent(in)  :: pos_i(3), pos_j(3)
        character(len=*), intent(in)  :: sym_i, sym_j
        integer,          intent(in)  :: ni, nj
        real(wp),         intent(out) :: hblk(:,:), sblk(:,:)

        integer  :: mu, nu
        real(wp) :: r, dir(3), x, y, z
        real(wp), allocatable :: h_ab(:), s_ab(:), h_ba(:), s_ba(:)

        hblk = 0.0_wp; sblk = 0.0_wp

        r = norm2(pos_j - pos_i)
        if (r < 1.0e-8_wp) return
        dir = (pos_j - pos_i) / r
        x = dir(1); y = dir(2); z = dir(3)

        h_ab = skf_get_hamiltonian(trim(sym_i), trim(sym_j), r)
        s_ab = skf_get_overlaps   (trim(sym_i), trim(sym_j), r)
        h_ba = skf_get_hamiltonian(trim(sym_j), trim(sym_i), r)
        s_ba = skf_get_overlaps   (trim(sym_j), trim(sym_i), r)

        do mu = 1, ni
            do nu = 1, nj
                hblk(mu, nu) = transform_sk(x, y, z, mu, nu, h_ab, h_ba)
                sblk(mu, nu) = transform_sk(x, y, z, mu, nu, s_ab, s_ba)
            end do
        end do

        deallocate(h_ab, s_ab, h_ba, s_ba)
    end subroutine off_block_at

end module matel
