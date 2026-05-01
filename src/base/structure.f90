!> Objet Structure : représentation complète du système moléculaire.
!>
!> Conteneur principal du calcul. Stocke la liste de tous les atomes du
!> système et la liste des molécules (par défaut une seule molécule
!> contenant tous les atomes), ainsi que les grandeurs dérivées :
!> centre de masse global, centres de masse par molécule, matrice des
!> distances atome-atome et matrice des distances molécule-molécule.
module structure_mod
    use kinds,        only: wp
    use atoms_mod,    only: atoms_t
    use molecule_mod, only: molecule_t, molecule_init, molecule_add_atom
    implicit none
    private

    type, public :: structure_t
        integer                       :: natoms      = 0
        type(atoms_t),    allocatable :: atoms(:)            ! (natoms)
        integer                       :: nmolecules  = 0
        type(molecule_t), allocatable :: molecules(:)        ! (nmolecules)
        real(wp)                      :: com(3)      = 0.0_wp
        real(wp),         allocatable :: com_mol(:,:)        ! (3, nmolecules)
        real(wp),         allocatable :: dist(:,:)           ! (natoms, natoms)
        real(wp),         allocatable :: dist_mol(:,:)       ! (nmolecules, nmolecules)
    contains
        procedure :: get_atoms     => structure_get_atoms
        procedure :: get_molecule  => structure_get_molecule
        procedure :: set_com       => structure_set_com
        procedure :: set_com_mol   => structure_set_com_mol
        procedure :: set_dist      => structure_set_dist
        procedure :: get_dist      => structure_get_dist
        procedure :: set_dist_mol  => structure_set_dist_mol
        procedure :: get_dist_mol  => structure_get_dist_mol
    end type structure_t

    public :: structure_init, structure_set_atom

contains

    subroutine structure_init(s, natoms)
        type(structure_t), intent(out) :: s
        integer,           intent(in)  :: natoms
        s%natoms = natoms
        allocate(s%atoms(natoms))
        s%nmolecules = 1
        allocate(s%molecules(1))
        call molecule_init(s%molecules(1), 1)
    end subroutine structure_init

    !> Affecte un atome (1-based) et l'enregistre dans sa molécule.
    subroutine structure_set_atom(s, idx, a)
        type(structure_t), intent(inout) :: s
        integer,           intent(in)    :: idx
        type(atoms_t),     intent(in)    :: a
        integer :: mid
        s%atoms(idx) = a
        mid = a%molecule_id
        if (mid < 1 .or. mid > s%nmolecules) mid = 1
        call molecule_add_atom(s%molecules(mid), idx)
    end subroutine structure_set_atom

    subroutine structure_get_atoms(self, out_atoms)
        class(structure_t), intent(in)  :: self
        type(atoms_t), allocatable, intent(out) :: out_atoms(:)
        out_atoms = self%atoms
    end subroutine structure_get_atoms

    subroutine structure_get_molecule(self, out_molecules)
        class(structure_t), intent(in)  :: self
        type(molecule_t), allocatable, intent(out) :: out_molecules(:)
        out_molecules = self%molecules
    end subroutine structure_get_molecule

    !> Centre de masse du système entier.
    subroutine structure_set_com(self)
        class(structure_t), intent(inout) :: self
        real(wp) :: total_mass
        integer  :: i
        self%com   = 0.0_wp
        total_mass = 0.0_wp
        do i = 1, self%natoms
            self%com   = self%com + self%atoms(i)%mass * self%atoms(i)%position
            total_mass = total_mass + self%atoms(i)%mass
        end do
        if (total_mass > 0.0_wp) self%com = self%com / total_mass
    end subroutine structure_set_com

    !> Centres de masse de chaque molécule, stockés dans `com_mol(:, m)`.
    subroutine structure_set_com_mol(self)
        class(structure_t), intent(inout) :: self
        integer :: m
        if (allocated(self%com_mol)) deallocate(self%com_mol)
        allocate(self%com_mol(3, self%nmolecules))
        do m = 1, self%nmolecules
            call self%molecules(m)%set_com(self%atoms)
            self%com_mol(:, m) = self%molecules(m)%com
        end do
    end subroutine structure_set_com_mol

    !> Matrice des distances atome-atome (symétrique, diagonale nulle).
    subroutine structure_set_dist(self)
        class(structure_t), intent(inout) :: self
        integer :: i, j
        real(wp) :: d
        if (allocated(self%dist)) deallocate(self%dist)
        allocate(self%dist(self%natoms, self%natoms))
        self%dist = 0.0_wp
        do i = 1, self%natoms
            do j = i + 1, self%natoms
                d = norm2(self%atoms(i)%position - self%atoms(j)%position)
                self%dist(i, j) = d
                self%dist(j, i) = d
            end do
        end do
    end subroutine structure_set_dist

    function structure_get_dist(self, index_1, index_2) result(d)
        class(structure_t), intent(in) :: self
        integer,            intent(in) :: index_1, index_2
        real(wp) :: d
        d = self%dist(index_1, index_2)
    end function structure_get_dist

    !> Matrice des distances entre centres de masse de molécules.
    !> Requiert `set_com_mol` au préalable.
    subroutine structure_set_dist_mol(self)
        class(structure_t), intent(inout) :: self
        integer  :: i, j
        real(wp) :: d
        if (allocated(self%dist_mol)) deallocate(self%dist_mol)
        allocate(self%dist_mol(self%nmolecules, self%nmolecules))
        self%dist_mol = 0.0_wp
        do i = 1, self%nmolecules
            do j = i + 1, self%nmolecules
                d = norm2(self%com_mol(:, i) - self%com_mol(:, j))
                self%dist_mol(i, j) = d
                self%dist_mol(j, i) = d
            end do
        end do
    end subroutine structure_set_dist_mol

    function structure_get_dist_mol(self, index_1, index_2) result(d)
        class(structure_t), intent(in) :: self
        integer,            intent(in) :: index_1, index_2
        real(wp) :: d
        d = self%dist_mol(index_1, index_2)
    end function structure_get_dist_mol
end module structure_mod
