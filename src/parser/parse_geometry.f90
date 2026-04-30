!> Lecture / écriture des fichiers de géométrie au format xyz étendu.
!>
!> Format (`geometry.dat`) :
!>     ligne 1   : nombre d'atomes 'natoms'
!>     ligne 2   : commentaire (ignoré)
!>     ligne 3+  : symbol  x  y  z  charge
!>
!> Coordonnées attendues en angström, converties en bohr en interne.
!> La colonne `charge` est la charge atomique initiale portée par l'atome.
!> Par défaut, tous les atomes appartiennent à la molécule numéro 1.
module parse_geometry
    use kinds,         only: wp
    use units,         only: ang_to_bohr, bohr_to_ang
    use atoms_mod,     only: atoms_t, atoms_set
    use structure_mod, only: structure_t, structure_init, structure_set_atom
    use errors,        only: fatal
    implicit none
    private

    public :: read_geometry, write_geometry

contains

    subroutine read_geometry(filename, struct)
        character(len=*),  intent(in)  :: filename
        type(structure_t), intent(out) :: struct

        integer :: u, ios, natoms, i
        character(len=512) :: line
        character(len=8) :: sym
        real(wp) :: x, y, z, q
        type(atoms_t) :: a

        open(newunit=u, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) call fatal("parse_geometry", "cannot open: "//trim(filename))

        read(u, *, iostat=ios) natoms
        if (ios /= 0 .or. natoms < 1) &
            call fatal("parse_geometry", "bad atom count in "//trim(filename))

        read(u, '(a)', iostat=ios) line     ! commentaire ignoré
        if (ios /= 0) call fatal("parse_geometry", "missing comment line")

        call structure_init(struct, natoms)

        do i = 1, natoms
            q = 0.0_wp
            read(u, *, iostat=ios) sym, x, y, z, q
            if (ios /= 0) then
                read(u, *, iostat=ios) sym, x, y, z
                q = 0.0_wp
                if (ios /= 0) call fatal("parse_geometry", "bad atom line")
            end if
            call atoms_set(a, sym, &
                           [ ang_to_bohr(x), ang_to_bohr(y), ang_to_bohr(z) ], &
                           q, 1)
            call structure_set_atom(struct, i, a)
        end do

        close(u)
    end subroutine read_geometry


    subroutine write_geometry(filename, struct, comment)
        character(len=*),  intent(in) :: filename
        type(structure_t), intent(in) :: struct
        character(len=*),  intent(in) :: comment

        integer :: u, ios, i
        real(wp) :: x, y, z

        open(newunit=u, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) call fatal("parse_geometry", "cannot write: "//trim(filename))

        write(u, '(i0)') struct%natoms
        write(u, '(a)')  trim(comment)
        do i = 1, struct%natoms
            x = bohr_to_ang(struct%atoms(i)%position(1))
            y = bohr_to_ang(struct%atoms(i)%position(2))
            z = bohr_to_ang(struct%atoms(i)%position(3))
            write(u, '(a4, 4(1x,f14.8))') struct%atoms(i)%symbol, x, y, z, &
                                          struct%atoms(i)%charge
        end do
        close(u)
    end subroutine write_geometry
end module parse_geometry
