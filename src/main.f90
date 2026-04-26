!> nextDFTB — programme principal.
!>
!> Point d'entrée de l'application. Pour l'instant, squelette qui affiche
!> la version et les constantes ; les phases suivantes ajouteront la
!> lecture d'un fichier d'entrée puis la chaîne SKF → H/S → diagonalisation
!> → SCC → sortie.
program nextdftb
    use constants, only: BOHR_TO_ANGSTROM, HARTREE_TO_EV
    implicit none

    character(len=*), parameter :: version = '0.1.0'

    write(*, '(a)')       '================================================'
    write(*, '(a,1x,a)')  ' nextDFTB', version
    write(*, '(a)')       ' Fortran 2018 — DFTB compute engine'
    write(*, '(a)')       '================================================'
    write(*, '(a)')       ''
    write(*, '(a,f20.10)') ' Bohr  -> Angstrom : ', BOHR_TO_ANGSTROM
    write(*, '(a,f20.10)') ' Ha    -> eV       : ', HARTREE_TO_EV
    write(*, '(a)')       ''
    write(*, '(a)')       ' [phase 0] infra en place. Phases suivantes :'
    write(*, '(a)')       '   1. lecture SKF'
    write(*, '(a)')       '   2. construction H / S'
    write(*, '(a)')       '   3. diagonalisation LAPACK'
    write(*, '(a)')       '   4. boucle SCC (Mulliken)'
    write(*, '(a)')       '   5. CLI + sortie'
end program nextdftb
