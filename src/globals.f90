!> Globals : module réservé aux constantes globales véritablement partagées.
!>
!> Le plan interdit les variables globales implicites. Ce module est
!> volontairement réduit à des constantes nominales. Les structures
!> d'input (basis, calcul, sortie) sont définies dans `parser/`.
module globals
    implicit none
    private
end module globals
