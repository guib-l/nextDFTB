!> Précision numérique unifiée du code nextDFTB.
module kinds
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    implicit none
    private

    integer, parameter, public :: wp = real64
    integer, parameter, public :: ip = int64
    integer, parameter, public :: i4 = int32
end module kinds
