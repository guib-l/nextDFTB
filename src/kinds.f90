!> nextDFTB — kind parameters.
!>
!> Single source of truth for numerical precision and index width.
module kinds
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    implicit none
    private

    integer, parameter, public :: wp = real64    !< working real precision
    integer, parameter, public :: ip = int64     !< internal index kind
    integer, parameter, public :: i4 = int32     !< small int (status codes)
end module kinds
