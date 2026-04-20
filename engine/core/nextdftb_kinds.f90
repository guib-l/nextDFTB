!> nextDFTB — kind parameters (Fortran 2018)
!>
!> Single source of truth for numerical precision and index width.
!> Working precision is real64 (IEEE double). Indices are int64 to match
!> the c_int64_t used across the ABI.
module nextdftb_kinds
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: iso_c_binding,   only: c_double, c_int, c_int64_t
    implicit none
    private

    integer, parameter, public :: wp = real64    !< working real precision
    integer, parameter, public :: ip = int64     !< internal index kind
    integer, parameter, public :: i4 = int32     !< small int (status codes, log levels)

    public :: c_double, c_int, c_int64_t
end module nextdftb_kinds
