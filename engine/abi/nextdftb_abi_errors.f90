!> nextDFTB — Fortran side of the shared error state.
!>
!> The canonical thread-local error state lives in C++ (core_cpp/src/error.cpp).
!> This module mirrors the status/severity codes of abi/include/nextdftb/errors.h
!> and provides a convenience wrapper so Fortran kernels can publish an error
!> by calling the C ABI setter nextdftb_set_error().
module nextdftb_abi_errors
    use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
    implicit none
    private

    ! Status codes — keep in sync with nextdftb_status_t in errors.h
    integer(c_int), parameter, public :: NEXTDFTB_OK                      =   0
    integer(c_int), parameter, public :: NEXTDFTB_ERR_INVALID_ARG         =  -1
    integer(c_int), parameter, public :: NEXTDFTB_ERR_NOT_INITIALIZED     =  -2
    integer(c_int), parameter, public :: NEXTDFTB_ERR_ALREADY_INITIALIZED =  -3
    integer(c_int), parameter, public :: NEXTDFTB_ERR_ALLOCATION          =  -4
    integer(c_int), parameter, public :: NEXTDFTB_ERR_DIMENSION_MISMATCH  =  -5
    integer(c_int), parameter, public :: NEXTDFTB_ERR_NUMERICAL           =  -6
    integer(c_int), parameter, public :: NEXTDFTB_ERR_FATAL               =  -7
    integer(c_int), parameter, public :: NEXTDFTB_ERR_IO                  =  -8

    integer(c_int), parameter, public :: NEXTDFTB_SEV_RECOVERABLE = 0
    integer(c_int), parameter, public :: NEXTDFTB_SEV_FATAL       = 1

    interface
        subroutine c_nextdftb_set_error(code, severity, layer, func, msg) &
                bind(C, name="nextdftb_set_error")
            import :: c_int, c_char
            integer(c_int),         value,  intent(in) :: code
            integer(c_int),         value,  intent(in) :: severity
            character(kind=c_char),         intent(in) :: layer(*)
            character(kind=c_char),         intent(in) :: func(*)
            character(kind=c_char),         intent(in) :: msg(*)
        end subroutine c_nextdftb_set_error
    end interface

    public :: fortran_set_error

contains

    !> Publish an error from the Fortran layer. Layer is hardcoded to "fortran".
    subroutine fortran_set_error(code, severity, func, msg)
        integer(c_int),   intent(in) :: code
        integer(c_int),   intent(in) :: severity
        character(len=*), intent(in) :: func
        character(len=*), intent(in) :: msg

        call c_nextdftb_set_error(code, severity,                       &
                                   "fortran"  // c_null_char,            &
                                   trim(func) // c_null_char,            &
                                   trim(msg)  // c_null_char)
    end subroutine fortran_set_error

end module nextdftb_abi_errors
