!> nextDFTB — Fortran logging bridge to the unified C ABI sink.
!>
!> The actual sink (file handle + formatter + mutex) is implemented in C++
!> (core_cpp/src/logger.cpp). This module exposes four convenience wrappers
!> that Fortran kernels call in-line.
module nextdftb_abi_logging
    use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
    implicit none
    private

    integer(c_int), parameter, public :: NEXTDFTB_LOG_DEBUG = 0
    integer(c_int), parameter, public :: NEXTDFTB_LOG_INFO  = 1
    integer(c_int), parameter, public :: NEXTDFTB_LOG_WARN  = 2
    integer(c_int), parameter, public :: NEXTDFTB_LOG_ERROR = 3

    interface
        subroutine c_nextdftb_log(level, layer, func, msg) &
                bind(C, name="nextdftb_log")
            import :: c_int, c_char
            integer(c_int),         value, intent(in) :: level
            character(kind=c_char),        intent(in) :: layer(*)
            character(kind=c_char),        intent(in) :: func(*)
            character(kind=c_char),        intent(in) :: msg(*)
        end subroutine c_nextdftb_log
    end interface

    public :: log_debug, log_info, log_warn, log_error

contains

    subroutine log_emit(level, func, msg)
        integer(c_int),   intent(in) :: level
        character(len=*), intent(in) :: func
        character(len=*), intent(in) :: msg
        call c_nextdftb_log(level, "fortran"//c_null_char, &
                             trim(func)//c_null_char,      &
                             trim(msg)//c_null_char)
    end subroutine log_emit

    subroutine log_debug(func, msg)
        character(len=*), intent(in) :: func, msg
        call log_emit(NEXTDFTB_LOG_DEBUG, func, msg)
    end subroutine log_debug

    subroutine log_info(func, msg)
        character(len=*), intent(in) :: func, msg
        call log_emit(NEXTDFTB_LOG_INFO, func, msg)
    end subroutine log_info

    subroutine log_warn(func, msg)
        character(len=*), intent(in) :: func, msg
        call log_emit(NEXTDFTB_LOG_WARN, func, msg)
    end subroutine log_warn

    subroutine log_error(func, msg)
        character(len=*), intent(in) :: func, msg
        call log_emit(NEXTDFTB_LOG_ERROR, func, msg)
    end subroutine log_error

end module nextdftb_abi_logging
