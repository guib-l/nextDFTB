!> nextDFTB — BIND(C) facade exposing the Fortran layer as C ABI symbols.
!>
!> Every procedure here:
!>   * has BIND(C, name="nextdftb_...") so the linker emits the exact symbol
!>     the public header abi/include/nextdftb/abi.h declares
!>   * accepts only iso_c_binding types (no allocatables, no derived types,
!>     no assumed-shape arrays)
!>   * receives shared buffers as type(c_ptr) + length and rehydrates them
!>     via c_f_pointer
!>   * returns an integer(c_int) status — NEXTDFTB_OK on success, negative
!>     on error (with the error state published via nextdftb_set_error)
module nextdftb_abi
    use, intrinsic :: iso_c_binding
    !$ use omp_lib, only: omp_set_num_threads

    use nextdftb_kinds,     only: wp, ip
    use nextdftb_runtime,   only: is_initialized
    use nextdftb_test_kernel, only: test_compute
    use nextdftb_abi_errors, only: fortran_set_error,                  &
                                    NEXTDFTB_OK,                        &
                                    NEXTDFTB_ERR_INVALID_ARG,           &
                                    NEXTDFTB_ERR_NOT_INITIALIZED,       &
                                    NEXTDFTB_ERR_ALREADY_INITIALIZED,   &
                                    NEXTDFTB_SEV_RECOVERABLE
    use nextdftb_abi_logging, only: log_info
    implicit none
    private

contains

    !----------------------------------------------------------------
    ! Runtime lifecycle
    !----------------------------------------------------------------

    function nextdftb_init(num_threads) result(status) &
            bind(C, name="nextdftb_init")
        integer(c_int), value, intent(in) :: num_threads
        integer(c_int)                    :: status

        if (is_initialized) then
            call fortran_set_error(NEXTDFTB_ERR_ALREADY_INITIALIZED,      &
                                    NEXTDFTB_SEV_RECOVERABLE,              &
                                    "nextdftb_init",                       &
                                    "Fortran runtime already initialized")
            status = NEXTDFTB_ERR_ALREADY_INITIALIZED
            return
        end if

        !$ if (num_threads > 0) call omp_set_num_threads(num_threads)

        is_initialized = .true.
        call log_info("nextdftb_init", "Fortran runtime initialized")
        status = NEXTDFTB_OK
    end function nextdftb_init

    function nextdftb_finalize() result(status) &
            bind(C, name="nextdftb_finalize")
        integer(c_int) :: status

        if (.not. is_initialized) then
            status = NEXTDFTB_OK   ! idempotent
            return
        end if

        is_initialized = .false.
        call log_info("nextdftb_finalize", "Fortran runtime finalized")
        status = NEXTDFTB_OK
    end function nextdftb_finalize

    function nextdftb_is_initialized() result(flag) &
            bind(C, name="nextdftb_is_initialized")
        integer(c_int) :: flag
        if (is_initialized) then
            flag = 1_c_int
        else
            flag = 0_c_int
        end if
    end function nextdftb_is_initialized

    !----------------------------------------------------------------
    ! Infrastructure smoke test
    !----------------------------------------------------------------

    function nextdftb_test(out_ptr) result(status) &
            bind(C, name="nextdftb_test")
        type(c_ptr), value, intent(in) :: out_ptr
        integer(c_int)                 :: status

        real(wp), pointer :: out => null()

        if (.not. is_initialized) then
            call fortran_set_error(NEXTDFTB_ERR_NOT_INITIALIZED,          &
                                    NEXTDFTB_SEV_RECOVERABLE,              &
                                    "nextdftb_test",                       &
                                    "call nextdftb_init() first")
            status = NEXTDFTB_ERR_NOT_INITIALIZED
            return
        end if
        if (.not. c_associated(out_ptr)) then
            call fortran_set_error(NEXTDFTB_ERR_INVALID_ARG,              &
                                    NEXTDFTB_SEV_RECOVERABLE,              &
                                    "nextdftb_test",                       &
                                    "null output pointer")
            status = NEXTDFTB_ERR_INVALID_ARG
            return
        end if

        call c_f_pointer(out_ptr, out)
        call test_compute(out)
        status = NEXTDFTB_OK
    end function nextdftb_test

end module nextdftb_abi
