!> Helpers d'assertion minimalistes pour les tests unitaires nextDFTB.
!>
!> Aucun framework externe : chaque programme de test appelle ces routines
!> et termine via `test_summary` qui retourne un code de sortie non nul si
!> au moins une assertion a échoué.
module test_utils
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use kinds, only: wp
    implicit none
    private

    integer, save :: n_checks = 0
    integer, save :: n_fail   = 0

    public :: assert_true, assert_eq_int, assert_close, assert_eq_str
    public :: test_summary

contains

    subroutine assert_true(cond, label)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: label
        n_checks = n_checks + 1
        if (.not. cond) then
            n_fail = n_fail + 1
            write(error_unit, '(a,a)') "  [FAIL] ", trim(label)
        end if
    end subroutine assert_true

    subroutine assert_eq_int(actual, expected, label)
        integer,          intent(in) :: actual, expected
        character(len=*), intent(in) :: label
        n_checks = n_checks + 1
        if (actual /= expected) then
            n_fail = n_fail + 1
            write(error_unit, '(a,a,a,i0,a,i0)') "  [FAIL] ", trim(label), &
                " got=", actual, " expected=", expected
        end if
    end subroutine assert_eq_int

    subroutine assert_close(actual, expected, tol, label)
        real(wp),         intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: label
        n_checks = n_checks + 1
        if (abs(actual - expected) > tol) then
            n_fail = n_fail + 1
            write(error_unit, '(a,a,a,es14.6,a,es14.6,a,es10.2)') &
                "  [FAIL] ", trim(label), " got=", actual, &
                " expected=", expected, " tol=", tol
        end if
    end subroutine assert_close

    subroutine assert_eq_str(actual, expected, label)
        character(len=*), intent(in) :: actual, expected, label
        n_checks = n_checks + 1
        if (trim(actual) /= trim(expected)) then
            n_fail = n_fail + 1
            write(error_unit, '(a,a,a,a,a,a,a)') "  [FAIL] ", trim(label), &
                " got='", trim(actual), "' expected='", trim(expected), "'"
        end if
    end subroutine assert_eq_str

    subroutine test_summary(suite_name)
        character(len=*), intent(in) :: suite_name
        write(output_unit, '(a,a,a,i0,a,i0,a)') &
            "[", trim(suite_name), "] ", n_checks - n_fail, "/", n_checks, " passed"
        if (n_fail > 0) error stop 1
    end subroutine test_summary
end module test_utils
