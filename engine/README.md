# `engine/` — Fortran 2018 compute layer

This directory holds the Fortran side of nextDFTB. It is split by
responsibility, not by file type:

```
engine/
├── core/   DFTB program: kinds, constants, runtime state, numerical
│           kernels, and the optional standalone Fortran driver. Pure
│           Fortran modules with no iso_c_binding exposure.
└── abi/    BIND(C) facade: translates the kernels in core/ into the
            C ABI symbols declared in abi/include/nextdftb/. No physics,
            no heavy computation — only iso_c_binding marshalling,
            argument validation, and error publication.
```

The dependency direction is strict: `abi/` depends on `core/`; `core/`
never depends on `abi/`. A change to the numerical code should not
require touching the ABI, and a change to the exported C symbols should
not require touching the kernels.

---

## Adding a new feature

### 1. Write the kernel in `core/` (no C ABI yet)

Add a new module in `engine/core/`, e.g. `nextdftb_saxpy.f90`:

```fortran
module nextdftb_saxpy_kernel
    use nextdftb_kinds, only: wp, ip
    implicit none
    private

    public :: saxpy_compute

contains

    subroutine saxpy_compute(n, a, x, y)
        integer(ip), intent(in)    :: n
        real(wp),    intent(in)    :: a
        real(wp),    intent(in)    :: x(n)
        real(wp),    intent(inout) :: y(n)
        integer(ip) :: i

        !$omp parallel do default(none) shared(a, x, y, n) private(i) schedule(static)
        do i = 1_ip, n
            y(i) = a * x(i) + y(i)
        end do
        !$omp end parallel do
    end subroutine saxpy_compute

end module nextdftb_saxpy_kernel
```

Rules for `core/` modules:

- Work in `wp` (real64) and `ip` (int64) from `nextdftb_kinds`.
- Use assumed-shape / explicit-shape Fortran arrays — *no `c_ptr`*.
- Parallelize with OpenMP when it pays off. Declare data scope
  explicitly (`default(none)`).
- Use constants from `nextdftb_constants` rather than magic numbers.
- Do not call `c_nextdftb_set_error` or `c_nextdftb_log` directly —
  leave error publication and logging to `abi/`.

### 2. Declare the C symbol in `abi/include/nextdftb/abi.h`

Pick a stable name with the `nextdftb_` prefix and an `int` status
return. Pointers go through `c_ptr` + length on the Fortran side.

### 3. Wire the BIND(C) facade in `abi/`

Either extend `engine/abi/nextdftb_abi.f90` or add a sibling module.
The facade must:

- Validate every pointer (`c_associated`) and length.
- Call `fortran_set_error(NEXTDFTB_ERR_*, …)` on any validation
  failure and return the matching status.
- Rehydrate buffers with `c_f_pointer(ptr, array, [n])`.
- Delegate to the kernel from `core/`. It should not contain a loop.

```fortran
function nextdftb_saxpy(n, a, x_ptr, y_ptr) result(status) &
        bind(C, name="nextdftb_saxpy")
    use nextdftb_kinds,         only: wp, ip
    use nextdftb_runtime,       only: is_initialized
    use nextdftb_saxpy_kernel,  only: saxpy_compute
    use nextdftb_abi_errors,    only: fortran_set_error, NEXTDFTB_OK, &
                                      NEXTDFTB_ERR_INVALID_ARG,        &
                                      NEXTDFTB_ERR_NOT_INITIALIZED,    &
                                      NEXTDFTB_SEV_RECOVERABLE
    integer(c_int64_t), value, intent(in) :: n
    real(c_double),     value, intent(in) :: a
    type(c_ptr),        value, intent(in) :: x_ptr, y_ptr
    integer(c_int)                        :: status

    real(wp), pointer :: x(:), y(:)

    if (.not. is_initialized) then
        call fortran_set_error(NEXTDFTB_ERR_NOT_INITIALIZED,     &
                               NEXTDFTB_SEV_RECOVERABLE,          &
                               "nextdftb_saxpy",                  &
                               "call nextdftb_init() first")
        status = NEXTDFTB_ERR_NOT_INITIALIZED
        return
    end if
    if (n <= 0 .or. .not. c_associated(x_ptr) .or. .not. c_associated(y_ptr)) then
        call fortran_set_error(NEXTDFTB_ERR_INVALID_ARG,          &
                               NEXTDFTB_SEV_RECOVERABLE,          &
                               "nextdftb_saxpy",                  &
                               "invalid size or null pointer")
        status = NEXTDFTB_ERR_INVALID_ARG
        return
    end if

    call c_f_pointer(x_ptr, x, [n])
    call c_f_pointer(y_ptr, y, [n])
    call saxpy_compute(int(n, ip), a, x, y)
    status = NEXTDFTB_OK
end function nextdftb_saxpy
```

### 4. Register the source in `engine/CMakeLists.txt`

Add the new file to either `NEXTDFTB_FORTRAN_CORE_SOURCES` (kernel) or
`NEXTDFTB_FORTRAN_ABI_SOURCES` (facade). Leave the build order to
CMake's Fortran module dependency tracking.

### 5. Consume the symbol upstream

Once the C symbol exists, add the C++ wrapper in `core_cpp/`, expose it
via pybind11 in `bindings/python/abi/module.cpp`, and add a Python
wrapper in `bindings/python/api/nextdftb/`.

---

## Error and status codes

All ABI functions must return an `integer(c_int)` status. The codes are
defined in one place (`abi/include/nextdftb/errors.h`) and mirrored for
Fortran in `engine/abi/nextdftb_abi_errors.f90`. Always publish the
details through `fortran_set_error` before returning the status — the
C++ layer rethrows the thread-local error as a typed exception.

## Logging

From inside the engine, use `log_debug/info/warn/error` from
`nextdftb_abi_logging`. These route to the same file sink as the C++
and Python layers, so a single log captures events from every layer.

## Precision and kinds

Never hard-code `real(8)` or `integer(4)`. Always go through
`nextdftb_kinds` (`wp`, `ip`, `i4`). The `c_double`, `c_int`, and
`c_int64_t` constants are re-exported from the same module.

## Tests

Unit tests live under `tests/unit/`. Add a smoke test for each new
ABI entry point that exercises:

1. happy path,
2. `NEXTDFTB_ERR_NOT_INITIALIZED` (call before `nextdftb_init`),
3. `NEXTDFTB_ERR_INVALID_ARG` (null pointer, bad size).
