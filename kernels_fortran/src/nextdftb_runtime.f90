!> nextDFTB — Fortran runtime state.
!>
!> Holds the global "initialized" flag set by nextdftb_init() and cleared
!> by nextdftb_finalize(). Private working memory of the Fortran layer
!> (kernel scratch buffers, LAPACK workspaces, etc.) is allocated from
!> modules that the kernels themselves own.
module nextdftb_runtime
    implicit none
    private

    logical, save, public :: is_initialized = .false.

end module nextdftb_runtime
