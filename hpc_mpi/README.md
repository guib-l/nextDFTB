# `hpc_mpi/` — distributed computing layer

**Status: DESIGN ONLY. Implementation deferred.**

This directory freezes the public interface of the planned MPI layer
(section 6.2 of the architecture document) so that orchestrator code can
be written against it today. No MPI call is made from anywhere in the
current build: when `NEXTDFTB_ENABLE_MPI=OFF` (default) this directory
contributes only headers; when `NEXTDFTB_ENABLE_MPI=ON` it builds stubs
that throw `std::logic_error`.

## Interfaces

| Header                                    | Responsibility                                  |
|-------------------------------------------|-------------------------------------------------|
| `include/nextdftb/hpc/mpi_topology.hpp`   | Cartesian topology (1D/2D/3D) + `MPI_Cart_create` |
| `include/nextdftb/hpc/halo_exchange.hpp`  | Non-blocking halo exchange with compute/comm overlap |
| `include/nextdftb/hpc/reductions.hpp`     | `MPI_Allreduce` wrappers for scalar and vector   |

## Planned decisions (to be locked in at implementation time)

- MPI thread-safety: `MPI_THREAD_FUNNELED` vs `MPI_THREAD_MULTIPLE` —
  deferred.
- Overlap strategy: `MPI_Isend`/`MPI_Irecv` + interior-first compute,
  then `MPI_Waitall` + boundary compute.
- Halo packing and reductions are candidates for C++ reimplementation of
  their Fortran counterparts (section 7) because the ownership of the
  MPI state lives in C++.
