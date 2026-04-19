// nextDFTB — HPC/MPI: Cartesian topology descriptor (DESIGN ONLY).
//
// Implementation is deferred per section 6.2 of the architecture document.
// The goal of this header is to lock in the public interface now so that
// downstream code (orchestrator, kernels targeted for C++ reimplementation)
// can be written against a stable API.
//
// No <mpi.h> include here on purpose: this header is consumable from units
// that are never built with MPI. The concrete MPI_Comm handle lives in the
// implementation translation unit and is opaque to clients.

#pragma once

#include <array>
#include <cstdint>

namespace nextdftb::hpc {

struct CartesianTopology {
    int rank  = 0;
    int size  = 1;
    int ndims = 0;
    std::array<int, 3> dims     {};   // grid extents per dimension
    std::array<int, 3> coords   {};   // this rank's coords inside the grid
    std::array<int, 3> periodic {};   // 0 or 1 per dimension
};

// Build a Cartesian topology from MPI_COMM_WORLD. `ndims` must be 1..3.
// The implementation will call MPI_Dims_create + MPI_Cart_create.
CartesianTopology make_cartesian_topology(int                        ndims,
                                          const std::array<int, 3>&  periodic);

// Tear down the topology and free the associated communicator.
void destroy_cartesian_topology(CartesianTopology& topo);

} // namespace nextdftb::hpc
