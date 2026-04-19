// nextDFTB — HPC/MPI: halo exchange (DESIGN ONLY).
//
// Non-blocking pattern (MPI_Isend / MPI_Irecv) with explicit overlap of
// computation and communication:
//
//     start_exchange(plan, field);       // posts all Isend / Irecv
//     // caller runs interior compute here
//     wait_exchange(plan);               // MPI_Waitall on request array
//     // caller runs boundary compute here
//
// Packing and unpacking of halo buffers is a candidate for C++
// reimplementation of the corresponding Fortran kernels (see section 7
// of the architecture document).

#pragma once

#include <cstdint>

#include "nextdftb/hpc/mpi_topology.hpp"

namespace nextdftb::hpc {

struct HaloExchangePlan {
    const CartesianTopology* topology   = nullptr;
    std::int64_t             halo_width = 1;
    void*                    impl       = nullptr;  // opaque: request array,
                                                    // packed buffers, datatypes
};

HaloExchangePlan make_halo_exchange_plan(const CartesianTopology& topo,
                                          std::int64_t             halo_width);

void destroy_halo_exchange_plan(HaloExchangePlan& plan);

// Pack the neighbour-facing slabs of `field_with_halos` and post the
// MPI_Isend / MPI_Irecv calls. Does not block.
void start_exchange(const HaloExchangePlan& plan, double* field_with_halos);

// MPI_Waitall on the pending requests, then unpack received slabs.
void wait_exchange(const HaloExchangePlan& plan);

} // namespace nextdftb::hpc
