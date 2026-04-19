// nextDFTB — HPC/MPI: collective reductions (DESIGN ONLY).
//
// All reductions go through the Cartesian topology's communicator. The
// planned implementation will fuse the local compute with an MPI_Allreduce
// (section 7 of the architecture document — candidate for C++
// reimplementation of a Fortran reduction kernel).

#pragma once

#include <cstdint>

#include "nextdftb/hpc/mpi_topology.hpp"

namespace nextdftb::hpc {

enum class ReductionOp { Sum, Max, Min };

// Reduce a scalar across the topology's communicator. Every rank gets the
// same result.
double allreduce_scalar(const CartesianTopology& topo,
                         double                   local_value,
                         ReductionOp              op);

// Element-wise reduction of a vector.
void allreduce_vector(const CartesianTopology& topo,
                      const double*            local_in,
                      double*                  global_out,
                      std::int64_t             n,
                      ReductionOp              op);

} // namespace nextdftb::hpc
