// nextDFTB — HPC/MPI stubs.
//
// Compiled only when NEXTDFTB_ENABLE_MPI=ON. Every entry point throws
// std::logic_error("not implemented") until the real MPI implementation
// lands. The point of shipping these is to make linker errors obvious if
// any code accidentally calls into the HPC layer before it is ready.

#include <stdexcept>
#include <string>

#include "nextdftb/hpc/halo_exchange.hpp"
#include "nextdftb/hpc/mpi_topology.hpp"
#include "nextdftb/hpc/reductions.hpp"

namespace nextdftb::hpc {

namespace {
[[noreturn]] void not_implemented(const char* who) {
    throw std::logic_error(std::string("nextdftb::hpc: ") + who
                           + " is not implemented yet "
                             "(HPC/MPI layer deferred)");
}
} // namespace

CartesianTopology make_cartesian_topology(int, const std::array<int, 3>&) {
    not_implemented("make_cartesian_topology");
}

void destroy_cartesian_topology(CartesianTopology&) {
    not_implemented("destroy_cartesian_topology");
}

HaloExchangePlan make_halo_exchange_plan(const CartesianTopology&,
                                          std::int64_t) {
    not_implemented("make_halo_exchange_plan");
}

void destroy_halo_exchange_plan(HaloExchangePlan&) {
    not_implemented("destroy_halo_exchange_plan");
}

void start_exchange(const HaloExchangePlan&, double*) {
    not_implemented("start_exchange");
}

void wait_exchange(const HaloExchangePlan&) {
    not_implemented("wait_exchange");
}

double allreduce_scalar(const CartesianTopology&, double, ReductionOp) {
    not_implemented("allreduce_scalar");
}

void allreduce_vector(const CartesianTopology&,
                       const double*,
                       double*,
                       std::int64_t,
                       ReductionOp) {
    not_implemented("allreduce_vector");
}

} // namespace nextdftb::hpc
