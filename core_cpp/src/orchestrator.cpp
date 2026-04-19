#include "nextdftb/orchestrator.hpp"

#include "nextdftb/abi.h"
#include "nextdftb/error.hpp"

namespace nextdftb {

Orchestrator::Orchestrator()  = default;
Orchestrator::~Orchestrator() = default;

void Orchestrator::init(int num_threads) {
    const int s = ::nextdftb_init(num_threads);
    throw_if_abi_error(s, "Orchestrator::init");
}

void Orchestrator::finalize() {
    const int s = ::nextdftb_finalize();
    throw_if_abi_error(s, "Orchestrator::finalize");
}

bool Orchestrator::initialized() const noexcept {
    return ::nextdftb_is_initialized() != 0;
}

} // namespace nextdftb
