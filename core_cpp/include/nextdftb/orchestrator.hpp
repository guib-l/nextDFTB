// nextDFTB — C++ orchestrator.
//
// Thin RAII wrapper over the Fortran runtime lifecycle. Simulation-level
// pipelines live in their own drivers (to be added per physics module).

#pragma once

namespace nextdftb {

class Orchestrator {
public:
    Orchestrator();
    ~Orchestrator();

    void init(int num_threads);
    void finalize();
    bool initialized() const noexcept;

    Orchestrator(const Orchestrator&)            = delete;
    Orchestrator& operator=(const Orchestrator&) = delete;
};

} // namespace nextdftb
