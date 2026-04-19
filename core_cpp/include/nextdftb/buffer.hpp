// nextDFTB — aligned RAII buffer for shared simulation data.
//
// C++ owns allocations for buffers that cross the ABI to Fortran. Every
// shared buffer is aligned to 64 bytes (cache line / AVX-512). This type
// is move-only; copy is forbidden so ownership semantics stay explicit.

#pragma once

#include <cstddef>
#include <cstdint>

namespace nextdftb {

class AlignedBuffer {
public:
    static constexpr std::size_t ALIGNMENT = 64;

    AlignedBuffer() = default;
    explicit AlignedBuffer(std::size_t n_doubles);
    ~AlignedBuffer();

    AlignedBuffer(const AlignedBuffer&)            = delete;
    AlignedBuffer& operator=(const AlignedBuffer&) = delete;

    AlignedBuffer(AlignedBuffer&& other) noexcept;
    AlignedBuffer& operator=(AlignedBuffer&& other) noexcept;

    double*       data()       noexcept { return data_; }
    const double* data() const noexcept { return data_; }
    std::size_t   size() const noexcept { return n_; }

    void fill(double v) noexcept;
    void resize(std::size_t n_doubles);

private:
    double*     data_ = nullptr;
    std::size_t n_    = 0;
};

} // namespace nextdftb
