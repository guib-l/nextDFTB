#include "nextdftb/buffer.hpp"

#include <cstdlib>
#include <cstring>
#include <new>

#if !defined(_WIN32)
#  include <cstdlib>   // posix_memalign on POSIX
#endif

namespace nextdftb {

namespace {

void* aligned_alloc_bytes(std::size_t n_bytes, std::size_t alignment) {
    std::size_t padded = ((n_bytes + alignment - 1) / alignment) * alignment;
    if (padded == 0) padded = alignment;
    void* p = nullptr;
#if defined(_WIN32)
    p = _aligned_malloc(padded, alignment);
#else
    if (posix_memalign(&p, alignment, padded) != 0) {
        p = nullptr;
    }
#endif
    if (!p) throw std::bad_alloc{};
    return p;
}

void aligned_free(void* p) noexcept {
    if (!p) return;
#if defined(_WIN32)
    _aligned_free(p);
#else
    std::free(p);
#endif
}

} // namespace

AlignedBuffer::AlignedBuffer(std::size_t n_doubles) : n_(n_doubles) {
    if (n_ == 0) return;
    data_ = static_cast<double*>(
        aligned_alloc_bytes(n_ * sizeof(double), ALIGNMENT));
    std::memset(data_, 0, n_ * sizeof(double));
}

AlignedBuffer::~AlignedBuffer() { aligned_free(data_); }

AlignedBuffer::AlignedBuffer(AlignedBuffer&& other) noexcept
    : data_(other.data_), n_(other.n_) {
    other.data_ = nullptr;
    other.n_    = 0;
}

AlignedBuffer& AlignedBuffer::operator=(AlignedBuffer&& other) noexcept {
    if (this != &other) {
        aligned_free(data_);
        data_       = other.data_;
        n_          = other.n_;
        other.data_ = nullptr;
        other.n_    = 0;
    }
    return *this;
}

void AlignedBuffer::fill(double v) noexcept {
    for (std::size_t i = 0; i < n_; ++i) data_[i] = v;
}

void AlignedBuffer::resize(std::size_t n_doubles) {
    if (n_doubles == n_) return;
    aligned_free(data_);
    data_ = nullptr;
    n_    = 0;
    if (n_doubles == 0) return;
    data_ = static_cast<double*>(
        aligned_alloc_bytes(n_doubles * sizeof(double), ALIGNMENT));
    n_ = n_doubles;
    std::memset(data_, 0, n_ * sizeof(double));
}

} // namespace nextdftb
