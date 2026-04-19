// nextDFTB — I/O abstraction (placeholder).
//
// The framework will eventually ship an HDF5 backend; this interface exists
// now so the rest of the code never takes a direct dependency on the
// output format. The only concrete backend delivered in this phase is the
// text writer (::nextdftb::io::TextWriter).
//
// Future backends (HDF5, NetCDF, …) implement the same Writer interface
// without any change to call sites.

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>

namespace nextdftb::io {

enum class Backend { Text /*, HDF5 (deferred) */ };

class Writer {
public:
    virtual ~Writer() = default;

    virtual void open(const std::string& path) = 0;
    virtual void close()                       = 0;

    // Scalar records.
    virtual void write_scalar(std::string_view name, double      value) = 0;
    virtual void write_scalar(std::string_view name, std::int64_t value) = 0;

    // 1-D float64 dataset.
    virtual void write_vector(std::string_view name,
                              const double*    data,
                              std::int64_t     n) = 0;

    virtual void flush() = 0;
};

// Factory — returns a concrete backend. Currently only Backend::Text
// is implemented; anything else throws std::logic_error.
std::unique_ptr<Writer> make_writer(Backend b);

} // namespace nextdftb::io
