#include "nextdftb/io_writer.hpp"

#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

namespace nextdftb::io {

namespace {

class TextWriter final : public Writer {
public:
    void open(const std::string& path) override {
        close();
        file_.open(path, std::ios::out | std::ios::trunc);
        if (!file_) {
            throw std::runtime_error("TextWriter: cannot open '" + path + "'");
        }
    }

    void close() override {
        if (file_.is_open()) file_.close();
    }

    void write_scalar(std::string_view name, double value) override {
        ensure_open();
        file_ << name << " " << value << "\n";
    }

    void write_scalar(std::string_view name, std::int64_t value) override {
        ensure_open();
        file_ << name << " " << value << "\n";
    }

    void write_vector(std::string_view name,
                       const double*    data,
                       std::int64_t     n) override {
        ensure_open();
        file_ << name << " " << n;
        for (std::int64_t i = 0; i < n; ++i) file_ << " " << data[i];
        file_ << "\n";
    }

    void flush() override {
        if (file_.is_open()) file_.flush();
    }

private:
    void ensure_open() {
        if (!file_.is_open()) {
            throw std::runtime_error("TextWriter: no file open");
        }
    }

    std::ofstream file_;
};

} // namespace

std::unique_ptr<Writer> make_writer(Backend b) {
    switch (b) {
        case Backend::Text: return std::make_unique<TextWriter>();
    }
    throw std::logic_error("nextdftb::io::make_writer: unknown backend");
}

} // namespace nextdftb::io
