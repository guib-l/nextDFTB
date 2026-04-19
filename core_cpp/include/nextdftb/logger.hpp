// nextDFTB — unified logger (C++ wrapper around the ABI sink).
//
// The sink is a thread-safe singleton: one file handle, one mutex. All
// four layers (Fortran, C++, Python, CLI) converge through
// nextdftb_log(level, layer, function, message).
//
// Log format emitted by the sink:
//   [YYYY-MM-DDThh:mm:ss.mmm] [LEVEL] [layer:function] message

#pragma once

#include <string>

#include "nextdftb/logging.h"

namespace nextdftb {

class Logger {
public:
    static Logger& instance();

    void open(const std::string& path, bool also_stdout);
    void close();
    void set_level(int level);
    int  get_level() const;

    void log(int level,
             const char* layer,
             const char* function,
             const char* message);

    Logger(const Logger&)            = delete;
    Logger& operator=(const Logger&) = delete;

private:
    Logger();
    ~Logger();

    struct Impl;
    Impl* p_ = nullptr;
};

// Convenience wrappers for C++ call sites.
inline void log_debug(const char* function, const std::string& message) {
    nextdftb_log(NEXTDFTB_LOG_DEBUG, "cpp", function, message.c_str());
}
inline void log_info(const char* function, const std::string& message) {
    nextdftb_log(NEXTDFTB_LOG_INFO, "cpp", function, message.c_str());
}
inline void log_warn(const char* function, const std::string& message) {
    nextdftb_log(NEXTDFTB_LOG_WARN, "cpp", function, message.c_str());
}
inline void log_error(const char* function, const std::string& message) {
    nextdftb_log(NEXTDFTB_LOG_ERROR, "cpp", function, message.c_str());
}

} // namespace nextdftb
