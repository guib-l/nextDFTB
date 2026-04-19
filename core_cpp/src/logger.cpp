#include "nextdftb/logger.hpp"
#include "nextdftb/errors.h"
#include "nextdftb/logging.h"

#include <chrono>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>

namespace nextdftb {

struct Logger::Impl {
    std::mutex    mu;
    std::ofstream file;
    bool          also_stdout = true;   // default: stdout-only sink until opened
    int           level       = NEXTDFTB_LOG_INFO;
};

Logger::Logger()  : p_(new Impl()) {}
Logger::~Logger() { delete p_; }

Logger& Logger::instance() {
    static Logger g;
    return g;
}

void Logger::open(const std::string& path, bool also_stdout) {
    std::lock_guard<std::mutex> lock(p_->mu);
    if (p_->file.is_open()) p_->file.close();
    p_->also_stdout = also_stdout;
    if (!path.empty()) {
        p_->file.open(path, std::ios::out | std::ios::app);
    }
}

void Logger::close() {
    std::lock_guard<std::mutex> lock(p_->mu);
    if (p_->file.is_open()) p_->file.close();
}

void Logger::set_level(int level) {
    std::lock_guard<std::mutex> lock(p_->mu);
    p_->level = level;
}

int Logger::get_level() const { return p_->level; }

namespace {

const char* level_name(int lvl) {
    switch (lvl) {
        case NEXTDFTB_LOG_DEBUG: return "DEBUG";
        case NEXTDFTB_LOG_INFO:  return "INFO ";
        case NEXTDFTB_LOG_WARN:  return "WARN ";
        case NEXTDFTB_LOG_ERROR: return "ERROR";
        default:                 return "?????";
    }
}

std::string timestamp_iso8601() {
    using namespace std::chrono;
    const auto now = system_clock::now();
    const auto t   = system_clock::to_time_t(now);
    const auto ms  = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

    std::tm tm_buf{};
#if defined(_WIN32)
    localtime_s(&tm_buf, &t);
#else
    localtime_r(&t, &tm_buf);
#endif
    char buf[64];
    std::snprintf(buf, sizeof(buf),
                  "%04d-%02d-%02dT%02d:%02d:%02d.%03lld",
                  tm_buf.tm_year + 1900,
                  tm_buf.tm_mon  + 1,
                  tm_buf.tm_mday,
                  tm_buf.tm_hour, tm_buf.tm_min, tm_buf.tm_sec,
                  static_cast<long long>(ms.count()));
    return std::string(buf);
}

} // namespace

void Logger::log(int level,
                 const char* layer,
                 const char* function,
                 const char* message) {
    if (level < p_->level) return;

    std::string line;
    line.reserve(200);
    line += '[';
    line += timestamp_iso8601();
    line += "] [";
    line += level_name(level);
    line += "] [";
    line += layer    ? layer    : "?";
    line += ':';
    line += function ? function : "?";
    line += "] ";
    line += message  ? message  : "";
    line += '\n';

    std::lock_guard<std::mutex> lock(p_->mu);
    if (p_->file.is_open()) {
        p_->file << line;
        p_->file.flush();
    }
    if (p_->also_stdout || !p_->file.is_open()) {
        std::cout << line;
        std::cout.flush();
    }
}

} // namespace nextdftb

// ---------------------------------------------------------------------
// C ABI shims
// ---------------------------------------------------------------------
extern "C" {

int nextdftb_log_open(const char* path, int also_stdout) {
    try {
        nextdftb::Logger::instance().open(path ? path : "", also_stdout != 0);
        return NEXTDFTB_OK;
    } catch (...) {
        nextdftb_set_error(NEXTDFTB_ERR_IO, NEXTDFTB_SEV_RECOVERABLE,
                            "cpp", "nextdftb_log_open",
                            "failed to open log file");
        return NEXTDFTB_ERR_IO;
    }
}

void nextdftb_log_close(void) {
    nextdftb::Logger::instance().close();
}

void nextdftb_log_set_level(int level) {
    nextdftb::Logger::instance().set_level(level);
}

int nextdftb_log_get_level(void) {
    return nextdftb::Logger::instance().get_level();
}

void nextdftb_log(int level,
                   const char* layer,
                   const char* function,
                   const char* message) {
    nextdftb::Logger::instance().log(level, layer, function, message);
}

} // extern "C"
