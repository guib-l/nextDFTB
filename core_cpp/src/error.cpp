#include "nextdftb/errors.h"
#include "nextdftb/error.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>

namespace {

struct ThreadErrorState {
    int         code      = NEXTDFTB_OK;
    int         severity  = NEXTDFTB_SEV_RECOVERABLE;
    std::string message;
    std::string location;   // "layer:function"
};

thread_local ThreadErrorState g_err;

void copy_to_c(const std::string& src, char* dst, std::int32_t dst_len) {
    if (!dst || dst_len <= 0) return;
    const std::size_t cap = static_cast<std::size_t>(dst_len - 1);
    const std::size_t n   = std::min(src.size(), cap);
    std::memcpy(dst, src.data(), n);
    dst[n] = '\0';
}

} // namespace

extern "C" {

int nextdftb_get_last_error_code(void) { return g_err.code; }

int nextdftb_get_last_error_severity(void) { return g_err.severity; }

void nextdftb_get_last_error_message(char* buffer, std::int32_t buffer_len) {
    copy_to_c(g_err.message, buffer, buffer_len);
}

void nextdftb_get_last_error_location(char* buffer, std::int32_t buffer_len) {
    copy_to_c(g_err.location, buffer, buffer_len);
}

void nextdftb_clear_error(void) {
    g_err.code     = NEXTDFTB_OK;
    g_err.severity = NEXTDFTB_SEV_RECOVERABLE;
    g_err.message.clear();
    g_err.location.clear();
}

void nextdftb_set_error(int         code,
                         int         severity,
                         const char* layer,
                         const char* function,
                         const char* message) {
    g_err.code     = code;
    g_err.severity = severity;
    g_err.message.assign(message ? message : "");
    g_err.location.clear();
    g_err.location.append(layer    ? layer    : "?");
    g_err.location.push_back(':');
    g_err.location.append(function ? function : "?");
}

} // extern "C"

namespace nextdftb {

Error::Error(int code, std::string layer, std::string function, std::string msg)
    : std::runtime_error(msg),
      code_(code),
      severity_(code == NEXTDFTB_ERR_FATAL ? NEXTDFTB_SEV_FATAL
                                            : NEXTDFTB_SEV_RECOVERABLE),
      layer_(std::move(layer)),
      function_(std::move(function)) {}

void Error::publish() const noexcept {
    nextdftb_set_error(code_,
                        severity_,
                        layer_.c_str(),
                        function_.c_str(),
                        what());
}

void throw_if_abi_error(int status, const char* fallback_function) {
    if (status == NEXTDFTB_OK) return;

    char loc[256] = {};
    char msg[512] = {};
    nextdftb_get_last_error_location(loc, sizeof(loc));
    nextdftb_get_last_error_message (msg, sizeof(msg));

    std::string location = loc[0] ? loc : std::string();
    std::string layer, function;
    auto colon = location.find(':');
    if (colon != std::string::npos) {
        layer    = location.substr(0, colon);
        function = location.substr(colon + 1);
    } else {
        layer    = "cpp";
        function = fallback_function ? fallback_function : "?";
    }

    throw Error(status,
                std::move(layer),
                std::move(function),
                msg[0] ? msg : "unknown ABI error");
}

} // namespace nextdftb
