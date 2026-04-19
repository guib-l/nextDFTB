// nextDFTB — C++ exception type + translation between ABI status codes
// and C++ exceptions.
//
// Rule of thumb:
//   * Fortran kernels publish errors via nextdftb_set_error and return a
//     negative status code.
//   * C++ catches the status at the ABI call site via throw_if_abi_error,
//     which pulls the last-error state and throws nextdftb::Error.
//   * pybind11 translates nextdftb::Error into RuntimeError / ValueError.

#pragma once

#include <stdexcept>
#include <string>

#include "nextdftb/errors.h"

namespace nextdftb {

class Error : public std::runtime_error {
public:
    Error(int code, std::string layer, std::string function, std::string msg);

    int                code()     const noexcept { return code_; }
    int                severity() const noexcept { return severity_; }
    const std::string& layer()    const noexcept { return layer_; }
    const std::string& function() const noexcept { return function_; }

    // Push this error onto the thread-local ABI error state so that C
    // callers and Python (via the getters) can read it back.
    void publish() const noexcept;

private:
    int         code_;
    int         severity_;
    std::string layer_;
    std::string function_;
};

// If status != NEXTDFTB_OK, read the last-error state from the ABI and
// throw nextdftb::Error. `fallback_function` is used if the ABI state is
// empty (e.g. early-return without set_error).
void throw_if_abi_error(int status, const char* fallback_function);

} // namespace nextdftb
