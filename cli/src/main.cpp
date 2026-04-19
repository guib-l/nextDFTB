// nextDFTB — standalone CLI.
//
// Single global executable that links every layer (core_cpp -> abi ->
// kernels_fortran) and can run the infrastructure smoke test without
// Python.
//
// Invocation:
//     nextdftb [OPTIONS]
//
// Options:
//     --threads N         OpenMP thread count (0 = leave default)
//     --log PATH          Log file path                        (default: stdout only)
//     --log-level LVL     DEBUG|INFO|WARN|ERROR                (default INFO)
//     --log-stdout        Also echo to stdout when --log is set
//     --output PATH       Write the final value to this file
//     --version           Print version and exit
//     --help              Print this usage and exit

#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

#include "nextdftb/abi.h"
#include "nextdftb/logger.hpp"

namespace {

struct CliOptions {
    int          threads      = 0;
    std::string  log_path;
    std::string  log_level    = "INFO";
    bool         log_stdout   = false;
    std::string  output_path;
    bool         show_help    = false;
    bool         show_version = false;
};

void print_usage(std::ostream& os) {
    os <<
        "Usage: nextdftb [OPTIONS]\n"
        "\n"
        "  --threads N         OpenMP thread count (0 = leave default)\n"
        "  --log PATH          Log file path (default: stdout only)\n"
        "  --log-level LVL     DEBUG|INFO|WARN|ERROR (default INFO)\n"
        "  --log-stdout        Also echo to stdout when --log is set\n"
        "  --output PATH       Write the final value to this file\n"
        "  --version           Print version and exit\n"
        "  --help              Print this usage and exit\n";
}

bool parse_int(const char* s, int& out) {
    if (!s || !*s) return false;
    char* end = nullptr;
    errno = 0;
    long long v = std::strtoll(s, &end, 10);
    if (errno || !end || *end != '\0') return false;
    out = static_cast<int>(v);
    return true;
}

int parse_log_level(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (s == "DEBUG") return NEXTDFTB_LOG_DEBUG;
    if (s == "INFO")  return NEXTDFTB_LOG_INFO;
    if (s == "WARN" || s == "WARNING") return NEXTDFTB_LOG_WARN;
    if (s == "ERROR") return NEXTDFTB_LOG_ERROR;
    return -1;
}

bool parse_args(int argc, char** argv, CliOptions& o) {
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        auto need_value = [&](const char* name, const char*& val) -> bool {
            if (i + 1 >= argc) {
                std::cerr << "nextdftb: missing value for " << name << "\n";
                return false;
            }
            val = argv[++i];
            return true;
        };

        if (a == "--help" || a == "-h") {
            o.show_help = true;
        } else if (a == "--version") {
            o.show_version = true;
        } else if (a == "--log-stdout") {
            o.log_stdout = true;
        } else if (a == "--threads") {
            const char* v = nullptr;
            if (!need_value("--threads", v)) return false;
            if (!parse_int(v, o.threads)) {
                std::cerr << "nextdftb: invalid --threads value\n"; return false;
            }
        } else if (a == "--log") {
            const char* v = nullptr;
            if (!need_value("--log", v)) return false;
            o.log_path = v;
        } else if (a == "--log-level") {
            const char* v = nullptr;
            if (!need_value("--log-level", v)) return false;
            o.log_level = v;
        } else if (a == "--output") {
            const char* v = nullptr;
            if (!need_value("--output", v)) return false;
            o.output_path = v;
        } else {
            std::cerr << "nextdftb: unknown option '" << a << "'\n";
            return false;
        }
    }
    return true;
}

void dump_last_abi_error(std::ostream& os) {
    char loc[256] = {};
    char msg[512] = {};
    nextdftb_get_last_error_location(loc, sizeof(loc));
    nextdftb_get_last_error_message (msg, sizeof(msg));
    os << "[nextdftb] error: code=" << nextdftb_get_last_error_code()
       << " severity=" << nextdftb_get_last_error_severity()
       << " at [" << loc << "]: " << msg << "\n";
}

} // namespace

int main(int argc, char** argv) {
    CliOptions opt;
    if (!parse_args(argc, argv, opt)) {
        print_usage(std::cerr);
        return 2;
    }
    if (opt.show_help) { print_usage(std::cout); return 0; }
    if (opt.show_version) {
        std::cout << "nextdftb " << nextdftb_version_string() << "\n";
        return 0;
    }

    // --- Logging -----------------------------------------------------
    const int lvl = parse_log_level(opt.log_level);
    if (lvl < 0) {
        std::cerr << "nextdftb: invalid --log-level '" << opt.log_level << "'\n";
        return 2;
    }
    nextdftb_log_set_level(lvl);

    if (!opt.log_path.empty()) {
        const int s = nextdftb_log_open(opt.log_path.c_str(),
                                          opt.log_stdout ? 1 : 0);
        if (s != NEXTDFTB_OK) {
            dump_last_abi_error(std::cerr);
            return 3;
        }
    }

    nextdftb::log_info("main",
                       "nextdftb CLI starting, version "
                       + std::string(nextdftb_version_string()));

    // --- Runtime init -----------------------------------------------
    int s = nextdftb_init(opt.threads);
    if (s != NEXTDFTB_OK) {
        dump_last_abi_error(std::cerr);
        return 4;
    }

    // --- Run smoke test ---------------------------------------------
    double out = 0.0;
    s = nextdftb_test(&out);
    if (s != NEXTDFTB_OK) {
        dump_last_abi_error(std::cerr);
        (void)nextdftb_finalize();
        return 5;
    }

    // --- Output -----------------------------------------------------
    std::cout << "test = " << out << "\n";
    if (!opt.output_path.empty()) {
        std::ofstream f(opt.output_path);
        if (!f) {
            std::cerr << "nextdftb: cannot open output file '"
                      << opt.output_path << "'\n";
            (void)nextdftb_finalize();
            return 6;
        }
        f << "test " << out << "\n";
    }

    nextdftb::log_info("main", "nextdftb CLI done");
    (void)nextdftb_finalize();
    nextdftb_log_close();
    return 0;
}
