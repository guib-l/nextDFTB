// nextDFTB — pybind11 bindings.
//
// Exposes the compiled engine to Python as `nextdftb._core`. The public
// Python surface is defined by `nextdftb/__init__.py` (pure-Python wrapper).
//
// Memory-layout strategy: the engine is Fortran-centric (column-major).
// Every array parameter is normalized with `py::array::f_style | forcecast`
// before crossing the ABI. If the incoming array is already F-contiguous
// and dtype=float64, no copy is made; otherwise pybind11 produces a
// transposed copy transparently so the Python user never has to specify
// `order='F'`.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>

#include "nextdftb/abi.h"
#include "nextdftb/error.hpp"
#include "nextdftb/orchestrator.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() =
        "nextDFTB low-level bindings (pybind11). "
        "The public API is exposed via the `nextdftb` Python package.";

    // ------------------------------------------------------------
    // Exception translation: nextdftb::Error → Python exceptions.
    // ------------------------------------------------------------
    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        }
        catch (const nextdftb::Error& e) {
            switch (e.code()) {
                case NEXTDFTB_ERR_INVALID_ARG:
                case NEXTDFTB_ERR_DIMENSION_MISMATCH:
                    PyErr_SetString(PyExc_ValueError, e.what());
                    return;
                case NEXTDFTB_ERR_ALLOCATION:
                    PyErr_SetString(PyExc_MemoryError, e.what());
                    return;
                case NEXTDFTB_ERR_NOT_IMPLEMENTED:
                    PyErr_SetString(PyExc_NotImplementedError, e.what());
                    return;
                case NEXTDFTB_ERR_IO:
                    PyErr_SetString(PyExc_IOError, e.what());
                    return;
                default:
                    PyErr_SetString(PyExc_RuntimeError, e.what());
                    return;
            }
        }
    });

    // ------------------------------------------------------------
    // Runtime lifecycle
    // ------------------------------------------------------------
    m.def("init", [](int num_threads) {
        const int s = nextdftb_init(num_threads);
        nextdftb::throw_if_abi_error(s, "init");
    }, py::arg("num_threads") = 0,
       "Initialize the Fortran runtime. num_threads<=0 leaves OpenMP default.");

    m.def("finalize", []() {
        const int s = nextdftb_finalize();
        nextdftb::throw_if_abi_error(s, "finalize");
    });

    m.def("is_initialized", []() { return nextdftb_is_initialized() != 0; });

    // ------------------------------------------------------------
    // Infrastructure smoke test
    // ------------------------------------------------------------
    m.def("test",
          []() -> double {
              double out = 0.0;
              const int s = nextdftb_test(&out);
              nextdftb::throw_if_abi_error(s, "test");
              return out;
          },
          "Run the infrastructure smoke test; returns a deterministic value.");

    // ------------------------------------------------------------
    // Version
    // ------------------------------------------------------------
    m.def("version", []() { return std::string(nextdftb_version_string()); });

    // ------------------------------------------------------------
    // Unified logging submodule (Python entry into the C ABI sink)
    // ------------------------------------------------------------
    py::module_ logmod = m.def_submodule("log", "Unified logging bridge.");
    logmod.attr("DEBUG") = static_cast<int>(NEXTDFTB_LOG_DEBUG);
    logmod.attr("INFO")  = static_cast<int>(NEXTDFTB_LOG_INFO);
    logmod.attr("WARN")  = static_cast<int>(NEXTDFTB_LOG_WARN);
    logmod.attr("ERROR") = static_cast<int>(NEXTDFTB_LOG_ERROR);

    logmod.def("open",
               [](const std::string& path, bool also_stdout) {
                   const int s = nextdftb_log_open(path.c_str(),
                                                     also_stdout ? 1 : 0);
                   nextdftb::throw_if_abi_error(s, "log.open");
               },
               py::arg("path"), py::arg("also_stdout") = false);

    logmod.def("close",     []() { nextdftb_log_close(); });
    logmod.def("set_level", [](int lvl) { nextdftb_log_set_level(lvl); });
    logmod.def("get_level", []() { return nextdftb_log_get_level(); });

    logmod.def("log",
               [](int level, const std::string& function, const std::string& message) {
                   nextdftb_log(level, "python", function.c_str(), message.c_str());
               },
               py::arg("level"), py::arg("function"), py::arg("message"));
}
