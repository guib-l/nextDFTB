# `bindings/python/` — Python bindings

Two sub-directories, split by responsibility:

```
python/
├── abi/   pybind11 C++ glue. Turns the C++ orchestrator and the C ABI
│         into the compiled extension `nextdftb._core`. Handles
│         exception translation, NumPy layout/dtype normalization,
│         submodule registration.
└── api/   Pure-Python package `nextdftb/`. Ergonomic wrappers,
          logging bridge, docstrings, type hints. No native code,
          no build step — just what end users import.
```

The dependency direction is one-way: `api/` imports from `abi/`
(via the compiled `_core` module). `abi/` never imports from `api/`.

End-user surface area:

```python
import nextdftb
nextdftb.init()
nextdftb.test()
nextdftb.finalize()
```

is defined by `api/nextdftb/__init__.py`. Everything else is an
implementation detail.

---

## Adding a new feature

Assume the Fortran side and the C++ wrapper are already in place
(see `engine/README.md` and `core_cpp/`). To expose the function to
Python:

### 1. Bind the symbol in `abi/module.cpp`

Every `m.def(...)` lives here. Responsibilities of the binding:

- Normalize NumPy arrays to F-contiguous `float64` with
  `py::array::f_style | py::array::forcecast`.
- Validate dimensions that the C++ layer will not re-check.
- Call the C++ wrapper from `nextdftb::...` (not the raw C ABI) so
  errors land as `nextdftb::Error` and are turned into
  `ValueError` / `MemoryError` / `RuntimeError` by the existing
  exception translator.
- Release the GIL (`py::call_guard<py::gil_scoped_release>()`) for
  any function that performs non-trivial work in native code.

```cpp
m.def("saxpy",
      [](double a,
         py::array_t<double, py::array::f_style | py::array::forcecast> x,
         py::array_t<double, py::array::f_style | py::array::forcecast> y) {
          if (x.size() != y.size()) {
              throw std::invalid_argument("saxpy: x and y must have the same length");
          }
          nextdftb::saxpy(a, x.mutable_data(), y.mutable_data(),
                          static_cast<std::int64_t>(x.size()));
      },
      py::arg("a"), py::arg("x"), py::arg("y"),
      py::call_guard<py::gil_scoped_release>(),
      "In-place y := a*x + y.");
```

Keep bindings thin. No control flow, no numerical work, no feature
flags. If it would live in `if`/`for`, it belongs in C++.

### 2. Re-export from `api/nextdftb/__init__.py`

Bring the new function into the public namespace:

```python
from ._core import saxpy
__all__ = [..., "saxpy"]
```

If the function needs ergonomic smoothing (default arguments,
NumPy reshapes, type coercion, docstring enrichment), write a small
wrapper in the package instead of re-exporting the raw binding:

```python
def saxpy(a: float, x: np.ndarray, y: np.ndarray) -> None:
    """In-place y := a*x + y."""
    _core.saxpy(float(a), np.ascontiguousarray(x, dtype=np.float64),
                          np.ascontiguousarray(y, dtype=np.float64))
```

Rule of thumb: one-liner re-exports live in `__init__.py`; anything
longer gets its own module under `api/nextdftb/` (e.g. `_saxpy.py`).

### 3. Staging and install

`configure_file` stages each `*.py` into the build tree. When you
add a new Python file under `api/nextdftb/`, append a matching
`configure_file(... COPYONLY)` entry to `bindings/CMakeLists.txt`
(look for the existing `__init__.py` / `_log_handler.py` lines).
The `install(DIRECTORY ...)` rule already globs `*.py` so install
does not need to be touched.

### 4. Tests

Python tests run under pytest from `tests/runners/`. Cover:

- the happy path,
- dtype/layout coercion (pass a C-order `float32` array and make
  sure it still works),
- error translation (pass a bad shape and check that a `ValueError`
  is raised, not a bare `RuntimeError`).

---

## Conventions

- **NumPy memory layout**: the engine is column-major. `abi/` is the
  only place allowed to force-copy — never do it in `api/`.
- **No Python-side math**. Loops, reductions, and array fills happen
  in Fortran. `api/` may reshape, broadcast, or validate, but must
  not iterate over large arrays.
- **Logging**: use the `nextdftb` logger and keep
  `install_python_logging_bridge()` as the single entry point so
  records flow into the unified C ABI sink.
- **Type hints**: add them to any function exposed in `__init__.py`.
- **Backwards compatibility**: everything in `__all__` is a public
  API. Renaming a binding is a breaking change — deprecate the old
  name by keeping it as an alias until a major release.
