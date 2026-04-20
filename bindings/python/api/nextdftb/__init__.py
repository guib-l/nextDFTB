"""nextDFTB — Python API.

Thin ergonomic wrapper over the compiled C++20 / Fortran 2018 engine.
All numerical work happens in native code; this layer only handles:

- Python exceptions
- Integration with the Python :mod:`logging` module via the unified sink
"""

from __future__ import annotations

from . import _core
from ._core import (
    init,
    finalize,
    is_initialized,
    test,
    version,
)
from ._log_handler import install_python_logging_bridge

# The `log` submodule of `_core` is re-exported so users can do:
#     import nextdftb
#     nextdftb.log.open("run.log", also_stdout=True)
log = _core.log

__all__ = [
    "init",
    "finalize",
    "is_initialized",
    "test",
    "version",
    "log",
    "install_python_logging_bridge",
]

__version__ = version()
