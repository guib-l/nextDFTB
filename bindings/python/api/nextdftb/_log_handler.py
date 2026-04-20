"""Bridge between :mod:`logging` and the unified nextDFTB C ABI sink.

Installing :class:`NextdftbHandler` on a Python logger causes every record
emitted on that logger to flow through the same file and formatter used by
the Fortran and C++ layers.
"""

from __future__ import annotations

import logging as _py_logging

from . import _core

_LEVEL_MAP = {
    _py_logging.DEBUG:    _core.log.DEBUG,
    _py_logging.INFO:     _core.log.INFO,
    _py_logging.WARNING:  _core.log.WARN,
    _py_logging.ERROR:    _core.log.ERROR,
    _py_logging.CRITICAL: _core.log.ERROR,
}


class NextdftbHandler(_py_logging.Handler):
    """Route Python log records into the unified C ABI log sink."""

    def emit(self, record: _py_logging.LogRecord) -> None:
        try:
            msg = self.format(record)
            level = _LEVEL_MAP.get(record.levelno, _core.log.INFO)
            func = record.funcName or "?"
            _core.log.log(level, func, msg)
        except Exception:  # noqa: BLE001
            self.handleError(record)


def install_python_logging_bridge(logger_name: str = "nextdftb") -> NextdftbHandler:
    """Attach a :class:`NextdftbHandler` to the named logger and return it."""
    handler = NextdftbHandler()
    handler.setFormatter(_py_logging.Formatter("%(message)s"))
    _py_logging.getLogger(logger_name).addHandler(handler)
    return handler
