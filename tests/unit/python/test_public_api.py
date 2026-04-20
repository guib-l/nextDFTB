"""Surface-level checks on the nextdftb public API."""

from __future__ import annotations

import nextdftb


_EXPECTED_EXPORTS = {
    "init",
    "finalize",
    "is_initialized",
    "test",
    "version",
    "log",
    "install_python_logging_bridge",
}


def test_all_advertises_expected_names():
    assert set(nextdftb.__all__) == _EXPECTED_EXPORTS


def test_all_names_are_attributes():
    for name in nextdftb.__all__:
        assert hasattr(nextdftb, name), f"missing public attribute: {name}"


def test_callables_are_callable():
    for name in (
        "init",
        "finalize",
        "is_initialized",
        "test",
        "version",
        "install_python_logging_bridge",
    ):
        assert callable(getattr(nextdftb, name))


def test_log_submodule_exposes_required_symbols():
    for name in ("DEBUG", "INFO", "WARN", "ERROR",
                 "open", "close", "set_level", "get_level", "log"):
        assert hasattr(nextdftb.log, name), f"log.{name} missing"
