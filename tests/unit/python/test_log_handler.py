"""Tests for the Python logging → C ABI bridge."""

from __future__ import annotations

import logging

import nextdftb
from nextdftb._log_handler import NextdftbHandler, _LEVEL_MAP


def test_install_returns_handler_instance(isolated_logger):
    handler = nextdftb.install_python_logging_bridge(isolated_logger.name)
    try:
        assert isinstance(handler, NextdftbHandler)
        assert handler in isolated_logger.handlers
    finally:
        isolated_logger.removeHandler(handler)


def test_install_uses_default_logger_name():
    handler = nextdftb.install_python_logging_bridge()
    try:
        default_logger = logging.getLogger("nextdftb")
        assert handler in default_logger.handlers
    finally:
        logging.getLogger("nextdftb").removeHandler(handler)


def test_handler_level_map_covers_standard_levels():
    assert _LEVEL_MAP[logging.DEBUG] == nextdftb.log.DEBUG
    assert _LEVEL_MAP[logging.INFO] == nextdftb.log.INFO
    assert _LEVEL_MAP[logging.WARNING] == nextdftb.log.WARN
    assert _LEVEL_MAP[logging.ERROR] == nextdftb.log.ERROR
    assert _LEVEL_MAP[logging.CRITICAL] == nextdftb.log.ERROR


def test_handler_emit_routes_message_to_sink(
    tmp_path, saved_log_level, isolated_logger
):
    log_path = tmp_path / "bridge.log"
    nextdftb.log.open(str(log_path), also_stdout=False)
    nextdftb.log.set_level(nextdftb.log.DEBUG)
    handler = nextdftb.install_python_logging_bridge(isolated_logger.name)
    try:
        isolated_logger.setLevel(logging.DEBUG)
        isolated_logger.info("bridge-marker-info")
        isolated_logger.error("bridge-marker-error")
    finally:
        isolated_logger.removeHandler(handler)
        nextdftb.log.close()

    content = log_path.read_text()
    assert "bridge-marker-info" in content
    assert "bridge-marker-error" in content


def test_handler_emit_uses_only_message_format(
    tmp_path, saved_log_level, isolated_logger
):
    """install_python_logging_bridge sets a formatter that emits only %(message)s."""
    log_path = tmp_path / "bridge.log"
    nextdftb.log.open(str(log_path), also_stdout=False)
    nextdftb.log.set_level(nextdftb.log.DEBUG)
    handler = nextdftb.install_python_logging_bridge(isolated_logger.name)
    try:
        isolated_logger.setLevel(logging.INFO)
        isolated_logger.info("clean-msg-%s", "interpolated")
    finally:
        isolated_logger.removeHandler(handler)
        nextdftb.log.close()

    content = log_path.read_text()
    assert "clean-msg-interpolated" in content


def test_handler_emit_swallows_unknown_level(
    tmp_path, saved_log_level, isolated_logger
):
    """Records with non-standard levelno must default to INFO without raising."""
    log_path = tmp_path / "bridge.log"
    nextdftb.log.open(str(log_path), also_stdout=False)
    nextdftb.log.set_level(nextdftb.log.DEBUG)
    handler = nextdftb.install_python_logging_bridge(isolated_logger.name)
    try:
        record = logging.LogRecord(
            name=isolated_logger.name,
            level=42,
            pathname=__file__,
            lineno=1,
            msg="custom-level-marker",
            args=(),
            exc_info=None,
            func="custom_fn",
        )
        handler.handle(record)
    finally:
        isolated_logger.removeHandler(handler)
        nextdftb.log.close()

    content = log_path.read_text()
    assert "custom-level-marker" in content


def test_log_submodule_reexport_is_same_object():
    assert nextdftb.log is nextdftb._core.log
