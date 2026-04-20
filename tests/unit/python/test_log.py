"""Tests for the unified logging submodule (nextdftb.log)."""

from __future__ import annotations

import nextdftb


def test_log_level_constants_distinct_and_ordered():
    levels = [
        nextdftb.log.DEBUG,
        nextdftb.log.INFO,
        nextdftb.log.WARN,
        nextdftb.log.ERROR,
    ]
    assert all(isinstance(level, int) for level in levels)
    assert levels == sorted(levels)
    assert len(set(levels)) == 4


def test_log_level_constants_have_expected_values():
    assert nextdftb.log.DEBUG == 0
    assert nextdftb.log.INFO == 1
    assert nextdftb.log.WARN == 2
    assert nextdftb.log.ERROR == 3


def test_set_and_get_level_roundtrip(saved_log_level):
    for lvl in (
        nextdftb.log.DEBUG,
        nextdftb.log.INFO,
        nextdftb.log.WARN,
        nextdftb.log.ERROR,
    ):
        nextdftb.log.set_level(lvl)
        assert nextdftb.log.get_level() == lvl


def test_log_emits_without_open_sink(saved_log_level):
    """Logging without an open file sink must be a no-op, not crash."""
    nextdftb.log.set_level(nextdftb.log.DEBUG)
    nextdftb.log.log(nextdftb.log.INFO, "test_func", "hello world")


def test_log_open_writes_to_file(tmp_path, saved_log_level):
    log_path = tmp_path / "engine.log"
    nextdftb.log.open(str(log_path), also_stdout=False)
    try:
        nextdftb.log.set_level(nextdftb.log.DEBUG)
        nextdftb.log.log(nextdftb.log.INFO, "test_func", "marker-payload")
    finally:
        nextdftb.log.close()

    assert log_path.exists()
    content = log_path.read_text()
    assert "marker-payload" in content


def test_log_close_is_idempotent():
    nextdftb.log.close()
    nextdftb.log.close()


def test_log_filters_below_threshold(tmp_path, saved_log_level):
    log_path = tmp_path / "engine.log"
    nextdftb.log.open(str(log_path), also_stdout=False)
    try:
        nextdftb.log.set_level(nextdftb.log.WARN)
        nextdftb.log.log(nextdftb.log.DEBUG, "fn", "should-be-filtered-out")
        nextdftb.log.log(nextdftb.log.ERROR, "fn", "should-survive")
    finally:
        nextdftb.log.close()

    content = log_path.read_text()
    assert "should-be-filtered-out" not in content
    assert "should-survive" in content
