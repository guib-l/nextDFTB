"""Tests for nextdftb.test() — the infrastructure smoke test."""

from __future__ import annotations

import math

import pytest

import nextdftb


# Closed form: sum_{i=1..1000} i = 1000 * 1001 / 2.
_EXPECTED_SMOKE = 500500.0


def test_smoke_returns_expected_value(runtime):
    assert math.isclose(runtime.test(), _EXPECTED_SMOKE, rel_tol=1e-12)


def test_smoke_returns_float(runtime):
    assert isinstance(runtime.test(), float)


def test_smoke_is_deterministic(runtime):
    assert runtime.test() == runtime.test()


def test_smoke_requires_initialization():
    """Calling test() before init() must not silently succeed."""
    if nextdftb.is_initialized():
        nextdftb.finalize()
    with pytest.raises(Exception):
        nextdftb.test()
