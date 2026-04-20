"""Tests for the runtime lifecycle: init / finalize / is_initialized."""

from __future__ import annotations

import pytest

import nextdftb


@pytest.fixture(autouse=True)
def _ensure_finalized():
    """Make sure the runtime is finalized between every test in this file."""
    if nextdftb.is_initialized():
        nextdftb.finalize()
    yield
    if nextdftb.is_initialized():
        nextdftb.finalize()


def test_initial_state_is_not_initialized():
    assert nextdftb.is_initialized() is False


def test_init_then_is_initialized_true():
    nextdftb.init()
    assert nextdftb.is_initialized() is True


def test_finalize_resets_initialized_state():
    nextdftb.init()
    nextdftb.finalize()
    assert nextdftb.is_initialized() is False


def test_init_with_explicit_thread_count():
    nextdftb.init(num_threads=1)
    assert nextdftb.is_initialized() is True


def test_init_with_zero_uses_default():
    nextdftb.init(num_threads=0)
    assert nextdftb.is_initialized() is True


def test_init_returns_none():
    assert nextdftb.init() is None


def test_finalize_returns_none():
    nextdftb.init()
    assert nextdftb.finalize() is None


def test_is_initialized_returns_bool():
    assert isinstance(nextdftb.is_initialized(), bool)
