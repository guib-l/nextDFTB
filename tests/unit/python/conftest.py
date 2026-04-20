"""Shared fixtures for the nextdftb Python API unit tests."""

from __future__ import annotations

import logging

import pytest

import nextdftb


@pytest.fixture
def runtime():
    """Initialize the native runtime for the duration of a test."""
    nextdftb.init()
    try:
        yield nextdftb
    finally:
        nextdftb.finalize()


@pytest.fixture
def saved_log_level():
    """Snapshot and restore the global log level around a test."""
    original = nextdftb.log.get_level()
    try:
        yield original
    finally:
        nextdftb.log.set_level(original)


@pytest.fixture
def isolated_logger():
    """Provide a logger with no inherited handlers; clean up afterwards."""
    name = "nextdftb.tests.isolated"
    logger = logging.getLogger(name)
    logger.handlers.clear()
    logger.propagate = False
    previous_level = logger.level
    try:
        yield logger
    finally:
        logger.handlers.clear()
        logger.setLevel(previous_level)
