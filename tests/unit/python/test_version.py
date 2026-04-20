"""Tests for the package version surface."""

from __future__ import annotations

import re

import nextdftb


_SEMVER_LIKE = re.compile(r"^\d+\.\d+(\.\d+)?([+-].+)?$")


def test_version_returns_non_empty_string():
    v = nextdftb.version()
    assert isinstance(v, str)
    assert v.strip() != ""


def test_version_matches_semver_like_pattern():
    assert _SEMVER_LIKE.match(nextdftb.version()) is not None


def test_dunder_version_matches_version_function():
    assert nextdftb.__version__ == nextdftb.version()


def test_version_is_stable_across_calls():
    assert nextdftb.version() == nextdftb.version()
