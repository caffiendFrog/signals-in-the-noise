"""Tests for signals_in_the_noise.analysis.statistics."""

import pytest

from signals_in_the_noise.analysis.statistics import fdr_to_stars


# ---------------------------------------------------------------------------
# Return-type contract
# ---------------------------------------------------------------------------


def test_fdr_to_stars_returns_string():
    assert isinstance(fdr_to_stars(0.05), str)


# ---------------------------------------------------------------------------
# Threshold boundary tests
# ---------------------------------------------------------------------------


def test_fdr_to_stars_below_0_01_returns_three_stars():
    assert fdr_to_stars(0.001) == " ***"


def test_fdr_to_stars_exactly_0_01_returns_two_stars():
    assert fdr_to_stars(0.01) == " **"


def test_fdr_to_stars_between_0_01_and_0_05_returns_two_stars():
    assert fdr_to_stars(0.03) == " **"


def test_fdr_to_stars_exactly_0_05_returns_one_star():
    assert fdr_to_stars(0.05) == " *"


def test_fdr_to_stars_between_0_05_and_0_1_returns_one_star():
    assert fdr_to_stars(0.07) == " *"


def test_fdr_to_stars_exactly_0_1_returns_empty_string():
    assert fdr_to_stars(0.1) == ""


def test_fdr_to_stars_above_0_1_returns_empty_string():
    assert fdr_to_stars(0.5) == ""


def test_fdr_to_stars_at_zero_returns_three_stars():
    assert fdr_to_stars(0.0) == " ***"


def test_fdr_to_stars_at_one_returns_empty_string():
    assert fdr_to_stars(1.0) == ""
