"""Tests for signals_in_the_noise.utils.log."""

import logging

from signals_in_the_noise.utils.log import get_logger


def test_get_logger_returns_logger_instance():
    result = get_logger("test.module")
    assert isinstance(result, logging.Logger)


def test_get_logger_uses_provided_name():
    result = get_logger("my.custom.name")
    assert result.name == "my.custom.name"


def test_get_logger_dunder_name_convention():
    result = get_logger(__name__)
    assert result.name == __name__


def test_get_logger_does_not_add_handlers():
    """Library loggers must not attach handlers — configuration belongs at the entry point."""
    result = get_logger("handler.check.module")
    assert len(result.handlers) == 0


def test_get_logger_same_name_returns_same_instance():
    first = get_logger("shared.logger")
    second = get_logger("shared.logger")
    assert first is second
