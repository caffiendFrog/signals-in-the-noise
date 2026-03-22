"""Tests for signals_in_the_noise.utils.logging_config."""

import logging

import pytest

from signals_in_the_noise.utils.logging_config import setup_logging


@pytest.fixture(autouse=True)
def reset_root_logger():
    """Remove any handlers added during a test so tests don't bleed into each other."""
    root = logging.getLogger()
    original_handlers = list(root.handlers)
    original_level = root.level
    yield
    root.handlers = original_handlers
    root.setLevel(original_level)


def test_setup_logging_attaches_handler():
    root = logging.getLogger()
    root.handlers.clear()
    setup_logging()
    assert len(root.handlers) == 1


def test_setup_logging_sets_info_level_by_default():
    root = logging.getLogger()
    root.handlers.clear()
    setup_logging()
    assert root.level == logging.INFO


def test_setup_logging_respects_custom_level():
    root = logging.getLogger()
    root.handlers.clear()
    setup_logging(level=logging.DEBUG)
    assert root.level == logging.DEBUG


def test_setup_logging_is_idempotent():
    root = logging.getLogger()
    root.handlers.clear()
    setup_logging()
    setup_logging()
    assert len(root.handlers) == 1


def test_setup_logging_handler_is_stream_handler():
    root = logging.getLogger()
    root.handlers.clear()
    setup_logging()
    assert isinstance(root.handlers[0], logging.StreamHandler)
