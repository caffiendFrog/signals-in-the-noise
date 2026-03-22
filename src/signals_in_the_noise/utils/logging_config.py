"""Logging configuration for the signals-in-the-noise package.

Call ``setup_logging()`` once at the start of any entry point (script,
notebook, or interactive session) to route package log records to the
console.  Library modules must never import or call this module — they
only call ``get_logger(__name__)`` and let the entry point decide how
records are handled.
"""

import logging


def setup_logging(level: int = logging.INFO) -> None:
    """Configure a console handler on the root logger.

    Safe to call multiple times; subsequent calls are no-ops once a
    handler is already attached to the root logger.

    Args:
        level: Minimum log level to emit.  Defaults to ``logging.INFO``.
            Pass ``logging.DEBUG`` for verbose output.
    """
    root = logging.getLogger()
    if root.handlers:
        return

    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
    )
    root.addHandler(handler)
    root.setLevel(level)
