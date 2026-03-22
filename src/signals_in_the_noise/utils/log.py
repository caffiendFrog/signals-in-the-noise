import logging


def get_logger(name: str) -> logging.Logger:
    """Return a module-level logger.

    Args:
        name: The logger name, typically ``__name__``.

    Returns:
        A Logger instance. Handler configuration is left to the application
        entry point; this function never calls ``basicConfig`` or attaches
        handlers.
    """
    return logging.getLogger(name)
