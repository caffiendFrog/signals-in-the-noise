import logging

logger = logging.getLogger(__name__)


def fdr_to_stars(fdr: float) -> str:
    """Convert an FDR q-value to a significance-star annotation string.

    Thresholds follow the conventional three-tier scheme:

    - ``< 0.01`` → ``' ***'``
    - ``< 0.05`` → ``' **'``
    - ``< 0.1``  → ``' *'``
    - ``>= 0.1`` → ``''`` (empty string, no annotation)

    Args:
        fdr: FDR-corrected q-value in the range ``[0, 1]``.

    Returns:
        A star-annotation string (with a leading space when non-empty) suitable
        for appending to axis tick labels or table cells.
    """
    if fdr < 0.01:
        return " ***"
    if fdr < 0.05:
        return " **"
    if fdr < 0.1:
        return " *"
    return ""
