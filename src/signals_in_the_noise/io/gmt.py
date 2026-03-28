import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def load_gmt(path: Path) -> list[str]:
    """Parse a single-entry GMT file and return the gene names.

    A GMT file has the format::

        PATHWAY_NAME<tab>description<tab>GENE1<tab>GENE2<tab>...

    The pathway name (token 0) and description (token 1) are skipped;
    only the gene identifiers are returned.

    Args:
        path: Path to the ``.gmt`` file.

    Returns:
        List of gene name strings extracted from the file.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file contains fewer than three whitespace-separated
            tokens (i.e. is missing gene entries).
    """
    path = Path(path)
    tokens = path.read_text(encoding="utf-8").split()
    if len(tokens) < 3:
        raise ValueError(
            f"GMT file {path.name!r} has fewer than 3 tokens; expected at least "
            "a pathway name, a description, and one gene."
        )
    genes = tokens[2:]
    logger.debug("loaded %d genes from %s", len(genes), path.name)
    return genes
