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


def combine_gmt_files(
    input_paths: list[Path],
    output_path: Path,
    overwrite: bool = False,
) -> Path:
    """Concatenate multiple GMT files into a single combined GMT file.

    Each input file is written to the output followed by a newline, preserving
    the original GMT content verbatim.  This mirrors the format expected by
    tools such as :mod:`gseapy` that accept a pre-merged gene-set file.

    When ``overwrite`` is ``False`` (the default) and ``output_path`` already
    exists, the function returns immediately without touching the file.  Set
    ``overwrite=True`` to force regeneration.

    Args:
        input_paths: Ordered sequence of paths to individual ``.gmt`` files.
        output_path: Destination path for the combined output file.  Parent
            directories must already exist.
        overwrite: If ``False``, skip writing when ``output_path`` exists.
            Defaults to ``False``.

    Returns:
        The resolved ``output_path``.

    Raises:
        FileNotFoundError: If any path in ``input_paths`` does not exist.
    """
    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        logger.debug("combined GMT already exists at %s, skipping", output_path)
        return output_path

    with output_path.open("w", encoding="utf-8") as out_fh:
        for gmt_path in input_paths:
            gmt_path = Path(gmt_path)
            out_fh.write(gmt_path.read_text(encoding="utf-8"))
            out_fh.write("\n")

    logger.debug("wrote combined GMT with %d entries to %s", len(input_paths), output_path)
    return output_path
