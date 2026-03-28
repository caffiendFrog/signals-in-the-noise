"""Tests for signals_in_the_noise.io.gmt."""

import pytest

from signals_in_the_noise.io.gmt import load_gmt


def _write_gmt(tmp_path, name: str, description: str, genes: list[str]) -> object:
    """Write a minimal GMT file and return its Path."""
    content = "\t".join([name, description] + genes)
    gmt_file = tmp_path / f"{name}.gmt"
    gmt_file.write_text(content, encoding="utf-8")
    return gmt_file


# ---------------------------------------------------------------------------
# Contract tests
# ---------------------------------------------------------------------------


def test_load_gmt_returns_list(tmp_path):
    path = _write_gmt(tmp_path, "HALLMARK_G2M_CHECKPOINT", "http://example.com", ["BRCA1", "TP53"])
    result = load_gmt(path)
    assert isinstance(result, list)


def test_load_gmt_returns_gene_strings(tmp_path):
    genes = ["BRCA1", "TP53", "MYC"]
    path = _write_gmt(tmp_path, "TEST_PATHWAY", "http://example.com", genes)
    result = load_gmt(path)
    assert all(isinstance(g, str) for g in result)


def test_load_gmt_skips_pathway_name_and_description(tmp_path):
    """GMT descriptions are URLs (no spaces); split() handles them as a single token."""
    genes = ["GENE_A", "GENE_B", "GENE_C"]
    path = _write_gmt(tmp_path, "PATHWAY_NAME", "http://example.com/pathway", genes)
    result = load_gmt(path)
    assert "PATHWAY_NAME" not in result
    assert "http://example.com/pathway" not in result


def test_load_gmt_returns_correct_genes(tmp_path):
    genes = ["BRCA1", "TP53", "MYC"]
    path = _write_gmt(tmp_path, "TEST_PATHWAY", "http://example.com", genes)
    result = load_gmt(path)
    assert result == genes


def test_load_gmt_handles_tab_separated_file(tmp_path):
    """Standard GMT files use tab separators; split() handles these identically."""
    content = "HALLMARK_G2M_CHECKPOINT\thttp://example.com\tBRCA1\tTP53\n"
    gmt_file = tmp_path / "test.gmt"
    gmt_file.write_text(content, encoding="utf-8")
    result = load_gmt(gmt_file)
    assert result == ["BRCA1", "TP53"]


def test_load_gmt_accepts_path_like_string(tmp_path):
    path = _write_gmt(tmp_path, "PATHWAY", "desc", ["GENE1"])
    result = load_gmt(str(path))
    assert result == ["GENE1"]


# ---------------------------------------------------------------------------
# Error tests
# ---------------------------------------------------------------------------


def test_load_gmt_raises_on_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_gmt(tmp_path / "nonexistent.gmt")


def test_load_gmt_raises_on_file_with_fewer_than_three_tokens(tmp_path):
    gmt_file = tmp_path / "bad.gmt"
    gmt_file.write_text("PATHWAY\tdescription", encoding="utf-8")
    with pytest.raises(ValueError, match="fewer than 3 tokens"):
        load_gmt(gmt_file)
