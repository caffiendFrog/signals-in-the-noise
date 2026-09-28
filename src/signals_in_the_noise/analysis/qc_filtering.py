"""Cohort tables and frozen parameters for the BRCA1 QC-filtering analysis.

Notebook 00 selects the cohort and writes these files. Later notebooks read
them. Nothing here tests a genotype difference.
"""

import json
import re
from pathlib import Path

import pandas as pd
from anndata import AnnData

from signals_in_the_noise.analysis.confounders import specimen_covariates
from signals_in_the_noise.config import get_data_path
from signals_in_the_noise.preprocessing.qc import (
    EXCLUSIVE_REASON_PRIORITY,
    NOISE_FLAG_COLUMNS,
    QcThresholds,
    exclusive_reason,
    modal_thresholds,
)
from signals_in_the_noise.preprocessing.specimens import Cohort, patient_id

BRCA1_CONDITION = "BRCA1 pre-neoplastic"
WT_CONDITION = "Normal"
GENOTYPE_OF_CONDITION = {BRCA1_CONDITION: "BRCA1", WT_CONDITION: "WT"}
COHORT_CONDITIONS = (BRCA1_CONDITION, WT_CONDITION)
PREMENOPAUSAL = "Pre"
CELL_POPULATION = "Total"

PRIMARY_ALPHA = 0.01
SECONDARY_ALPHA = 0.02
PRIMARY_CONFIDENCE = 0.99
SECONDARY_CONFIDENCE = 0.98
LOO_CONFIDENCE = 0.98
MDE_POWER = 0.8
MAD_N_MADS = 3.0
NEAR_THRESHOLD_FRACTION = 0.05
THINNING_SEEDS = (0, 1, 2, 3, 4)

OUTPUT_DIRNAME = "brca1_qc_filtering"
CELLS_FILENAME = "cells.csv"
DONORS_FILENAME = "donors.csv"
PARAMS_FILENAME = "analysis_params.json"
HANDOFF_FILENAME = "00_handoff.md"

EXPECTED_DONORS = frozenset(
    {
        "B1-0894",
        "B1-0033",
        "B1-0023",
        "B1-0090",
        "N-0092",
        "N-0019",
        "N-0093",
        "N-0230.17",
        "N-0064",
        "N-0233",
        "N-0169",
        "N-0123",
    }
)
"""Donor ids for the 4 BRCA1 and 8 premenopausal WT Total samples."""

_DONOR_PREFIX = {BRCA1_CONDITION: "B1", WT_CONDITION: "N"}
_SOURCE_PREFIX = re.compile(r"_(?:N|B1)-([A-Z]+)")

CELL_COLUMNS = (
    "barcode",
    "donor",
    "specimen_id",
    "genotype",
    "qc_status",
    "is_noise",
    *NOISE_FLAG_COLUMNS,
    "exclusive_reason",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
)

PAL_MODEL_VARIANTS = (
    {"dispersion": "pearson", "test": "f"},
    {"dispersion": "deviance", "test": "f"},
    {"dispersion": "pearson", "test": "chi_square"},
    {"dispersion": "deviance", "test": "chi_square"},
)


def output_directory() -> Path:
    """Directory for this analysis, created if needed."""
    path = get_data_path("GSE161529_analysis_cache") / OUTPUT_DIRNAME
    path.mkdir(parents=True, exist_ok=True)
    return path


def brca1_wt_cohort(adatas) -> Cohort:
    """Total-population cohort: every BRCA1 donor, and premenopausal WT donors only.

    Post-oophorectomy BRCA1 donors stay in. The menopause filter applies to
    Normal specimens only.
    """
    return Cohort.from_objects(
        adatas,
        COHORT_CONDITIONS,
        cell_population=CELL_POPULATION,
        normal_menopause=PREMENOPAUSAL,
    )


def donor_label(adata: AnnData) -> str:
    """Short donor id, for example ``B1-0894`` or ``N-0230.17``."""
    condition = str(adata.uns["cancer_type"])
    try:
        prefix = _DONOR_PREFIX[condition]
    except KeyError as error:
        raise ValueError(f"No donor prefix for cancer type {condition!r}.") from error
    return f"{prefix}-{patient_id(adata)}"


def genotype_of(adata: AnnData) -> str:
    """``BRCA1`` or ``WT`` for a cohort specimen."""
    condition = str(adata.uns["cancer_type"])
    try:
        return GENOTYPE_OF_CONDITION[condition]
    except KeyError as error:
        raise ValueError(f"No genotype for cancer type {condition!r}.") from error


def inferred_source_prefix(specimen_id: str) -> str:
    """Source token parsed from a file stem, for example ``PM``, ``MH``, or ``KCF``.

    The token is inferred from the sample name. It is not a recorded batch.
    """
    match = _SOURCE_PREFIX.search(specimen_id)
    if match is None:
        raise ValueError(f"No source prefix in specimen id {specimen_id!r}.")
    return match.group(1)


def assert_published_cell_counts(adata: AnnData) -> None:
    """Require the observed cell counts to match the published GEO counts.

    Raises:
        ValueError: If either the total or the retained count disagrees.
    """
    observed_total = int(adata.n_obs)
    observed_retained = int((adata.obs["is_noise"] == 0).sum())
    published_before = int(adata.uns["num_cells_before"])
    published_after = int(adata.uns["num_cells_after"])
    if observed_total == published_before and observed_retained == published_after:
        return
    raise ValueError(
        f"{donor_label(adata)} cell counts do not match the published GEO counts: "
        f"observed {observed_total} cells and {observed_retained} retained, "
        f"published {published_before} before filtering and {published_after} after."
    )


def cell_table(cohort: Cohort) -> pd.DataFrame:
    """One row per cell in the cohort.

    ``specimen_id`` is included because 10x barcodes repeat across donors.
    ``exclusive_reason`` is missing for retained cells. ``is_noise`` is kept
    so existing QC helpers can subset the table.
    """
    frames = []
    for specimen_id, adata in cohort.specimens.items():
        assert_published_cell_counts(adata)
        obs = adata.obs
        flags = obs.loc[:, list(NOISE_FLAG_COLUMNS)].astype(int)
        recorded_noise = obs["is_noise"].astype(int)
        if not recorded_noise.equals(flags.any(axis=1).astype(int)):
            raise ValueError(
                f"{donor_label(adata)} has is_noise values that disagree with the four QC flags."
            )
        frame = pd.DataFrame(
            {
                "barcode": obs.index.astype(str),
                "donor": donor_label(adata),
                "specimen_id": specimen_id,
                "genotype": genotype_of(adata),
                "qc_status": obs["is_noise"]
                .astype(int)
                .map({0: "retained", 1: "removed"})
                .to_numpy(),
                "is_noise": obs["is_noise"].astype(int).to_numpy(),
                "exclusive_reason": exclusive_reason(flags).to_numpy(),
                "total_counts": obs["total_counts"].to_numpy(),
                "n_genes_by_counts": obs["n_genes_by_counts"].to_numpy(),
                "pct_counts_mt": obs["pct_counts_mt"].to_numpy(),
            },
            index=obs.index,
        )
        for column in NOISE_FLAG_COLUMNS:
            frame[column] = flags[column].to_numpy()
        frames.append(frame.loc[:, CELL_COLUMNS])
    table = pd.concat(frames, ignore_index=True)
    table["exclusive_reason"] = table["exclusive_reason"].astype("string")
    return table


def donor_table(cohort: Cohort) -> pd.DataFrame:
    """One row per donor, indexed by specimen id.

    Starts from :func:`specimen_covariates` and adds the genotype, parity,
    inferred source, published GEO counts, and raw median depth and gene count
    over all barcodes.
    """
    table = specimen_covariates(cohort.specimens)
    extra = {}
    for specimen_id, adata in cohort.specimens.items():
        retained = int((adata.obs["is_noise"] == 0).sum())
        n_cells = int(adata.n_obs)
        extra[specimen_id] = {
            "donor": donor_label(adata),
            "genotype": genotype_of(adata),
            "parity": str(adata.uns["parity"]),
            "source_prefix": inferred_source_prefix(specimen_id),
            "num_cells_before": int(adata.uns["num_cells_before"]),
            "num_cells_after": int(adata.uns["num_cells_after"]),
            "n_retained": retained,
            "n_removed": n_cells - retained,
            "median_total_counts": float(adata.obs["total_counts"].median()),
            "median_n_genes_by_counts": float(adata.obs["n_genes_by_counts"].median()),
        }
    table = table.join(pd.DataFrame.from_dict(extra, orient="index"))
    table["removed_fraction"] = table["n_removed"] / table["n_cells"]
    return table.sort_values(["genotype", "donor"])


def validate_cohort(donors: pd.DataFrame) -> None:
    """Require the donor set, and that every WT donor is premenopausal.

    Raises:
        ValueError: If the donor ids differ from :data:`EXPECTED_DONORS`,
            a donor id is repeated, or a WT donor is not premenopausal.
    """
    found = set(donors["donor"])
    if found != set(EXPECTED_DONORS):
        missing = sorted(set(EXPECTED_DONORS) - found)
        extra = sorted(found - set(EXPECTED_DONORS))
        raise ValueError(f"Cohort donors do not match the plan. Missing {missing}, extra {extra}.")
    if donors["donor"].duplicated().any():
        raise ValueError("A donor id is repeated.")
    wt = donors.loc[donors["genotype"] == "WT", "menopause_status"]
    if not (wt == PREMENOPAUSAL).all():
        raise ValueError("A WT donor is not premenopausal.")


def _uniform_thresholds(cohort: Cohort) -> QcThresholds:
    return modal_thresholds(QcThresholds.from_uns(adata.uns) for adata in cohort.specimens.values())


def _depth_target(donors: pd.DataFrame) -> float:
    """Lowest per-donor median UMI count. Thinning can only reduce depth."""
    return float(donors["median_total_counts"].min())


def phase1_analysis_params(cohort: Cohort, donors: pd.DataFrame) -> dict:
    """Phase 1 parameters, including the values computed from this cohort.

    The uniform thresholds and the depth target are properties of the cohort,
    not test results. ``c`` is left null for notebook 02.
    """
    uniform = _uniform_thresholds(cohort)
    return {
        "alpha": {"primary": PRIMARY_ALPHA, "secondary": SECONDARY_ALPHA},
        "confidence": {
            "primary": PRIMARY_CONFIDENCE,
            "secondary": SECONDARY_CONFIDENCE,
            "loo": LOO_CONFIDENCE,
        },
        "tests": {
            "overall_removed_fraction": {
                "role": "primary",
                "alternative": "two-sided",
                "family": None,
                "alpha": PRIMARY_ALPHA,
            },
            "per_reason_nonexclusive": {
                "role": "secondary",
                "alternative": "two-sided",
                "family": "reasons",
                "reasons": [label for _, label in EXCLUSIVE_REASON_PRIORITY],
                "adjustment": "bh",
                "alpha": SECONDARY_ALPHA,
            },
            "exclusive_reasons": {"role": "descriptive"},
            "lp_baseline_and_tipping_point": {
                "role": "primary",
                "alternative": "greater",
                "contrast": "BRCA1 > WT",
                "alpha": PRIMARY_ALPHA,
            },
        },
        "loo": {
            "estimate": True,
            "interval": "two-sided",
            "level": LOO_CONFIDENCE,
            "applies_to_every_run": True,
        },
        "mde": {
            "method": "exact_enumeration",
            "n_assignments": 495,
            "power": MDE_POWER,
            "scale": "log_odds",
            "continuity": "(x + 0.5) / (n + 1)",
            "centering": "center each genotype on its own mean, then add the pooled mean",
            "report": "percentage_points",
            "cross_check": "normal_approximation",
            "optional_sensitivity": "resample donors with replacement within genotype",
        },
        "robustness": {
            "applies_to": ["uniform", "depth_equalized"],
            "excluded_from_verdict": ["mad"],
            "conclusion_unchanged_at_alpha": PRIMARY_ALPHA,
            "estimate_within": "half the primary MDE, on the log-odds scale",
            "depth_seed_spread": "if not small relative to MDE/2, report cannot_be_judged",
            "diagnostic": "spearman rank correlation of per-donor removed fractions",
        },
        "pal_proportion_check": {
            "published_p": 0.14,
            "pass": "exact donor-by-cluster counts, and 0.14 inside the variant p-value range",
            "model": "counts ~ donor + cluster + cluster:genotype",
            "family": "quasi-Poisson",
            "variants": list(PAL_MODEL_VARIANTS),
        },
        "sensitivity": {
            "uniform_thresholds": {
                "rule": "modal",
                "tie_break": "smallest modal value",
                **{name: float(getattr(uniform, name)) for name in QcThresholds.UNS_KEYS},
            },
            "mad": {
                "n_mads": MAD_N_MADS,
                "genes_and_counts_scale": "log1p",
                "mito_scale": "fraction",
                "genes_and_counts_sides": "both",
                "center": "median",
                "scope": "all barcodes of the sample",
                "in_robustness_verdict": False,
            },
            "depth_equalized": {
                "target_rule": "lowest per-donor median total_counts over all barcodes",
                "target_median_umi": _depth_target(donors),
                "qc_rule_after_thinning": "uniform",
                "seeds": list(THINNING_SEEDS),
                "generator": "numpy.random.default_rng",
            },
        },
        "near_threshold_band": {
            "fraction": NEAR_THRESHOLD_FRACTION,
            "population": "retained",
            "per": "sample and threshold",
            "thresholds": [label for _, label in EXCLUSIVE_REASON_PRIORITY],
            "pool": "all donors",
            "genotype_blind": True,
            "c": None,
            "c_definition": "largest enrichment of any mapped cell type",
        },
        "tipping_point": {
            "wt_removed_lp_share": "that donor's retained LP share",
            "brca1_removed_lp_share": "min(1, k times that donor's retained LP share)",
            "k_max": "1 / smallest BRCA1 retained LP share",
            "zones": [
                {"name": "plausible", "when": "k_star <= c", "phase_2": "go"},
                {
                    "name": "possible_but_implausible",
                    "when": "c < k_star <= k_max",
                    "phase_2": "confirmation",
                },
                {"name": "impossible", "when": "k_star > k_max", "phase_2": "optional"},
            ],
        },
        "donors": donors["donor"].tolist(),
        "batch": None,
    }


def cohort_handoff(donors: pd.DataFrame, params: dict, output_dir: Path) -> str:
    """Markdown summary to paste back after the notebook runs on AWS."""
    sensitivity = params["sensitivity"]
    uniform = sensitivity["uniform_thresholds"]
    depth = sensitivity["depth_equalized"]["target_median_umi"]
    shown = donors.reset_index().loc[
        :,
        [
            "donor",
            "genotype",
            "menopause_status",
            "parity",
            "source_prefix",
            "n_cells",
            "n_retained",
            "n_removed",
            "removed_fraction",
            "median_total_counts",
        ],
    ]
    shown = shown.copy()
    shown["removed_fraction"] = shown["removed_fraction"].map(lambda value: f"{value:.4f}")
    shown["median_total_counts"] = shown["median_total_counts"].map(lambda value: f"{value:.1f}")
    post = donors.loc[donors["menopause_status"] != PREMENOPAUSAL, ["donor", "menopause_status"]]
    post_lines = (
        "\n".join(f"- {row.donor}: {row.menopause_status}" for row in post.itertuples()) or "- none"
    )
    counts = donors["genotype"].value_counts()
    return "\n".join(
        [
            "# Notebook 00 handoff",
            "",
            "Paste this file back into the coding session. It records the cohort and the",
            "parameters frozen before any genotype test. It does not contain a test result.",
            "",
            f"Output directory: `{output_dir}`",
            "",
            "Files:",
            "",
            f"- `{CELLS_FILENAME}`: one row per cell",
            f"- `{DONORS_FILENAME}`: one row per donor",
            f"- `{PARAMS_FILENAME}`: Phase 1 section",
            "",
            f"Donors: {int(counts.get('BRCA1', 0))} BRCA1 and {int(counts.get('WT', 0))} WT.",
            "",
            _markdown_table(shown),
            "",
            "Menopause differs from premenopausal:",
            "",
            post_lines,
            "",
            "Source prefix is inferred from the sample name. No batch or processing date exists.",
            "",
            "Frozen uniform (modal) thresholds:",
            "",
            f"- mito_upper: {uniform['mito_upper']}",
            f"- genes_lower: {uniform['genes_lower']}",
            f"- genes_upper: {uniform['genes_upper']}",
            f"- total_upper: {uniform['total_upper']}",
            "",
            f"Frozen depth target (lowest per-donor median UMI count): {depth:.1f}",
            f"Thinning seeds: {list(THINNING_SEEDS)}",
            "",
            "Not computed in this notebook: the near-threshold constant c (notebook 02),",
            "and every genotype test (notebook 01).",
            "",
        ]
    )


def write_phase1_outputs(cohort: Cohort, output_dir: Path | None = None) -> str:
    """Write the cell table, donor table, Phase 1 parameters, and the handoff.

    An existing ``phase_2`` section in the parameters file is left in place.

    Returns:
        The handoff text.
    """
    output_dir = output_directory() if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cells = cell_table(cohort)
    donors = donor_table(cohort)
    validate_cohort(donors)
    params = phase1_analysis_params(cohort, donors)
    cells.to_csv(output_dir / CELLS_FILENAME, index=False)
    donors.to_csv(output_dir / DONORS_FILENAME)
    _write_phase1_params(output_dir / PARAMS_FILENAME, params)
    handoff = cohort_handoff(donors, params, output_dir)
    (output_dir / HANDOFF_FILENAME).write_text(handoff, encoding="utf-8")
    return handoff


def _write_phase1_params(path: Path, phase1: dict) -> None:
    existing: dict = {}
    if path.exists():
        existing = json.loads(path.read_text(encoding="utf-8"))
    existing["phase_1"] = phase1
    path.write_text(json.dumps(existing, indent=2) + "\n", encoding="utf-8")


def _markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    body = [
        "| " + " | ".join(str(row[column]) for column in frame.columns) + " |"
        for _, row in frame.iterrows()
    ]
    return "\n".join([header, separator, *body])
