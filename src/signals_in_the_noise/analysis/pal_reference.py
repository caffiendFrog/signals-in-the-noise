"""Pal's Fig 4B labels: reproduce P = 0.14, and freeze the near-threshold enrichment.

The proportion check uses the raw cluster ids. It passes only when the exported
cells are exactly the retained cells and 0.14 lies inside the range of four
quasi-Poisson variants. The enrichment ``c`` is a fold, on the same scale as
the later tipping-point multiplier, and it is computed without looking at genotype.
"""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.genmod.families import Poisson

from signals_in_the_noise.analysis.qc_filtering import (
    PAL_MODEL_VARIANTS,
    PARAMS_FILENAME,
    output_directory,
)

PAL_RETAINED_CELLS = 59766
PUBLISHED_P = 0.14
LABELS_FILENAME = "pal_norm_b1_labels.csv"
COUNTS_FILENAME = "02_cluster_counts.csv"
VARIANTS_FILENAME = "02_variants.csv"
ENRICHMENT_FILENAME = "02_enrichment.csv"
HANDOFF_FILENAME = "02_handoff.md"

_BARCODE = re.compile(r"([ACGT]{10,})-(\d+)$", re.IGNORECASE)
_LIBRARY = re.compile(r"(?:^|_)((?:N|B1)-[A-Z]+[0-9.]+)")
_THRESHOLD_COLUMNS = {
    "genes_lower": ("n_genes_by_counts", "qc_genes_lower"),
    "genes_upper": ("n_genes_by_counts", "qc_genes_upper"),
    "mito_upper": ("pct_counts_mt", "qc_mito_upper"),
    "total_upper": ("total_counts", "qc_total_upper"),
}


def normalize_barcode(value: str) -> str:
    """10x barcode ``ACGT...-1``, dropping any sample prefix Pal stored in front."""
    match = _BARCODE.search(str(value).strip())
    if match is None:
        raise ValueError(f"No 10x barcode in {value!r}.")
    return f"{match.group(1).upper()}-{match.group(2)}"


def cell_sample_name(cell_name: str, sample: str) -> str:
    """Sample id Pal stored in the cell name.

    ``orig.ident`` is only ``N`` or ``B1``, because Seurat keeps the first
    underscore-separated token. The cell name is ``N_0019_total_<barcode>``
    or ``B1_0023_<barcode>``. A bare 10x barcode keeps the sample column.
    """
    text = str(cell_name).strip()
    match = _BARCODE.search(text)
    if match is None or match.start() == 0:
        return str(sample)
    prefix = text[: match.start()].rstrip("_-")
    return prefix or str(sample)


def map_pal_samples(samples, donors: pd.DataFrame) -> pd.Series:
    """Map each Pal sample name to one donor id.

    A name matches a donor id (``N-0093`` or ``N_0019_total``) or the library
    token stored in the specimen id (``N-PM0095``). ``B1-MH0023`` and
    ``N-MH0023`` stay apart.
    """
    frame = _donor_frame(donors)
    library = {
        _library_token(specimen): donor
        for specimen, donor in zip(frame["specimen_id"], frame["donor"], strict=True)
    }
    donor_ids = {str(donor).upper(): donor for donor in frame["donor"]}
    mapped = []
    for sample in samples:
        key = _normalize_name(sample)
        found = []
        if key in donor_ids:
            found.append(donor_ids[key])
        for token, donor in library.items():
            if token is not None and token in key and donor not in found:
                found.append(donor)
        if len(found) != 1:
            raise ValueError(f"Pal sample {sample!r} matches donors {found}.")
        mapped.append(found[0])
    index = samples.index if isinstance(samples, pd.Series) else None
    return pd.Series(mapped, index=index, dtype="object")


def label_counts_match(
    labels: pd.DataFrame,
    cells: pd.DataFrame,
    donors: pd.DataFrame,
    *,
    expected_total: int = PAL_RETAINED_CELLS,
) -> dict:
    """Whether Pal's cells are exactly this cohort's retained cells.

    The comparison is the set of ``(donor, barcode)`` pairs, the per-donor
    retained counts, and the published total of 59,766.
    """
    exported = _label_keys(labels)
    retained = _retained_keys(cells)
    donor_frame = _donor_frame(donors).set_index("donor")
    per_donor = []
    for donor, exported_n in exported.groupby("donor").size().items():
        retained_n = int((retained["donor"] == donor).sum())
        published_n = (
            int(donor_frame.loc[donor, "n_retained"]) if donor in donor_frame.index else None
        )
        per_donor.append(
            {
                "donor": donor,
                "exported": int(exported_n),
                "retained": retained_n,
                "published_retained": published_n,
            }
        )
    exported_pairs = set(zip(exported["donor"], exported["barcode"], strict=True))
    retained_pairs = set(zip(retained["donor"], retained["barcode"], strict=True))
    only_exported = len(exported_pairs - retained_pairs)
    only_retained = len(retained_pairs - exported_pairs)
    totals_match = (
        len(exported) == expected_total
        and len(retained) == expected_total
        and only_exported == 0
        and only_retained == 0
        and all(
            row["exported"] == row["retained"] == row["published_retained"] for row in per_donor
        )
    )
    return {
        "match": totals_match,
        "n_exported": len(exported),
        "n_retained": len(retained),
        "expected_total": expected_total,
        "only_exported": only_exported,
        "only_retained": only_retained,
        "per_donor": pd.DataFrame(per_donor),
    }


def donor_cluster_counts(labels: pd.DataFrame) -> pd.DataFrame:
    """One row per donor and cluster, with zeros where a donor has none of that cluster."""
    required = {"donor", "genotype", "cluster"}
    missing = required - set(labels.columns)
    if missing:
        raise ValueError(f"Labels are missing columns {sorted(missing)}.")
    genotype = labels.groupby("donor")["genotype"].nunique()
    if genotype.max() != 1:
        raise ValueError("A donor has more than one genotype.")
    donor_genotype = labels.groupby("donor")["genotype"].first()
    clusters = sorted(labels["cluster"].astype(str).unique())
    donors = list(donor_genotype.index)
    full = pd.MultiIndex.from_product([donors, clusters], names=["donor", "cluster"])
    counts = (
        labels.assign(cluster=labels["cluster"].astype(str))
        .groupby(["donor", "cluster"], observed=True)
        .size()
        .reindex(full, fill_value=0)
        .rename("count")
        .reset_index()
    )
    counts["genotype"] = counts["donor"].map(donor_genotype)
    return counts


def quasi_poisson_interaction_p(counts: pd.DataFrame, *, dispersion: str, test: str) -> float:
    """P-value for ``cluster:genotype`` in a quasi-Poisson model.

    The mean model is ``count ~ donor + cluster + cluster:genotype``. Donor
    already contains genotype, so genotype is not a separate main effect.
    Dispersion is Pearson or deviance, taken from the full model. The test is
    an F test or a chi-square test on the deviance difference from the model
    without the interaction.
    """
    if dispersion not in {"pearson", "deviance"}:
        raise ValueError(f"Unknown dispersion {dispersion!r}.")
    if test not in {"f", "chi_square"}:
        raise ValueError(f"Unknown test {test!r}.")
    full = _fit_poisson(counts, interaction=True)
    reduced = _fit_poisson(counts, interaction=False)
    if full.df_resid < 1:
        raise ValueError("The full model has no residual degrees of freedom.")
    scale = (full.pearson_chi2 if dispersion == "pearson" else full.deviance) / full.df_resid
    df_diff = reduced.df_resid - full.df_resid
    dev_diff = max(0.0, float(reduced.deviance - full.deviance))
    if scale <= 0:
        return 1.0 if dev_diff == 0 else 0.0
    statistic = dev_diff / scale
    if test == "chi_square":
        return float(stats.chi2.sf(statistic, df_diff))
    f_statistic = statistic / df_diff
    return float(stats.f.sf(f_statistic, df_diff, full.df_resid))


def variant_p_values(counts: pd.DataFrame, variants=PAL_MODEL_VARIANTS) -> pd.DataFrame:
    """The four pre-specified quasi-Poisson variants, in the stored order."""
    rows = []
    for variant in variants:
        p_value = quasi_poisson_interaction_p(
            counts, dispersion=variant["dispersion"], test=variant["test"]
        )
        rows.append({**variant, "p_value": p_value})
    return pd.DataFrame(rows)


def published_p_in_range(p_values, published: float = PUBLISHED_P) -> bool:
    """True when ``published`` lies between the smallest and largest variant p-value."""
    values = [float(value) for value in p_values]
    if not values or not all(np.isfinite(values)):
        return False
    return min(values) <= published <= max(values)


def reference_passes(*, counts_match: bool, p_values, published: float = PUBLISHED_P) -> bool:
    """Both halves of the reproduction rule."""
    return bool(counts_match) and published_p_in_range(p_values, published)


def near_threshold_enrichment(
    cells: pd.DataFrame,
    labels: pd.DataFrame,
    thresholds: pd.DataFrame,
    *,
    fraction: float = 0.05,
) -> pd.DataFrame:
    """Fold enrichment of each cell type in each near-threshold band.

    For each donor and each of Pal's four thresholds, the band is the closest
    ``fraction`` of that donor's retained cells. Ties at the cutoff stay in the
    band. Bands and the retained background are pooled across donors. Genotype
    is not used. Enrichment is the type's share of the band divided by its share
    of all retained cells.
    """
    retained = _retained_with_types(cells, labels)
    threshold_frame = thresholds.set_index("donor")
    pieces = []
    for threshold, (metric, column) in _THRESHOLD_COLUMNS.items():
        if column not in threshold_frame.columns:
            raise ValueError(f"Donor table is missing {column}.")
        bands = []
        for donor, group in retained.groupby("donor", sort=True):
            distance = _threshold_distance(
                group[metric], float(threshold_frame.loc[donor, column]), metric
            )
            chosen = _closest_band(distance, fraction)
            bands.append(group.loc[chosen, ["donor", "cell_type"]])
        band = pd.concat(bands, ignore_index=True)
        pieces.append(_enrichment_table(threshold, band, retained))
    return pd.concat(pieces, ignore_index=True)


def largest_enrichment(enrichment: pd.DataFrame) -> dict:
    """The largest finite fold, and every type-threshold pair that reaches it."""
    finite = enrichment.replace([np.inf, -np.inf], np.nan).dropna(subset=["enrichment"])
    if finite.empty:
        raise ValueError("No finite enrichment was computed.")
    best = float(finite["enrichment"].max())
    winners = finite.loc[np.isclose(finite["enrichment"], best), ["threshold", "cell_type"]]
    return {
        "c": best,
        "winners": winners.drop_duplicates().reset_index(drop=True),
    }


def run_pal_reference(
    cells: pd.DataFrame,
    donors: pd.DataFrame,
    labels: pd.DataFrame,
    phase1: dict,
) -> dict:
    """Count check, variant p-values when the counts match, and the frozen enrichment."""
    prepared = prepare_labels(labels, donors)
    comparison = label_counts_match(prepared, cells, donors)
    published = float(phase1["pal_proportion_check"]["published_p"])
    variants = None
    if comparison["match"]:
        variants = variant_p_values(
            donor_cluster_counts(prepared), phase1["pal_proportion_check"]["variants"]
        )
        passed = reference_passes(
            counts_match=True, p_values=variants["p_value"], published=published
        )
    else:
        passed = False
    fraction = float(phase1["near_threshold_band"]["fraction"])
    enrichment = near_threshold_enrichment(
        cells, prepared, _threshold_table(donors), fraction=fraction
    )
    summary = largest_enrichment(enrichment)
    return {
        "labels": prepared,
        "comparison": comparison,
        "variants": variants,
        "passed": passed,
        "published_p": published,
        "enrichment": enrichment,
        "c": summary["c"],
        "winners": summary["winners"],
        "fraction": fraction,
        "cell_type_source": _cell_type_source(labels),
    }


def load_pal_labels(path: Path) -> pd.DataFrame:
    """Read the CSV written by ``scripts/export_pal_norm_b1_labels.R``."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {path}. Export Pal's Fig 4B object with "
            "Rscript scripts/export_pal_norm_b1_labels.R <SeuratObject_NormB1Total.rds> "
            f"{path}"
        )
    return pd.read_csv(path)


def prepare_labels(labels: pd.DataFrame, donors: pd.DataFrame) -> pd.DataFrame:
    """Normalize barcodes and attach donor ids and genotypes."""
    required = {"barcode", "sample", "cluster", "cell_type"}
    missing = required - set(labels.columns)
    if missing:
        raise ValueError(f"Pal export is missing columns {sorted(missing)}.")
    prepared = labels.copy()
    prepared["sample"] = [
        cell_sample_name(cell_name, sample)
        for cell_name, sample in zip(prepared["barcode"], prepared["sample"], strict=True)
    ]
    prepared["barcode"] = prepared["barcode"].map(normalize_barcode)
    prepared["donor"] = map_pal_samples(prepared["sample"], donors).to_numpy()
    donor_frame = _donor_frame(donors).set_index("donor")
    prepared["genotype"] = prepared["donor"].map(donor_frame["genotype"])
    if prepared["genotype"].isna().any():
        unknown = sorted(prepared.loc[prepared["genotype"].isna(), "donor"].unique())
        raise ValueError(f"No genotype for donors {unknown}.")
    prepared["cluster"] = prepared["cluster"].astype(str)
    prepared["cell_type"] = prepared["cell_type"].astype(str)
    return prepared


def write_pal_reference_outputs(
    result: dict,
    output_dir: Path | None = None,
) -> str:
    """Write the tables, freeze ``c``, and return the handoff.

    An existing Phase 2 section of the parameter file is left in place.
    """
    output_dir = output_directory() if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result["comparison"]["per_donor"].to_csv(output_dir / "02_count_check.csv", index=False)
    if result["variants"] is not None:
        result["variants"].to_csv(output_dir / VARIANTS_FILENAME, index=False)
    result["enrichment"].to_csv(output_dir / ENRICHMENT_FILENAME, index=False)
    counts = donor_cluster_counts(result["labels"])
    counts.to_csv(output_dir / COUNTS_FILENAME, index=False)
    _freeze_c(output_dir / PARAMS_FILENAME, result)
    handoff = pal_reference_handoff(result, output_dir)
    (output_dir / HANDOFF_FILENAME).write_text(handoff, encoding="utf-8")
    return handoff


def pal_reference_handoff(result: dict, output_dir: Path) -> str:
    """Markdown summary to paste back after the notebook runs."""
    comparison = result["comparison"]
    winners = result["winners"]
    winner_text = ", ".join(
        f"{row.cell_type} at {row.threshold}" for row in winners.itertuples(index=False)
    )
    lines = [
        "# Notebook 02 handoff",
        "",
        "Paste this file back into the coding session. This notebook reproduces Pal's",
        "retained-cell cluster result and freezes the near-threshold enrichment.",
        "It does not test whether removed cells are luminal progenitors.",
        "",
        f"Output directory: `{output_dir}`",
        "",
        "## Count match",
        "",
        "The check passes only when Pal's cells are exactly the retained cells:",
        "the same donor-barcode pairs, the same per-donor counts, and 59,766 cells.",
        "",
        f"- Match: {comparison['match']}",
        f"- Exported cells: {comparison['n_exported']}",
        f"- Retained cells: {comparison['n_retained']}",
        f"- Expected total: {comparison['expected_total']}",
        f"- Exported only: {comparison['only_exported']}",
        f"- Retained only: {comparison['only_retained']}",
        "",
        "## P = 0.14",
        "",
        "The model is a quasi-Poisson regression of donor-by-cluster counts on donor,",
        "cluster, and the cluster-by-genotype interaction. The p-value is for that",
        "interaction. Clusters are the raw Fig 4B ids. The published value has to lie",
        "inside the range of the four variants.",
        "",
    ]
    if result["variants"] is None:
        lines.append("The variants were not fit because the cell counts do not match.")
    else:
        shown = result["variants"].copy()
        shown["p_value"] = shown["p_value"].map(lambda value: f"{value:.4f}")
        lines.append(_markdown_table(shown))
        inside = published_p_in_range(result["variants"]["p_value"], result["published_p"])
        lines.extend(["", f"- {result['published_p']} inside the variant range: {inside}"])
    lines.extend(
        [
            "",
            f"- Reproduction passes: {result['passed']}",
            "",
            "## Near-threshold enrichment",
            "",
            "For each donor and each of Pal's four thresholds, the band is the closest",
            f"{result['fraction'] * 100:g}% of that donor's retained cells.",
            "Ties at the cutoff are kept.",
            "Enrichment pools the 12 donors and does not look at genotype. It is the",
            "cell type's share of the band divided by its share of all retained cells,",
            "a fold on the same scale as the tipping-point multiplier.",
            "",
            f"- c: {result['c']:.4f}",
            f"- Largest enrichment: {winner_text}",
            f"- Cell-type column source: {result['cell_type_source']}",
            "",
        ]
    )
    if result["cell_type_source"] == "ident" and _types_are_cluster_ids(result["labels"]):
        lines.extend(
            [
                "The export had no separate cell-type column. `c` uses the active identity,",
                "which matches the raw cluster id for every cell.",
                "",
            ]
        )
    return "\n".join(lines)


def _fit_poisson(counts: pd.DataFrame, *, interaction: bool):
    formula = "count ~ C(donor) + C(cluster)"
    if interaction:
        formula += " + C(cluster):C(genotype)"
    return smf.glm(formula, data=counts, family=Poisson()).fit(scale=1.0)


def _label_keys(labels: pd.DataFrame) -> pd.DataFrame:
    keys = pd.DataFrame(
        {
            "donor": labels["donor"].astype(str),
            "barcode": labels["barcode"].map(normalize_barcode),
        }
    )
    if keys.duplicated().any():
        raise ValueError("A donor-barcode pair is repeated in the Pal export.")
    return keys


def _retained_keys(cells: pd.DataFrame) -> pd.DataFrame:
    retained = cells.loc[cells["is_noise"].astype(int) == 0, ["donor", "barcode"]].copy()
    retained["barcode"] = retained["barcode"].map(normalize_barcode)
    retained["donor"] = retained["donor"].astype(str)
    if retained.duplicated().any():
        raise ValueError("A retained donor-barcode pair is repeated.")
    return retained.reset_index(drop=True)


def _retained_with_types(cells: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    retained = cells.loc[cells["is_noise"].astype(int) == 0].copy()
    retained["barcode"] = retained["barcode"].map(normalize_barcode)
    retained["donor"] = retained["donor"].astype(str)
    typed = labels[["donor", "barcode", "cell_type"]].copy()
    typed["barcode"] = typed["barcode"].map(normalize_barcode)
    typed["donor"] = typed["donor"].astype(str)
    merged = retained.merge(typed, on=["donor", "barcode"], how="left", validate="one_to_one")
    if merged["cell_type"].isna().any():
        raise ValueError("A retained cell has no Pal cell type.")
    return merged


def _threshold_distance(values: pd.Series, threshold: float, metric: str) -> np.ndarray:
    observed = values.to_numpy(dtype=float)
    if metric == "pct_counts_mt":
        observed = observed / 100.0
    return np.abs(observed - threshold)


def _closest_band(distance: np.ndarray, fraction: float) -> np.ndarray:
    n_cells = len(distance)
    n_band = max(1, int(round(fraction * n_cells)))
    order = np.argsort(distance, kind="mergesort")
    cutoff = distance[order[n_band - 1]]
    return distance <= cutoff


def _enrichment_table(threshold: str, band: pd.DataFrame, retained: pd.DataFrame) -> pd.DataFrame:
    band_counts = band["cell_type"].value_counts()
    retained_counts = retained["cell_type"].value_counts()
    n_band = len(band)
    n_retained = len(retained)
    rows = []
    for cell_type, retained_n in retained_counts.items():
        band_n = int(band_counts.get(cell_type, 0))
        share_band = band_n / n_band
        share_retained = retained_n / n_retained
        rows.append(
            {
                "threshold": threshold,
                "cell_type": cell_type,
                "n_band": band_n,
                "n_band_total": n_band,
                "n_retained": int(retained_n),
                "n_retained_total": n_retained,
                "share_band": share_band,
                "share_retained": share_retained,
                "enrichment": share_band / share_retained,
            }
        )
    return pd.DataFrame(rows)


def _threshold_table(donors: pd.DataFrame) -> pd.DataFrame:
    frame = _donor_frame(donors)
    columns = ["donor", *[column for _, column in _THRESHOLD_COLUMNS.values()]]
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"Donor table is missing {missing}.")
    return frame[columns]


def _donor_frame(donors: pd.DataFrame) -> pd.DataFrame:
    frame = donors.copy()
    if "specimen_id" not in frame.columns:
        frame = frame.reset_index(names="specimen_id")
    if "donor" not in frame.columns:
        raise ValueError("Donor table is missing a donor column.")
    return frame


def _library_token(specimen_id: str) -> str | None:
    match = _LIBRARY.search(str(specimen_id).upper())
    if match is None:
        return None
    return match.group(1)


def _normalize_name(value: str) -> str:
    text = str(value).strip().upper().replace("_", "-")
    if text.endswith("-TOTAL"):
        text = text[: -len("-TOTAL")]
    return text


def _cell_type_source(labels: pd.DataFrame) -> str:
    if "cell_type_source" not in labels.columns or labels["cell_type_source"].nunique() != 1:
        return "not recorded"
    return str(labels["cell_type_source"].iloc[0])


def _types_are_cluster_ids(labels: pd.DataFrame) -> bool:
    return bool((labels["cell_type"].astype(str) == labels["cluster"].astype(str)).all())


def _freeze_c(path: Path, result: dict) -> None:
    existing = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    phase1 = existing.setdefault("phase_1", {})
    band = phase1.setdefault("near_threshold_band", {})
    band["c"] = result["c"]
    band["c_winners"] = result["winners"].to_dict(orient="records")
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
