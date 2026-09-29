"""Run a set of arms under one shared protocol so their results are directly comparable."""

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import pandas as pd

from signals_in_the_noise.modeling.datasets import (
    Arm,
    build_cell_table,
    feature_columns,
    specimen_sizes,
)
from signals_in_the_noise.modeling.evaluation import (
    ModelFactory,
    PermutationResult,
    balanced_specimen_draws,
    evaluate_draws,
    evaluate_held_out,
    paired_comparison,
    permutation_test,
)
from signals_in_the_noise.preprocessing.specimens import (
    CONDITION_COLUMN,
    PATIENT_COLUMN,
    SPECIMEN_COLUMN,
)

logger = logging.getLogger(__name__)

PUBLISHED_PROTOCOL: dict[str, int] = {"min_reference_cells": 50, "n_minority": 6, "n_draws": 100}
"""Settings of the published ``iscb-scs/10-a`` classifier: specimens need at
least 50 PBS-labelled noise cells, and each of 100 draws uses 6 Normal
specimens matched by cell count with ER+ specimens."""


class InsufficientSpecimensError(ValueError):
    """Too few specimens pass the reference-arm cell minimum to compare conditions."""


def reference_specimens(reference_table: pd.DataFrame, min_cells: int) -> list[str]:
    """Specimens contributing at least ``min_cells`` rows to the reference arm's table."""
    sizes = specimen_sizes(reference_table)
    return sorted(sizes.index[sizes >= min_cells])


@dataclass
class ArmSuiteResult:
    """Outputs of :func:`run_arm_suite`.

    Attributes:
        arms: Arms in run order; the first is the reference.
        specimens: Specimens included in every arm.
        tables: Cell table per arm, restricted to ``specimens``.
        full: Per-arm held-out metrics over all included specimens.
        permutations: Per-arm specimen-label permutation result.
        draws: Shared balanced specimen draws.
        draw_metrics: Long table of per-arm, per-draw metrics.
    """

    arms: list[Arm]
    specimens: list[str]
    tables: dict[str, pd.DataFrame]
    full: pd.DataFrame
    permutations: dict[str, PermutationResult]
    draws: list[list[str]]
    draw_metrics: pd.DataFrame = field(repr=False)

    @property
    def reference(self) -> str:
        return self.arms[0].name

    def summary(self) -> pd.DataFrame:
        """One row per arm: full-cohort metrics, permutation p-value and mean draw metrics."""
        permutation = pd.DataFrame(
            {
                name: {
                    "permutation_metric": result.metric,
                    "permutation_null_mean": result.null_mean,
                    "permutation_null_q95": result.null_q95,
                    "permutation_p": result.p_value,
                }
                for name, result in self.permutations.items()
            }
        ).T
        draws = self.draw_metrics.groupby("arm").mean(numeric_only=True).add_prefix("draws_mean_")
        descriptions = pd.Series(
            {arm.name: arm.description for arm in self.arms}, name="description"
        )
        summary = self.full.join(permutation).join(draws).join(descriptions)
        return summary.loc[[arm.name for arm in self.arms]]

    def paired_vs_reference(self, metric: str = "cell_accuracy") -> pd.DataFrame:
        """Paired draw-by-draw comparison of the reference arm against each other arm.

        ``mean_diff`` is reference minus arm, so positive values favour the reference.
        """
        by_arm = {name: frame.set_index("draw") for name, frame in self.draw_metrics.groupby("arm")}
        reference = by_arm[self.reference]
        rows = {
            arm.name: paired_comparison(reference, by_arm[arm.name], metric)
            for arm in self.arms[1:]
        }
        return pd.DataFrame.from_dict(rows, orient="index")


def run_arm_suite(
    specimen_obs: Mapping[str, pd.DataFrame],
    conditions: Mapping[str, str],
    arms: Sequence[Arm],
    *,
    positive_label: str,
    negative_label: str,
    patients: Mapping[str, str] | None = None,
    noise_column: str = "is_noise",
    min_reference_cells: int = PUBLISHED_PROTOCOL["min_reference_cells"],
    n_minority: int = PUBLISHED_PROTOCOL["n_minority"],
    n_draws: int = PUBLISHED_PROTOCOL["n_draws"],
    n_permutations: int = 200,
    model_factory: ModelFactory | None = None,
    seed: int = 0,
) -> ArmSuiteResult:
    """Evaluate every arm on the same specimens, draws and permutation seeds.

    Specimens are included when the reference arm (``arms[0]``) yields at
    least ``min_reference_cells`` cells for them, matching the published
    filter; that set is then fixed for every arm so differences between arms
    cannot come from different specimens. Balanced draws are computed once
    from reference-arm cell counts and shared, making per-draw comparisons paired.

    When ``patients`` is given, each fold holds out a whole patient, while
    metrics and label permutations stay per specimen.
    """
    if not arms:
        raise ValueError("At least one arm is required.")
    wanted = {positive_label, negative_label}
    specimen_obs = {sid: obs for sid, obs in specimen_obs.items() if conditions[sid] in wanted}

    tables = {
        arm.name: build_cell_table(
            specimen_obs, conditions, arm, patients=patients, noise_column=noise_column, seed=seed
        )
        for arm in arms
    }
    reference_sizes = specimen_sizes(tables[arms[0].name])
    specimens = reference_specimens(tables[arms[0].name], min_reference_cells)
    tables = {
        name: t[t[SPECIMEN_COLUMN].isin(specimens)].reset_index(drop=True)
        for name, t in tables.items()
    }
    included_conditions = pd.Series({sid: conditions[sid] for sid in specimens})
    counts = included_conditions.value_counts()
    if len(counts) < 2:
        raise InsufficientSpecimensError(
            f"Need both conditions among specimens with >= {min_reference_cells} "
            f"'{arms[0].name}' cells; got {counts.to_dict()}."
        )
    logger.info("included specimens per condition: %s", counts.to_dict())

    minority = (
        negative_label if counts[negative_label] <= counts[positive_label] else positive_label
    )
    majority = positive_label if minority == negative_label else negative_label
    draws = balanced_specimen_draws(
        reference_sizes[specimens],
        included_conditions,
        minority=minority,
        majority=majority,
        n_minority=n_minority,
        n_draws=n_draws,
        seed=seed,
    )

    cv_kwargs = {
        "positive_label": positive_label,
        "model_factory": model_factory,
        "group": SPECIMEN_COLUMN if patients is None else PATIENT_COLUMN,
    }
    full, permutations, draw_frames = {}, {}, []
    for arm in arms:
        logger.info("evaluating arm %s", arm.name)
        table = tables[arm.name]
        columns = feature_columns(table)
        full[arm.name] = evaluate_held_out(table, columns, **cv_kwargs)
        permutations[arm.name] = permutation_test(
            table, columns, n_permutations=n_permutations, seed=seed, **cv_kwargs
        )
        per_draw = evaluate_draws(table, columns, draws, **cv_kwargs).reset_index()
        per_draw.insert(0, "arm", arm.name)
        draw_frames.append(per_draw)

    return ArmSuiteResult(
        arms=list(arms),
        specimens=specimens,
        tables=tables,
        full=pd.DataFrame.from_dict(full, orient="index"),
        permutations=permutations,
        draws=draws,
        draw_metrics=pd.concat(draw_frames, ignore_index=True),
    )


def compare_feature_sets(
    table: pd.DataFrame,
    feature_sets: Mapping[str, tuple[Sequence[str], Sequence[str] | None]],
    *,
    positive_label: str,
    n_permutations: int = 1000,
    model_factory: ModelFactory | None = None,
    group: str = SPECIMEN_COLUMN,
    seed: int = 0,
) -> pd.DataFrame:
    """Held-out evaluation + permutation for several feature sets on one specimen-level table.

    Args:
        table: One row per specimen with a ``condition`` column; the specimen
            identifier may be the index.
        feature_sets: Name to ``(feature_columns, residualize_on)``; set
            ``residualize_on`` to ``None`` for no adjustment.
        positive_label: Positive condition.
        n_permutations: Specimen-label permutations per feature set.
        model_factory: Model builder; defaults to logistic regression.
        group: Fold unit, e.g. ``"patient_id"`` to hold out whole patients.
        seed: Permutation seed, shared across feature sets.

    Returns:
        One row per feature set with specimen metrics and permutation p-value.
    """
    rows = {}
    for name, (columns, residualize_on) in feature_sets.items():
        kwargs = {
            "positive_label": positive_label,
            "model_factory": model_factory,
            "residualize_on": residualize_on,
            "target": CONDITION_COLUMN,
            "group": group,
        }
        metrics = evaluate_held_out(table, columns, **kwargs)
        result = permutation_test(
            table, columns, n_permutations=n_permutations, seed=seed, **kwargs
        )
        rows[name] = {
            "features": ", ".join(columns),
            "adjusted_for": ", ".join(residualize_on or []),
            "specimen_auc": metrics["specimen_auc"],
            "specimen_balanced_accuracy": metrics["specimen_balanced_accuracy"],
            "permutation_null_mean": result.null_mean,
            "permutation_p": result.p_value,
        }
    return pd.DataFrame.from_dict(rows, orient="index")
