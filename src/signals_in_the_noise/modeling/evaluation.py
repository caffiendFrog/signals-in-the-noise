"""Grouped cross-validation with specimen-level inference.

Cells from one specimen are not independent, so every held-out fold is a whole
group (a specimen, or a patient when patients contribute several specimens),
headline metrics are computed per specimen, and significance comes from
permuting condition labels across specimens (never across cells).
"""

import logging
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.preprocessing import StandardScaler

from signals_in_the_noise.preprocessing.specimens import CONDITION_COLUMN, SPECIMEN_COLUMN

logger = logging.getLogger(__name__)

ModelFactory = Callable[[], object]
"""Zero-argument callable returning a fresh, unfitted sklearn classifier."""


def logistic_regression_factory(C: float = 10.0) -> ModelFactory:
    """Standardized, class-balanced logistic regression (the published primary model)."""

    def factory() -> Pipeline:
        return make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=1000, class_weight="balanced", C=C, random_state=43),
        )

    return factory


def random_forest_factory() -> ModelFactory:
    """Shallow class-balanced random forest with the published hyperparameters."""

    def factory() -> RandomForestClassifier:
        return RandomForestClassifier(
            max_depth=5, min_samples_leaf=20, class_weight="balanced", random_state=43
        )

    return factory


def _column_or_index(table: pd.DataFrame, name: str) -> np.ndarray:
    if name in table.columns:
        return table[name].to_numpy()
    return table.index.get_level_values(name).to_numpy()


def _positive_score(model, X: np.ndarray, positive_label: str) -> np.ndarray:
    classes = list(model.classes_)
    if positive_label not in classes:
        raise ValueError(f"positive_label {positive_label!r} not among trained classes {classes}.")
    if hasattr(model, "predict_proba"):
        return model.predict_proba(X)[:, classes.index(positive_label)]
    decision = model.decision_function(X)
    return decision if classes[1] == positive_label else -decision


def held_out_predictions(
    table: pd.DataFrame,
    feature_columns: Sequence[str],
    *,
    positive_label: str,
    model_factory: ModelFactory | None = None,
    target: str = CONDITION_COLUMN,
    group: str = SPECIMEN_COLUMN,
    residualize_on: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Predict every row with a model trained on all other groups (leave-one-group-out).

    Args:
        table: Cell- or specimen-level table; ``group`` may be a column or index level.
        feature_columns: Model inputs.
        positive_label: Class whose probability is returned as ``score``.
        model_factory: Builds a fresh model per fold. Defaults to
            :func:`logistic_regression_factory`.
        target: Label column.
        group: Fold unit; each fold holds out one group. Use the patient
            column when patients contribute several specimens.
        residualize_on: Optional covariate columns. Within each fold the
            features are replaced by their residuals from a linear regression
            on these covariates, fit on the training specimens only.

    Returns:
        DataFrame indexed like ``table`` with ``specimen_id`` (plus ``group``
        when it differs), ``y_true``, ``y_pred`` and ``score`` columns.
    """
    factory = model_factory or logistic_regression_factory()
    X = table[list(feature_columns)].to_numpy(dtype=float)
    y = table[target].to_numpy()
    groups = _column_or_index(table, group)
    covariates = None if not residualize_on else table[list(residualize_on)].to_numpy(dtype=float)

    y_pred = np.empty(len(table), dtype=object)
    score = np.empty(len(table), dtype=float)
    for train_idx, test_idx in LeaveOneGroupOut().split(X, y, groups):
        X_train, X_test = X[train_idx], X[test_idx]
        if covariates is not None:
            regression = LinearRegression().fit(covariates[train_idx], X_train)
            X_train = X_train - regression.predict(covariates[train_idx])
            X_test = X_test - regression.predict(covariates[test_idx])
        model = factory().fit(X_train, y[train_idx])
        y_pred[test_idx] = model.predict(X_test)
        score[test_idx] = _positive_score(model, X_test, positive_label)

    columns = {SPECIMEN_COLUMN: _column_or_index(table, SPECIMEN_COLUMN)}
    if group != SPECIMEN_COLUMN:
        columns[group] = groups
    columns.update({"y_true": y, "y_pred": y_pred, "score": score})
    return pd.DataFrame(columns, index=table.index.set_names([None] * table.index.nlevels))


def summarize_predictions(predictions: pd.DataFrame, *, positive_label: str) -> dict[str, float]:
    """Cell-level and specimen-level metrics for :func:`held_out_predictions` output.

    Specimen-level calls are a majority vote of cell predictions (ties go to
    ``positive_label``); ``specimen_auc`` ranks specimens by their mean score.
    Metrics are always per specimen, even when folds hold out whole patients.
    """
    y_true = predictions["y_true"].to_numpy()
    y_pred = predictions["y_pred"].to_numpy()
    per_specimen = predictions.groupby(SPECIMEN_COLUMN).agg(
        y_true=("y_true", "first"),
        positive_vote=("y_pred", lambda v: float(np.mean(v == positive_label))),
        score=("score", "mean"),
    )
    specimen_true = per_specimen["y_true"] == positive_label
    specimen_vote = per_specimen["positive_vote"] >= 0.5
    both_classes = specimen_true.nunique() == 2
    return {
        "n_specimens": float(len(per_specimen)),
        "n_cells": float(len(predictions)),
        "cell_accuracy": float(np.mean(y_true == y_pred)),
        "cell_balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "specimen_balanced_accuracy": float(balanced_accuracy_score(specimen_true, specimen_vote)),
        "specimen_auc": (
            float(roc_auc_score(specimen_true, per_specimen["score"])) if both_classes else np.nan
        ),
    }


def evaluate_held_out(
    table: pd.DataFrame, feature_columns: Sequence[str], *, positive_label: str, **cv_kwargs
) -> dict[str, float]:
    """:func:`held_out_predictions` followed by :func:`summarize_predictions`."""
    predictions = held_out_predictions(
        table, feature_columns, positive_label=positive_label, **cv_kwargs
    )
    return summarize_predictions(predictions, positive_label=positive_label)


@dataclass(frozen=True)
class PermutationResult:
    """Observed metric, its specimen-label permutation null, and the p-value."""

    metric: str
    observed: float
    null: np.ndarray
    p_value: float

    @property
    def null_mean(self) -> float:
        return float(np.nanmean(self.null))

    @property
    def null_q95(self) -> float:
        return float(np.nanquantile(self.null, 0.95))


def permutation_test(
    table: pd.DataFrame,
    feature_columns: Sequence[str],
    *,
    positive_label: str,
    n_permutations: int = 200,
    metric: str = "specimen_auc",
    seed: int = 0,
    **cv_kwargs,
) -> PermutationResult:
    """Permute condition labels across specimens and re-run the full held-out evaluation.

    All cells of a specimen keep a shared label, so the null respects the
    specimen as the unit of replication. Folds still follow ``group`` (e.g.
    patient) inside each permutation. ``p = (1 + #{null >= observed}) / (1 + n)``.
    """
    target = cv_kwargs.get("target", CONDITION_COLUMN)
    observed = evaluate_held_out(
        table, feature_columns, positive_label=positive_label, **cv_kwargs
    )[metric]

    units = _column_or_index(table, SPECIMEN_COLUMN)
    unit_labels = pd.Series(table[target].to_numpy(), index=units).groupby(level=0).first()
    rng = np.random.default_rng(seed)
    null = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        shuffled = dict(zip(unit_labels.index, rng.permutation(unit_labels.to_numpy())))
        permuted = table.assign(**{target: pd.Series(units, index=table.index).map(shuffled)})
        null[i] = evaluate_held_out(
            permuted, feature_columns, positive_label=positive_label, **cv_kwargs
        )[metric]

    p_value = float((1 + np.sum(null >= observed)) / (1 + n_permutations))
    return PermutationResult(metric=metric, observed=float(observed), null=null, p_value=p_value)


def balanced_specimen_draws(
    sizes: pd.Series,
    conditions: pd.Series,
    *,
    minority: str,
    majority: str,
    n_minority: int,
    n_draws: int,
    seed: int = 0,
) -> list[list[str]]:
    """Reproduce the published cell-balanced resampling of specimens.

    Each draw picks ``n_minority`` minority specimens, then adds majority
    specimens in random order (skipping any larger than the minority cell
    total) until their cell count reaches the minority cell total.

    Args:
        sizes: Cells per specimen, indexed by specimen identifier.
        conditions: Condition per specimen, same index as ``sizes``.
        minority: Condition sampled first.
        majority: Condition matched to the minority's cell count.
        n_minority: Minority specimens per draw (capped at the number available).
        n_draws: Number of draws.
        seed: Random seed.

    Returns:
        One sorted list of specimen identifiers per draw.
    """
    rng = np.random.default_rng(seed)
    minority_ids = sorted(conditions.index[conditions == minority])
    majority_ids = sorted(conditions.index[conditions == majority])
    k = min(n_minority, len(minority_ids))

    draws = []
    for _ in range(n_draws):
        chosen_minority = [str(s) for s in rng.choice(minority_ids, size=k, replace=False)]
        target_cells = int(sizes[chosen_minority].sum())
        chosen_majority, total = [], 0
        for specimen_id in rng.permutation(majority_ids):
            if sizes[specimen_id] > target_cells:
                continue
            chosen_majority.append(str(specimen_id))
            total += int(sizes[specimen_id])
            if total >= target_cells:
                break
        draws.append(sorted(chosen_minority) + sorted(chosen_majority))
    return draws


def evaluate_draws(
    table: pd.DataFrame,
    feature_columns: Sequence[str],
    draws: Sequence[Sequence[str]],
    *,
    positive_label: str,
    **cv_kwargs,
) -> pd.DataFrame:
    """Run :func:`evaluate_held_out` on each draw of specimen identifiers; one row per draw."""
    specimens = pd.Series(_column_or_index(table, SPECIMEN_COLUMN), index=table.index)
    rows = []
    for draw in draws:
        subset = table[specimens.isin(draw).to_numpy()]
        rows.append(
            evaluate_held_out(subset, feature_columns, positive_label=positive_label, **cv_kwargs)
        )
    frame = pd.DataFrame(rows)
    frame.index.name = "draw"
    return frame


def paired_comparison(
    reference: pd.DataFrame, other: pd.DataFrame, metric: str, *, ci: float = 0.95
) -> dict[str, float]:
    """Compare two arms evaluated on the same draws (rows aligned by draw index).

    Returns both means, the mean paired difference ``reference - other`` with a
    percentile interval over draws, and the fraction of draws where
    ``reference`` beats ``other``.
    """
    diff = (reference[metric] - other[metric]).dropna()
    tail = (1 - ci) / 2
    return {
        "mean_reference": float(reference[metric].mean()),
        "mean_other": float(other[metric].mean()),
        "mean_diff": float(diff.mean()),
        "diff_ci_low": float(diff.quantile(tail)),
        "diff_ci_high": float(diff.quantile(1 - tail)),
        "fraction_reference_greater": float((diff > 0).mean()),
    }
