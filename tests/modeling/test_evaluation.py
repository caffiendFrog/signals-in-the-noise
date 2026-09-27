"""Tests for signals_in_the_noise.modeling.evaluation."""

import numpy as np
import pandas as pd
import pytest
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.svm import LinearSVC

from signals_in_the_noise.modeling.evaluation import (
    balanced_specimen_draws,
    evaluate_draws,
    evaluate_held_out,
    held_out_predictions,
    logistic_regression_factory,
    paired_comparison,
    permutation_test,
    random_forest_factory,
    summarize_predictions,
)


def _table(
    n_per_condition: int = 5, cells: int = 20, signal: float = 3.0, seed: int = 0
) -> pd.DataFrame:
    """Cell table where ``x`` is shifted by ``signal`` in ER specimens."""
    rng = np.random.default_rng(seed)
    rows = []
    for condition in ("ER", "N"):
        for i in range(n_per_condition):
            specimen = f"{condition}-{i}"
            x = rng.normal(signal if condition == "ER" else 0.0, 1.0, cells)
            rows.append(pd.DataFrame({"specimen_id": specimen, "condition": condition, "x": x}))
    return pd.concat(rows, ignore_index=True)


class _RecordingClassifier(ClassifierMixin, BaseEstimator):
    """Records the specimen codes (column 0) seen during each fit."""

    fits: list = []

    def fit(self, X, y):
        self.classes_ = np.unique(y)
        _RecordingClassifier.fits.append(set(X[:, 0]))
        return self

    def predict(self, X):
        return np.full(len(X), self.classes_[0])

    def predict_proba(self, X):
        return np.tile([1.0, 0.0], (len(X), 1))


def test_never_trains_on_the_held_out_specimen():
    table = _table()
    codes = {s: i for i, s in enumerate(table["specimen_id"].unique())}
    table["code"] = table["specimen_id"].map(codes).astype(float)
    _RecordingClassifier.fits = []

    predictions = held_out_predictions(
        table, ["code"], positive_label="ER", model_factory=_RecordingClassifier
    )

    assert len(_RecordingClassifier.fits) == len(codes)
    for specimen, code in codes.items():
        fits_without = [seen for seen in _RecordingClassifier.fits if code not in seen]
        assert len(fits_without) == 1
        assert fits_without[0] == set(codes.values()) - {code}
    assert predictions.index.equals(table.index)
    assert list(predictions.columns) == ["specimen_id", "y_true", "y_pred", "score"]


def test_patient_grouping_holds_out_all_specimens_of_a_patient():
    table = _table()
    patients = {"ER-0": "p0", "ER-1": "p0"}
    table["patient_id"] = table["specimen_id"].map(lambda s: patients.get(s, s))
    codes = {s: i for i, s in enumerate(table["specimen_id"].unique())}
    table["code"] = table["specimen_id"].map(codes).astype(float)
    _RecordingClassifier.fits = []

    predictions = held_out_predictions(
        table, ["code"], positive_label="ER", model_factory=_RecordingClassifier, group="patient_id"
    )

    assert len(_RecordingClassifier.fits) == 9
    assert any(
        {codes["ER-0"], codes["ER-1"]}.isdisjoint(seen) for seen in _RecordingClassifier.fits
    )
    assert not any(
        (codes["ER-0"] in seen) != (codes["ER-1"] in seen) for seen in _RecordingClassifier.fits
    )
    assert {"specimen_id", "patient_id"} <= set(predictions.columns)


def test_metrics_stay_per_specimen_for_mixed_condition_patients():
    """A patient with one ER and one Normal specimen must not collapse to one label."""
    table = _table()
    table["patient_id"] = table["specimen_id"].replace({"ER-0": "shared", "N-0": "shared"})
    metrics = evaluate_held_out(table, ["x"], positive_label="ER", group="patient_id")
    assert metrics["n_specimens"] == 10
    assert metrics["specimen_auc"] == 1.0
    result = permutation_test(
        table, ["x"], positive_label="ER", n_permutations=5, group="patient_id"
    )
    assert result.observed == 1.0


def test_informative_feature_separates_specimens():
    metrics = evaluate_held_out(_table(), ["x"], positive_label="ER")
    assert metrics["specimen_auc"] == 1.0
    assert metrics["specimen_balanced_accuracy"] == 1.0
    assert metrics["cell_accuracy"] > 0.9
    assert metrics["n_specimens"] == 10
    assert metrics["n_cells"] == 200


def test_uninformative_feature_is_not_perfect():
    metrics = evaluate_held_out(_table(signal=0.0), ["x"], positive_label="ER")
    assert metrics["specimen_auc"] < 0.9
    assert metrics["cell_accuracy"] < 0.65


@pytest.mark.parametrize(
    "factory", [random_forest_factory(), lambda: LinearSVC(class_weight="balanced")]
)
def test_score_is_oriented_towards_positive_label_for_other_models(factory):
    metrics = evaluate_held_out(_table(), ["x"], positive_label="ER", model_factory=factory)
    assert metrics["specimen_auc"] == 1.0


def test_score_orientation_when_positive_label_sorts_first():
    table = _table().replace({"condition": {"ER": "A-positive", "N": "B-negative"}})
    for factory in (logistic_regression_factory(), lambda: LinearSVC()):
        metrics = evaluate_held_out(
            table, ["x"], positive_label="A-positive", model_factory=factory
        )
        assert metrics["specimen_auc"] == 1.0


def test_unknown_positive_label_raises():
    with pytest.raises(ValueError, match="positive_label"):
        held_out_predictions(_table(), ["x"], positive_label="missing")


def test_residualizing_on_the_generating_covariate_removes_the_signal():
    table = _table(n_per_condition=8, signal=0.0)
    specimen_index = table["specimen_id"].str[-1].astype(float)
    table["depth"] = np.where(table["condition"] == "ER", 2.0, 0.0) + 0.1 * specimen_index
    table["x"] += 3 * table["depth"]
    assert evaluate_held_out(table, ["x"], positive_label="ER")["specimen_auc"] == 1.0
    adjusted = evaluate_held_out(table, ["x"], positive_label="ER", residualize_on=["depth"])
    assert adjusted["specimen_auc"] < 0.8


def test_residualizing_keeps_signal_independent_of_the_covariate():
    table = _table()
    table["depth"] = table["specimen_id"].str[-1].astype(float)
    adjusted = evaluate_held_out(table, ["x"], positive_label="ER", residualize_on=["depth"])
    assert adjusted["specimen_auc"] == 1.0


def test_group_may_be_the_index_for_specimen_tables():
    specimen_table = _table(cells=1).set_index("specimen_id")
    metrics = evaluate_held_out(specimen_table, ["x"], positive_label="ER")
    assert metrics["n_specimens"] == 10
    assert metrics["specimen_auc"] == 1.0


def test_summarize_majority_vote_ties_go_to_positive():
    predictions = pd.DataFrame(
        {
            "specimen_id": ["a", "a", "b", "b"],
            "y_true": ["ER", "ER", "N", "N"],
            "y_pred": ["ER", "N", "N", "N"],
            "score": [0.6, 0.4, 0.2, 0.3],
        }
    )
    metrics = summarize_predictions(predictions, positive_label="ER")
    assert metrics["specimen_balanced_accuracy"] == 1.0
    assert metrics["cell_accuracy"] == 0.75
    assert metrics["specimen_auc"] == 1.0


def test_summarize_auc_is_nan_with_one_class():
    predictions = pd.DataFrame(
        {
            "specimen_id": ["a", "b"],
            "y_true": ["ER", "ER"],
            "y_pred": ["ER", "N"],
            "score": [0.9, 0.1],
        }
    )
    assert np.isnan(summarize_predictions(predictions, positive_label="ER")["specimen_auc"])


def test_permutation_test_flags_real_signal():
    result = permutation_test(_table(), ["x"], positive_label="ER", n_permutations=40, seed=1)
    assert result.observed == 1.0
    assert result.null.shape == (40,)
    assert result.p_value <= 0.1
    assert result.null_mean < 0.8


def test_permutation_test_does_not_flag_null_data():
    result = permutation_test(
        _table(signal=0.0), ["x"], positive_label="ER", n_permutations=40, seed=1
    )
    assert result.p_value > 0.05


def test_permutation_test_is_deterministic_given_seed():
    a = permutation_test(_table(), ["x"], positive_label="ER", n_permutations=10, seed=3)
    b = permutation_test(_table(), ["x"], positive_label="ER", n_permutations=10, seed=3)
    np.testing.assert_array_equal(a.null, b.null)


def test_permutation_keeps_labels_constant_within_specimen():
    """Every permuted table must still assign one label per specimen, and labels must move."""
    fits = []

    class _Spy(_RecordingClassifier):
        def fit(self, X, y):
            fits.append(pd.DataFrame({"code": X[:, 0], "y": y}))
            return super().fit(X, y)

    table = _table()
    table["code"] = table["specimen_id"].factorize()[0].astype(float)
    permutation_test(table, ["code"], positive_label="ER", n_permutations=3, model_factory=_Spy)

    true_labels = table.groupby("code")["condition"].first()
    for fit in fits:
        per_specimen = fit.groupby("code")["y"]
        assert (per_specimen.nunique() == 1).all()
    permuted = [
        (fit.groupby("code")["y"].first() != true_labels[fit["code"].unique()]).any()
        for fit in fits
    ]
    assert any(permuted)


def _sizes_and_conditions():
    sizes = pd.Series(
        {"n1": 50, "n2": 60, "n3": 70, "e1": 40, "e2": 45, "e3": 500, "e4": 30, "e5": 90}
    )
    conditions = pd.Series({k: ("Normal" if k.startswith("n") else "ER") for k in sizes.index})
    return sizes, conditions


def test_balanced_draws_follow_published_protocol():
    sizes, conditions = _sizes_and_conditions()
    draws = balanced_specimen_draws(
        sizes, conditions, minority="Normal", majority="ER", n_minority=2, n_draws=25, seed=0
    )
    assert len(draws) == 25
    for draw in draws:
        minority = [s for s in draw if conditions[s] == "Normal"]
        majority = [s for s in draw if conditions[s] == "ER"]
        target = sizes[minority].sum()
        assert len(minority) == 2
        assert all(sizes[s] <= target for s in majority)
        eligible_total = sizes[
            [s for s in conditions.index[conditions == "ER"] if sizes[s] <= target]
        ].sum()
        assert sizes[majority].sum() >= min(target, eligible_total)
        assert all(isinstance(s, str) for s in draw)


def test_balanced_draws_cap_minority_and_are_reproducible():
    sizes, conditions = _sizes_and_conditions()
    kwargs = {"minority": "Normal", "majority": "ER", "n_minority": 10, "n_draws": 5}
    a = balanced_specimen_draws(sizes, conditions, seed=4, **kwargs)
    b = balanced_specimen_draws(sizes, conditions, seed=4, **kwargs)
    assert a == b
    assert all(sum(conditions[s] == "Normal" for s in draw) == 3 for draw in a)


def test_evaluate_draws_returns_one_row_per_draw():
    table = _table()
    draws = [["ER-0", "ER-1", "N-0", "N-1"], ["ER-2", "ER-3", "ER-4", "N-2", "N-3"]]
    frame = evaluate_draws(table, ["x"], draws, positive_label="ER")
    assert frame.index.name == "draw"
    assert frame["n_specimens"].tolist() == [4, 5]


def test_paired_comparison():
    a = pd.DataFrame({"m": [0.8, 0.7, 0.9]})
    b = pd.DataFrame({"m": [0.5, 0.75, 0.6]})
    result = paired_comparison(a, b, "m")
    assert result["mean_diff"] == pytest.approx(np.mean([0.3, -0.05, 0.3]))
    assert result["fraction_reference_greater"] == pytest.approx(2 / 3)
    assert result["mean_reference"] == pytest.approx(0.8)
    assert result["diff_ci_low"] <= result["mean_diff"] <= result["diff_ci_high"]
