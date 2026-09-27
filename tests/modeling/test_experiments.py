"""Tests for signals_in_the_noise.modeling.experiments.

The synthetic worlds check that the arm suite can tell apart the two
explanations it exists to separate: a signal specific to discarded cells, and
a sample-level technical difference that touches every cell.
"""

import numpy as np
import pandas as pd
import pytest

from signals_in_the_noise.analysis.noise_phenotypes import ISCB_SCS_PBS_THRESHOLDS
from signals_in_the_noise.modeling.datasets import specimen_composition, standard_arms
from signals_in_the_noise.modeling.experiments import (
    InsufficientSpecimensError,
    compare_feature_sets,
    reference_specimens,
    run_arm_suite,
)
from signals_in_the_noise.preprocessing.qc import (
    QcThresholds,
    flag_noise,
    modal_thresholds,
    relabel_noise,
)
from tests.conftest import make_specimen_obs, make_world

POSITIVE, NEGATIVE = "ER+ tumour", "Normal"


def _run(world_factory, kind: str, **kwargs):
    specimen_obs, conditions = world_factory(kind)
    defaults = {"min_reference_cells": 5, "n_draws": 3, "n_permutations": 5, "n_minority": 4}
    defaults.update(kwargs)
    return run_arm_suite(
        specimen_obs,
        conditions,
        standard_arms(ISCB_SCS_PBS_THRESHOLDS),
        positive_label=POSITIVE,
        negative_label=NEGATIVE,
        **defaults,
    )


@pytest.fixture(scope="module")
def specific_result():
    return _run(make_world, "specific")


@pytest.fixture(scope="module")
def confounded_result():
    return _run(make_world, "confounded")


def test_run_arm_suite_raises_when_too_few_specimens_pass_minimum(world_factory):
    with pytest.raises(InsufficientSpecimensError, match="min|>="):
        _run(world_factory, "null", min_reference_cells=10**6)


def test_specific_world_signal_lives_only_in_discarded_cells(specific_result):
    auc = specific_result.full["specimen_auc"]
    assert auc["noise_pbs"] == 1.0
    assert auc["retained_pbs"] <= 0.75
    assert auc["retained_pbs_matched"] <= 0.75
    assert auc["retained_qc"] <= 0.75


def test_confounded_world_signal_also_appears_in_retained_cells(confounded_result):
    auc = confounded_result.full["specimen_auc"]
    assert auc["noise_pbs"] == 1.0
    assert auc["retained_pbs"] == 1.0
    assert auc["retained_pbs_matched"] == 1.0
    assert auc["retained_qc"] == 1.0


def test_rank_features_ignore_pure_location_shift(world_factory):
    """A shift with no change in within-specimen structure is invisible to ranks and PBS."""
    specimen_obs, conditions = world_factory("null")
    for sid, obs in specimen_obs.items():
        if conditions[sid] == POSITIVE:
            obs["log1p_total_counts"] += 2.0
    result = run_arm_suite(
        specimen_obs,
        conditions,
        standard_arms(ISCB_SCS_PBS_THRESHOLDS),
        positive_label=POSITIVE,
        negative_label=NEGATIVE,
        min_reference_cells=5,
        n_draws=2,
        n_permutations=3,
        n_minority=4,
    )
    auc = result.full["specimen_auc"]
    assert auc["retained_qc"] == 1.0
    assert auc["noise_qc"] == 1.0
    assert auc["noise_qc_rank"] <= 0.75
    assert auc["noise_pbs"] <= 0.75


def _threshold_world() -> tuple[dict[str, pd.DataFrame], dict[str, str], dict[str, QcThresholds]]:
    """Identical cells in both conditions; only the per-specimen gene threshold differs.

    ER+ specimens get a lowered gene threshold (Pal's low-coverage adjustment),
    so fewer low-gene cells are discarded there.
    """
    specimen_obs, conditions, thresholds = {}, {}, {}
    for condition, genes_lower in ((POSITIVE, np.expm1(6 - 0.4 * 1.6)), (NEGATIVE, np.expm1(5.84))):
        for i in range(6):
            sid = f"{condition[:2]}-{i}"
            obs = make_specimen_obs(0, 1500, seed=len(specimen_obs)).drop(columns="is_noise")
            thresholds[sid] = QcThresholds(
                mito_upper=float(1 / (1 + np.exp(1.0))),
                genes_lower=float(genes_lower),
                genes_upper=np.inf,
                total_upper=np.inf,
            )
            obs["is_noise"] = flag_noise(obs, thresholds[sid])["is_noise"]
            specimen_obs[sid], conditions[sid] = obs, condition
    return specimen_obs, conditions, thresholds


def test_threshold_driven_noise_membership_mimics_a_discarded_cell_signal():
    """Per-sample thresholds alone make noise_pbs look specific; a uniform rule removes it."""
    specimen_obs, conditions, thresholds = _threshold_world()
    uniform = modal_thresholds(thresholds.values())
    specimen_obs = relabel_noise(specimen_obs, lambda obs: uniform, "is_noise_uniform")
    arms = standard_arms(ISCB_SCS_PBS_THRESHOLDS, names=["noise_pbs", "retained_pbs_matched"])
    kwargs = {"min_reference_cells": 5, "n_draws": 1, "n_permutations": 1, "n_minority": 4}

    published = run_arm_suite(
        specimen_obs, conditions, arms, positive_label=POSITIVE, negative_label=NEGATIVE, **kwargs
    ).full["specimen_auc"]
    assert published["noise_pbs"] == 1.0
    assert published["retained_pbs_matched"] <= 0.75

    relabelled = run_arm_suite(
        specimen_obs,
        conditions,
        arms,
        positive_label=POSITIVE,
        negative_label=NEGATIVE,
        noise_column="is_noise_uniform",
        **kwargs,
    ).full["specimen_auc"]
    assert relabelled["noise_pbs"] <= 0.75


def test_suite_shares_specimens_and_draws_across_arms(specific_result):
    specimens = set(specific_result.specimens)
    for table in specific_result.tables.values():
        assert set(table["specimen_id"]) == specimens
    draw_counts = specific_result.draw_metrics.groupby("arm")["draw"].nunique()
    assert (draw_counts == len(specific_result.draws)).all()
    assert specific_result.reference == "noise_pbs"


def test_summary_and_paired_tables(specific_result):
    summary = specific_result.summary()
    assert list(summary.index) == [arm.name for arm in specific_result.arms]
    assert {"specimen_auc", "permutation_p", "draws_mean_cell_accuracy", "description"} <= set(
        summary.columns
    )

    paired = specific_result.paired_vs_reference("cell_accuracy")
    assert "noise_pbs" not in paired.index
    assert paired.loc["retained_pbs", "mean_diff"] > 0
    assert paired.loc["retained_pbs", "mean_reference"] == pytest.approx(
        specific_result.draw_metrics.query("arm == 'noise_pbs'")["cell_accuracy"].mean()
    )


def test_suite_with_patients_holds_out_patients_but_scores_specimens(world_factory):
    specimen_obs, conditions = world_factory("specific")
    patients = {sid: sid for sid in specimen_obs}
    patients.update({"ER-00": "shared", "ER-01": "shared", "N-00": "mixed", "ER-02": "mixed"})
    result = run_arm_suite(
        specimen_obs,
        conditions,
        standard_arms(ISCB_SCS_PBS_THRESHOLDS)[:2],
        positive_label=POSITIVE,
        negative_label=NEGATIVE,
        patients=patients,
        min_reference_cells=5,
        n_draws=2,
        n_permutations=3,
        n_minority=4,
    )
    assert "patient_id" in result.tables["noise_pbs"].columns
    assert result.full.loc["noise_pbs", "n_specimens"] == 12
    assert result.full.loc["noise_pbs", "specimen_auc"] == 1.0


def test_reference_specimens_keeps_specimens_at_or_above_minimum():
    table = pd.DataFrame({"specimen_id": ["a"] * 3 + ["b"] * 5 + ["c"] * 1})
    assert reference_specimens(table, 3) == ["a", "b"]


def test_min_reference_cells_filters_specimens(world_factory):
    specimen_obs, conditions = world_factory("specific")
    specimen_obs["N-tiny"] = specimen_obs["N-00"].iloc[:3].copy()
    conditions["N-tiny"] = NEGATIVE
    result = run_arm_suite(
        specimen_obs,
        conditions,
        standard_arms(ISCB_SCS_PBS_THRESHOLDS)[:2],
        positive_label=POSITIVE,
        negative_label=NEGATIVE,
        min_reference_cells=5,
        n_draws=1,
        n_permutations=1,
        n_minority=4,
    )
    assert "N-tiny" not in result.specimens


def test_suite_ignores_other_conditions_and_requires_both(world_factory):
    specimen_obs, conditions = world_factory("specific")
    conditions = {k: ("HER2+ tumour" if v == NEGATIVE else v) for k, v in conditions.items()}
    with pytest.raises(InsufficientSpecimensError, match="ER"):
        run_arm_suite(
            specimen_obs,
            conditions,
            standard_arms()[:1],
            positive_label=POSITIVE,
            negative_label=NEGATIVE,
            min_reference_cells=1,
        )


def test_compare_feature_sets_detects_confounded_composition(world_factory):
    """PBS composition driven entirely by a technical covariate should not survive adjustment."""
    specimen_obs, conditions = world_factory("null", n_per_condition=8)
    rng = np.random.default_rng(0)
    table = pd.DataFrame(
        {
            "condition": [conditions[s] for s in specimen_obs],
            "depth": [
                (1.0 if conditions[s] == POSITIVE else 0.0) + rng.normal(0, 0.3)
                for s in specimen_obs
            ],
        },
        index=pd.Index(list(specimen_obs), name="specimen_id"),
    )
    table["pbs-1"] = 0.2 * table["depth"] + rng.normal(0, 0.05, len(table))
    result = compare_feature_sets(
        table,
        {"pbs": (["pbs-1"], None), "pbs_adjusted": (["pbs-1"], ["depth"])},
        positive_label=POSITIVE,
        n_permutations=30,
    )
    assert result.loc["pbs", "specimen_auc"] > 0.85
    assert result.loc["pbs", "permutation_p"] < 0.1
    assert result.loc["pbs_adjusted", "specimen_auc"] < 0.8
    assert result.loc["pbs_adjusted", "adjusted_for"] == "depth"


def test_specimen_composition_feeds_compare_feature_sets(specific_result):
    composition = specimen_composition(specific_result.tables["noise_pbs"])
    result = compare_feature_sets(
        composition,
        {"pbs": (["pbs-1", "pbs-2", "pbs-3"], None)},
        positive_label=POSITIVE,
        n_permutations=5,
    )
    assert result.loc["pbs", "specimen_auc"] == 1.0
