"""Tests for signals_in_the_noise.modeling.datasets."""

import numpy as np
import pandas as pd
import pytest

from signals_in_the_noise.analysis.noise_phenotypes import (
    ISCB_SCS_PBS_THRESHOLDS,
    PBS_LABELS,
    pbs_label,
    pbs_masks,
)
from signals_in_the_noise.modeling.datasets import (
    CORE_ARM_NAMES,
    META_COLUMNS,
    Arm,
    build_cell_table,
    feature_columns,
    pbs_features,
    qc_features,
    rank_qc_features,
    specimen_composition,
    specimen_sizes,
    standard_arms,
)
from signals_in_the_noise.preprocessing.qc import QC_METRICS


def test_pbs_features_one_hot_matches_pbs_label_and_drops_unlabelled(specimen_obs_factory):
    cells = specimen_obs_factory(noise_coupled=True).query("is_noise == 1")
    features = pbs_features(ISCB_SCS_PBS_THRESHOLDS)(cells)
    labels = pbs_label(pbs_masks(cells, ISCB_SCS_PBS_THRESHOLDS))

    assert list(features.columns) == list(PBS_LABELS)
    assert (features.sum(axis=1) == 1).all()
    assert features.index.equals(labels.index[labels != "none"])
    assert (features.idxmax(axis=1) == labels.loc[features.index]).all()


def test_pbs_features_is_invariant_to_monotone_shift_of_metrics(specimen_obs_factory):
    cells = specimen_obs_factory(noise_coupled=True)
    shifted = cells.assign(**{m: cells[m] * 2 + 5 for m in QC_METRICS})
    build = pbs_features(ISCB_SCS_PBS_THRESHOLDS)
    pd.testing.assert_frame_equal(build(cells), build(shifted))


def test_qc_and_rank_features(specimen_obs_factory):
    cells = specimen_obs_factory()
    assert list(qc_features(cells).columns) == list(QC_METRICS)
    ranks = rank_qc_features(cells)
    assert ranks.min().min() > 0 and ranks.max().max() == 1.0
    shifted = cells.assign(**{m: cells[m] + 100 for m in QC_METRICS})
    pd.testing.assert_frame_equal(rank_qc_features(shifted), ranks)


def _two_specimens(factory):
    specimen_obs = {
        "b": factory(n_noise=30, n_retained=80, seed=1),
        "a": factory(n_noise=20, n_retained=50, seed=2),
    }
    return specimen_obs, {"a": "Normal", "b": "ER+ tumour"}


def test_build_cell_table_layout_and_population(specimen_obs_factory):
    specimen_obs, conditions = _two_specimens(specimen_obs_factory)
    table = build_cell_table(specimen_obs, conditions, Arm("x", "retained", qc_features))

    assert list(table.columns) == ["specimen_id", "condition", "cell_id", *QC_METRICS]
    assert table["specimen_id"].unique().tolist() == ["a", "b"]
    assert specimen_sizes(table).to_dict() == {"a": 50, "b": 80}
    assert (table.loc[table["specimen_id"] == "b", "condition"] == "ER+ tumour").all()
    assert isinstance(table.index, pd.RangeIndex)
    assert feature_columns(table) == list(QC_METRICS)


def test_build_cell_table_adds_patient_column_when_given(specimen_obs_factory):
    specimen_obs, conditions = _two_specimens(specimen_obs_factory)
    table = build_cell_table(
        specimen_obs, conditions, Arm("x", "noise", qc_features), patients={"a": "p1", "b": "p1"}
    )
    assert list(table.columns) == [*META_COLUMNS, *QC_METRICS]
    assert set(table["patient_id"]) == {"p1"}
    assert feature_columns(table) == list(QC_METRICS)
    composition = specimen_composition(table)
    assert list(composition.columns[:2]) == ["patient_id", "condition"]


def test_build_cell_table_size_match_subsamples_to_noise_count(specimen_obs_factory):
    specimen_obs, conditions = _two_specimens(specimen_obs_factory)
    arm = Arm("x", "retained", qc_features, size_match=True)
    table = build_cell_table(specimen_obs, conditions, arm, seed=3)

    assert specimen_sizes(table).to_dict() == {"a": 20, "b": 30}
    retained_ids = set(specimen_obs["a"].query("is_noise == 0").index)
    assert set(table.loc[table["specimen_id"] == "a", "cell_id"]) <= retained_ids
    again = build_cell_table(specimen_obs, conditions, arm, seed=3)
    pd.testing.assert_frame_equal(table, again)


def test_build_cell_table_uses_alternative_noise_column(specimen_obs_factory):
    specimen_obs, conditions = _two_specimens(specimen_obs_factory)
    specimen_obs = {k: v.assign(is_noise_alt=1) for k, v in specimen_obs.items()}
    table = build_cell_table(
        specimen_obs, conditions, Arm("x", "noise", qc_features), noise_column="is_noise_alt"
    )
    assert specimen_sizes(table).to_dict() == {"a": 70, "b": 110}


def test_build_cell_table_empty_input():
    table = build_cell_table({}, {}, Arm("x", "all", qc_features))
    assert table.empty
    assert list(table.columns) == ["specimen_id", "condition", "cell_id"]


def test_specimen_composition_gives_subtype_fractions(specimen_obs_factory):
    specimen_obs, conditions = _two_specimens(specimen_obs_factory)
    table = build_cell_table(
        specimen_obs, conditions, Arm("x", "all", pbs_features(ISCB_SCS_PBS_THRESHOLDS))
    )
    composition = specimen_composition(table)

    assert composition.index.name == "specimen_id"
    assert composition.loc["b", "condition"] == "ER+ tumour"
    np.testing.assert_allclose(composition[list(PBS_LABELS)].sum(axis=1), 1.0)


def test_standard_arms_reference_is_published_noise_pbs():
    arms = standard_arms(ISCB_SCS_PBS_THRESHOLDS)
    names = [arm.name for arm in arms]
    assert names[0] == "noise_pbs"
    assert arms[0].population == "noise"
    assert len(set(names)) == len(names)
    assert {"retained_pbs", "retained_pbs_matched", "noise_qc_rank", "retained_qc"} <= set(names)


def test_standard_arms_subset_keeps_standard_order():
    arms = standard_arms(names=reversed(CORE_ARM_NAMES))
    assert [arm.name for arm in arms] == list(CORE_ARM_NAMES)
    with pytest.raises(ValueError, match="Unknown"):
        standard_arms(names=["nope"])


@pytest.mark.parametrize("population", ["noise", "retained", "all"])
def test_standard_arm_populations_are_valid(population):
    assert population in {arm.population for arm in standard_arms()}
