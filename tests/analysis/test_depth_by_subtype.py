"""The per-type depth label follows the permutation p-values, and a one-specimen type is skipped."""

import pandas as pd

from signals_in_the_noise.analysis.depth_by_subtype import (
    _one_contrast,
    depth_reading,
    other_conditions,
)


def test_every_cancer_type_except_normal_is_a_contrast():
    conditions = {
        "a": "ER+ tumour",
        "b": "Normal",
        "c": "HER2+ tumour",
        "d": "Normal",
        "e": "PR+ tumour",
    }
    assert other_conditions(conditions) == ["ER+ tumour", "HER2+ tumour", "PR+ tumour"]


def test_pbs_separation_that_vanishes_after_thinning_is_labelled_depth():
    assert depth_reading(0.01, 0.40) == "collapses after thinning"
    assert depth_reading(0.01, 0.02) == "remains after thinning"
    assert depth_reading(0.20, 0.30) == "no specimen-level PBS separation"


def test_a_single_specimen_type_is_skipped_before_any_matrix_work():
    row = _one_contrast(
        {"a": None, "b": None, "c": None},
        {"a": "PR+ tumour", "b": "Normal", "c": "Normal"},
        {"a": "p1", "b": "p2", "c": "p3"},
        positive="PR+ tumour",
        reference="Normal",
        n_mads=3.0,
        n_permutations=10,
        seed=0,
        min_per_condition=2,
    )
    assert row["status"] == "skipped"
    assert row["n_positive"] == 1
    assert row["reading"] == "not compared"
    assert pd.isna(row["noise_pbs_p_original"])
