"""What effect on Δ can the main test detect with 4 BRCA1 vs 8 WT donors?

Everything here is computed from WT donors only, so the framing of the main
test (hypothesis test or estimation) can be fixed and recorded before any
BRCA1 Δ is seen.

The main test compares mean Δ between genotypes with an exact permutation test.
That test depends on the data only through differences between donors, so for
normally distributed Δ its power depends on the effect only in units of the SD
of Δ; the SD estimated from WT donors therefore sets the minimum detectable
effect (MDE). The rescued share of epithelium bounds how large an effect can be
at all (``|Δ| <= w``, see
:func:`~signals_in_the_noise.analysis.rescue_arms.delta_by_donor`), which is
what makes an MDE judgeable as plausible or not.
"""

import itertools
import json
import logging
import math
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from signals_in_the_noise.analysis.statistics import (
    Alternative,
    bootstrap_mean_difference_ci,
    pooled_t_interval,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SdEstimate:
    """Sample SD with a chi-square confidence interval (assumes normal values)."""

    sd: float
    low: float
    high: float
    n: int
    confidence: float


def sd_with_ci(values: Sequence[float], *, confidence: float = 0.95) -> SdEstimate:
    """Sample SD (``ddof=1``) of ``values`` with its chi-square confidence interval."""
    x = np.asarray(values, dtype=float)
    if x.size < 2:
        raise ValueError("Need at least two values to estimate an SD.")
    df = x.size - 1
    sd = float(x.std(ddof=1))
    tail = (1.0 - confidence) / 2.0
    return SdEstimate(
        sd=sd,
        low=sd * math.sqrt(df / stats.chi2.ppf(1.0 - tail, df)),
        high=sd * math.sqrt(df / stats.chi2.ppf(tail, df)),
        n=int(x.size),
        confidence=confidence,
    )


@cache
def relabeling_masks(n_test: int, n_reference: int) -> np.ndarray:
    """Every assignment of ``n_test`` of the pooled donors to the test group.

    Row 0 is the observed labeling (the first ``n_test`` columns are the test group).
    """
    n_total = n_test + n_reference
    combinations = list(itertools.combinations(range(n_total), n_test))
    masks = np.zeros((len(combinations), n_total), dtype=bool)
    for row, members in enumerate(combinations):
        masks[row, list(members)] = True
    return masks


def permutation_p_values(
    data: np.ndarray, n_test: int, *, alternative: Alternative = "two-sided"
) -> np.ndarray:
    """Exact permutation p-values for the difference in means, for many datasets at once.

    Same test as :func:`~signals_in_the_noise.analysis.statistics.exact_permutation_test`
    with its default statistic, vectorised over datasets.

    Args:
        data: ``(n_datasets, n_donors)`` array; the first ``n_test`` columns are
            the test group.
        n_test: Size of the test group.
        alternative: ``"greater"``, ``"less"`` or ``"two-sided"``.
    """
    data = np.atleast_2d(np.asarray(data, dtype=float))
    n_reference = data.shape[1] - n_test
    masks = relabeling_masks(n_test, n_reference).astype(float)
    test_sums = data @ masks.T
    totals = data.sum(axis=1, keepdims=True)
    statistics = test_sums / n_test - (totals - test_sums) / n_reference
    observed = statistics[:, :1]
    tol = 1e-12 * np.maximum(1.0, np.abs(observed))
    if alternative == "greater":
        extreme = statistics >= observed - tol
    elif alternative == "less":
        extreme = statistics <= observed + tol
    elif alternative == "two-sided":
        extreme = np.abs(statistics) >= np.abs(observed) - tol
    else:
        raise ValueError(f"Unknown alternative {alternative!r}.")
    return extreme.mean(axis=1)


def normal_errors(
    n_simulations: int, n_donors: int, *, sd: float = 1.0, seed: int | None = None
) -> np.ndarray:
    """``(n_simulations, n_donors)`` independent normal donor-level errors."""
    return np.random.default_rng(seed).normal(0.0, sd, size=(n_simulations, n_donors))


def resampled_errors(
    residuals: Sequence[float], n_simulations: int, n_donors: int, *, seed: int | None = None
) -> np.ndarray:
    """Donor-level errors resampled from observed values instead of a normal distribution.

    ``residuals`` are centred and scaled by ``sqrt(n / (n - 1))`` so the
    resampled errors have the sample SD of ``residuals`` rather than the
    smaller plug-in SD.
    """
    r = np.asarray(residuals, dtype=float)
    r = (r - r.mean()) * math.sqrt(r.size / (r.size - 1))
    rng = np.random.default_rng(seed)
    return r[rng.integers(0, r.size, size=(n_simulations, n_donors))]


def simulate_power(
    effects: Sequence[float],
    errors: np.ndarray,
    n_test: int,
    *,
    alpha: float = 0.05,
    alternative: Alternative = "two-sided",
) -> np.ndarray:
    """Power of the exact permutation test for each effect on the test group.

    The same ``errors`` are reused for every effect (common random numbers),
    which keeps the power curve smooth.

    Args:
        effects: Shifts added to the test group (first ``n_test`` columns).
        errors: Output of :func:`normal_errors` or :func:`resampled_errors`.
        n_test: Size of the test group.
        alpha: Rejection threshold; a dataset counts as detected when ``p <= alpha``.
        alternative: Direction of the test.
    """
    power = np.empty(len(effects))
    for i, effect in enumerate(effects):
        data = errors.copy()
        data[:, :n_test] += effect
        power[i] = np.mean(permutation_p_values(data, n_test, alternative=alternative) <= alpha)
    return power


def minimum_detectable_effect(
    effects: Sequence[float], power: Sequence[float], *, target_power: float = 0.8
) -> float:
    """Smallest effect whose power reaches ``target_power``, interpolated linearly; NaN if never."""
    effects = np.asarray(effects, dtype=float)
    power = np.asarray(power, dtype=float)
    reached = np.flatnonzero(power >= target_power)
    if reached.size == 0:
        return float("nan")
    i = int(reached[0])
    if i == 0:
        return float(effects[0])
    return float(np.interp(target_power, power[i - 1 : i + 1], effects[i - 1 : i + 1]))


def t_interval_half_width(
    sd: float, n_test: int, n_reference: int, *, confidence: float = 0.95
) -> float:
    """Half-width of the pooled-variance t interval for a difference in means with this SD."""
    df = n_test + n_reference - 2
    return float(
        stats.t.ppf(0.5 + confidence / 2.0, df) * sd * math.sqrt(1.0 / n_test + 1.0 / n_reference)
    )


def simulate_ci_coverage(
    n_test: int,
    n_reference: int,
    *,
    confidence: float = 0.95,
    n_simulations: int = 2000,
    n_resamples: int = 2000,
    seed: int | None = None,
) -> pd.DataFrame:
    """Coverage and mean half-width (in SD units) of the donor bootstrap and pooled t intervals.

    Donor-level values are standard normal with no true difference; coverage
    does not depend on the true difference or the SD.
    """
    rng = np.random.default_rng(seed)
    methods = {
        "donor bootstrap (percentile)": lambda a, b, s: bootstrap_mean_difference_ci(
            a, b, confidence=confidence, n_resamples=n_resamples, seed=s
        ),
        "pooled-variance t": lambda a, b, s: pooled_t_interval(a, b, confidence=confidence),
    }
    covered = {name: 0 for name in methods}
    widths = {name: 0.0 for name in methods}
    for _ in range(n_simulations):
        a = rng.normal(size=n_test)
        b = rng.normal(size=n_reference)
        child_seed = int(rng.integers(2**32))
        for name, interval in methods.items():
            low, high = interval(a, b, child_seed)
            covered[name] += low <= 0.0 <= high
            widths[name] += (high - low) / 2.0
    return pd.DataFrame(
        {
            "method": list(methods),
            "coverage": [covered[m] / n_simulations for m in methods],
            "nominal": confidence,
            "mean_half_width_sd_units": [widths[m] / n_simulations for m in methods],
        }
    )


def ceiling_effect(
    rescued_share: Sequence[float],
    reference_lp_fraction: Sequence[float],
    delta: Sequence[float],
) -> float:
    """Effect on mean Δ if every rescued epithelial cell in a WT-like donor were LP.

    ``mean(w * (1 - p)) - mean(Δ)``: the largest effect the rescued cells can
    produce, using WT donors as stand-ins for BRCA1 donors' ``w`` and ``p``.
    """
    w = np.asarray(rescued_share, dtype=float)
    p = np.asarray(reference_lp_fraction, dtype=float)
    return float(np.mean(w * (1.0 - p)) - np.mean(np.asarray(delta, dtype=float)))


def fold_enrichment_effect(
    rescued_share: Sequence[float], rescued_lp_fraction: Sequence[float], fold: float
) -> float:
    """Effect on mean Δ if rescued epithelial cells were LP ``fold`` times as often as in WT.

    ``mean(w * (min(1, fold * r) - r))``. Donors without rescued epithelial
    cells (``w = 0``, ``r`` undefined) contribute zero.
    """
    w = np.asarray(rescued_share, dtype=float)
    r = np.nan_to_num(np.asarray(rescued_lp_fraction, dtype=float))
    return float(np.mean(w * (np.minimum(1.0, fold * r) - r)))


@dataclass(frozen=True)
class DesignDecision:
    """Framing of the main test, fixed from WT donors before BRCA1 Δ is computed.

    Attributes:
        n_test_donors: BRCA1 donors in the main test.
        n_reference_donors: WT donors the SD of Δ is estimated from.
        delta_sd: SD of Δ across WT donors.
        delta_sd_low: Lower confidence bound of ``delta_sd``.
        delta_sd_high: Upper confidence bound of ``delta_sd``.
        alpha: Significance level of the main test.
        alternative: Direction of the main test.
        target_power: Power defining the MDE.
        mde_sd_units: MDE in units of the SD of Δ (normal errors).
        mde: MDE on mean Δ (BRCA1 - WT) at ``delta_sd``.
        mde_low: MDE at ``delta_sd_low``.
        mde_high: MDE at ``delta_sd_high``.
        mde_resampled: MDE with errors resampled from the WT Δ values.
        ceiling_effect: Effect if every rescued epithelial cell were LP.
        plausible_effect: Effect if rescued cells were LP ``plausible_fold``
            times as often as in WT.
        plausible_fold: Enrichment defining ``plausible_effect``.
        exclusion_bound: Expected 95% CI half-width; a null result excludes
            effects larger than about this.
        identity_method: How cell identity was assigned for the WT Δ values.
    """

    n_test_donors: int
    n_reference_donors: int
    delta_sd: float
    delta_sd_low: float
    delta_sd_high: float
    alpha: float
    alternative: str
    target_power: float
    mde_sd_units: float
    mde: float
    mde_low: float
    mde_high: float
    mde_resampled: float
    ceiling_effect: float
    plausible_effect: float
    plausible_fold: float
    exclusion_bound: float
    identity_method: str

    @property
    def framing(self) -> str:
        """``"estimation"`` when the MDE exceeds the plausible effect, else ``"hypothesis_test"``.

        An undefined MDE (target power never reached) counts as exceeding it.
        """
        return "estimation" if not self.mde <= self.plausible_effect else "hypothesis_test"

    @property
    def ceiling_detectable(self) -> bool:
        """Whether even the largest possible effect reaches the target power."""
        return self.mde <= self.ceiling_effect

    def summary(self) -> str:
        """Human-readable verdict, meant to lead the report of the main test."""
        test = (
            f"{self.n_test_donors} vs {self.n_reference_donors} donors, "
            f"alpha={self.alpha:g} {self.alternative}, {self.target_power:.0%} power"
        )
        mde = (
            f"MDE {self.mde:.4f} on mean Δ "
            f"(range {self.mde_low:.4f}-{self.mde_high:.4f} over the SD CI; "
            f"{self.mde_sd_units:.2f} SD of Δ)"
        )
        effects = (
            f"plausible effect {self.plausible_effect:.4f} "
            f"(rescued cells LP at {self.plausible_fold:g}x the WT rate), "
            f"ceiling {self.ceiling_effect:.4f} (every rescued epithelial cell LP)"
        )
        if self.framing == "estimation":
            verdict = (
                "ESTIMATION. The design cannot detect a plausible effect"
                + ("" if self.ceiling_detectable else ", or even the largest possible one")
                + ". Report the BRCA1 - WT difference in mean Δ with its 95% CI, not a "
                f"verdict from the p-value; a null result excludes effects larger than "
                f"about {self.exclusion_bound:.4f}."
            )
        else:
            verdict = "HYPOTHESIS TEST. The design can detect a plausible effect."
        return f"{verdict}\n{test}: {mde}; {effects}."

    def to_json(self, path: Path) -> None:
        """Persist the decision so Step 5 reads it instead of recomputing it after unblinding."""
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            **asdict(self),
            "framing": self.framing,
            "ceiling_detectable": self.ceiling_detectable,
        }
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    @classmethod
    def from_json(cls, path: Path) -> "DesignDecision":
        """Load a decision written by :meth:`to_json`."""
        payload = json.loads(path.read_text(encoding="utf-8"))
        return cls(**{name: payload[name] for name in cls.__dataclass_fields__})
