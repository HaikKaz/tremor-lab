"""The Aki b-value: exact arithmetic, then recovery of a known b from synthetic
draws."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tremor_lab import constants
from tremor_lab.bvalue import (
    b_stability,
    b_value_aki,
    b_value_tinti,
    mc_b_stability,
)

DATA = Path(__file__).parent / "data"
LN10 = np.log(10)


def gutenberg_richter_sample(b, mc, n, seed, dm=0.1):
    """Magnitudes complete from the lower edge of the Mc bin, reported on the grid.

    The exceedance law 10 ** (-b (M - Mc)) is an exponential of rate b ln10. A
    catalogue complete at Mc contains events from the lower edge of that bin upwards,
    which is what the half-bin offset in the estimator assumes.
    """
    rng = np.random.default_rng(seed)
    draws = mc - dm / 2 + rng.exponential(1.0 / (b * LN10), n)
    return np.round(draws, 1) if dm else draws


def test_b_value_of_a_single_magnitude_bin():
    # Every event at Mc exactly: mean - (Mc - dm/2) = 0.05, so b = 1 / (ln10 * 0.05).
    result = b_value_aki([3.0, 3.0, 3.0, 3.0], mc=3.0)
    assert result.b == pytest.approx(1.0 / (LN10 * 0.05))
    assert result.sigma == pytest.approx(0.0)
    assert result.n == 4


def test_events_below_the_completeness_bin_are_discarded():
    result = b_value_aki([2.0, 2.9, 2.95, 3.0, 3.1], mc=3.0)
    assert result.n == 3


def test_the_half_bin_offset_can_be_switched_off():
    unbinned = b_value_aki([3.0, 3.0, 3.0], mc=2.95, dm=0.0)
    assert unbinned.b == pytest.approx(1.0 / (LN10 * 0.05))


def test_too_few_events_is_an_error_not_a_silent_nan():
    with pytest.raises(ValueError, match="at least two events"):
        b_value_aki([3.0], mc=3.0)


def test_standard_error_coefficient_can_be_overridden(monkeypatch):
    mags = gutenberg_richter_sample(b=1.0, mc=3.0, n=2000, seed=0)
    doubled = b_value_aki(mags, mc=3.0, shi_bolt_k=4.60)
    assert doubled.sigma == pytest.approx(2 * b_value_aki(mags, mc=3.0).sigma)
    monkeypatch.setattr(constants, "SHI_BOLT_K", 4.60)
    assert b_value_aki(mags, mc=3.0).sigma == pytest.approx(doubled.sigma)


def test_standard_error_shrinks_as_the_catalogue_grows():
    small = b_value_aki(gutenberg_richter_sample(1.0, 3.0, 1000, 0), mc=3.0)
    large = b_value_aki(gutenberg_richter_sample(1.0, 3.0, 100000, 0), mc=3.0)
    assert large.sigma < small.sigma
    # Shi and Bolt's sigma carries the 1/sqrt(n) of a standard error.
    assert large.sigma == pytest.approx(small.sigma / 10, rel=0.15)


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("b_true", [0.7, 0.85, 1.0, 1.2])
def test_a_known_b_is_recovered_from_binned_magnitudes(b_true, seed):
    mags = gutenberg_richter_sample(b_true, mc=3.0, n=20000, seed=seed)
    assert b_value_aki(mags, mc=3.0).b == pytest.approx(b_true, abs=0.02)


@pytest.mark.parametrize("seed", range(5))
def test_a_known_b_is_recovered_from_unbinned_magnitudes(seed):
    mags = gutenberg_richter_sample(1.0, mc=3.0, n=20000, seed=seed, dm=0.0)
    assert b_value_aki(mags, mc=3.0, dm=0.0).b == pytest.approx(1.0, abs=0.02)


def test_the_true_b_lies_within_a_few_standard_errors():
    result = b_value_aki(gutenberg_richter_sample(0.85, 3.0, 20000, 3), mc=3.0)
    assert abs(result.b - 0.85) < 3 * result.sigma


def test_the_stability_curve_is_flat_for_a_single_gutenberg_richter_population():
    # Drawn from one law, so b must not drift with the threshold; this is the
    # control that gives the curve its meaning on real data.
    mags = gutenberg_richter_sample(b=1.0, mc=3.0, n=200000, seed=0)
    curve = b_stability(mags, dm=0.1, min_events=500)
    assert curve.b.size > 10
    assert curve.b.max() - curve.b.min() < 0.08
    assert curve.b.mean() == pytest.approx(1.0, abs=0.02)


def test_the_stability_curve_reports_one_entry_per_usable_threshold():
    mags = gutenberg_richter_sample(b=1.0, mc=3.0, n=5000, seed=1)
    curve = b_stability(mags, thresholds=[3.0, 3.5, 4.0], dm=0.1, min_events=1)
    assert curve.thresholds.tolist() == [3.0, 3.5, 4.0]
    assert (np.diff(curve.n) < 0).all()


def test_b_stability_finds_completeness_for_a_catalogue_missing_small_events():
    # Complete from 3.0, then everything below 3.5 thinned away: the slope only
    # settles above the true completeness.
    rng = np.random.default_rng(3)
    mags = gutenberg_richter_sample(b=1.0, mc=3.0, n=200000, seed=2)
    keep = (mags >= 3.5) | (rng.random(mags.size) < 0.25)
    mc = mc_b_stability(mags[keep], dm=0.1, min_events=200)
    assert mc is not None
    assert 3.3 <= mc <= 3.8


def test_one_anomalously_small_event_does_not_move_the_stability_grid():
    # The grid is anchored to the mode of the incremental distribution, not to
    # the smallest magnitude, which one stray event would drag below the data.
    mags = gutenberg_richter_sample(b=1.0, mc=3.0, n=20000, seed=0)
    clean = b_stability(mags)
    with_outlier = b_stability(np.append(mags, 0.2))
    assert with_outlier.thresholds.tolist() == clean.thresholds.tolist()
    assert mc_b_stability(np.append(mags, 0.2)) == mc_b_stability(mags)


def test_the_exact_binned_estimator_reproduces_the_independent_package():
    # seismostats (ETH Zurich) computes the Tinti and Mulargia estimator and
    # returns 0.846804 for this catalogue at Mc 3.5. Reproducing an independent
    # implementation to six decimals is the strongest correctness evidence here,
    # so the number is pinned rather than left to a script nobody runs.
    mags = pd.read_csv(DATA / "kahramanmaras_180d.csv")["mw"].to_numpy()
    assert b_value_tinti(mags, 3.5).b == pytest.approx(0.846804, abs=1e-6)
    assert b_value_tinti(mags, 3.5).n == 1529


def test_the_half_bin_form_approaches_the_exact_binned_estimator():
    # Aki with Utsu's half-bin offset is the first-order approximation of the
    # exact binned MLE, so the gap must shrink with the bin width. If it stopped
    # shrinking, one of the two would be wrong.
    mags = pd.read_csv(DATA / "kahramanmaras_180d.csv")["mw"].to_numpy()
    gaps = [
        abs(b_value_aki(mags, 3.5, dm=dm).b - b_value_tinti(mags, 3.5, dm=dm).b)
        for dm in (0.2, 0.1, 0.05, 0.01)
    ]
    assert gaps == sorted(gaps, reverse=True)
    assert gaps[-1] < 1e-4


def test_the_two_b_estimators_agree_on_the_reference_catalogue_to_half_a_percent():
    mags = pd.read_csv(DATA / "kahramanmaras_180d.csv")["mw"].to_numpy()
    aki = b_value_aki(mags, 3.5).b
    tinti = b_value_tinti(mags, 3.5).b
    assert abs(aki / tinti - 1) < 0.005
    assert b_value_tinti(mags, 3.5).n == b_value_aki(mags, 3.5).n


def test_the_standard_error_uses_the_finite_sample_denominator():
    """Shi and Bolt divide by n(n-1), not by n squared.

    The difference is a factor sqrt(n/(n-1)): invisible at the 1,529 events of
    the reference catalogue, where it moves sigma by six millionths, and 1 per
    cent at fifty events. Every other test that touches sigma either compares a
    ratio in which the denominator cancels, or allows a tolerance wider than the
    difference, so replacing n(n-1) with n squared left the suite green.
    """
    rng = np.random.default_rng(0)
    mags = np.round(2.0 + rng.exponential(1 / (1.0 * np.log(10)), 50), 1)
    estimate = b_value_aki(mags, 2.0)
    n = estimate.n
    mean = mags[mags >= 2.0 - 0.05].mean()
    above = mags[mags >= 2.0 - 0.05]
    expected = (
        2.30 * estimate.b**2 * np.sqrt(((above - mean) ** 2).sum() / (n * (n - 1)))
    )
    assert estimate.sigma == pytest.approx(expected, rel=1e-12)
    # And it is not the n-squared form, which the suite could not previously see.
    wrong = 2.30 * estimate.b**2 * np.sqrt(((above - mean) ** 2).sum() / (n * n))
    assert estimate.sigma != pytest.approx(wrong, rel=1e-6)


def test_a_threshold_is_judged_usable_by_the_events_the_estimator_will_use():
    """The count that gates a threshold must match the sample taken at it.

    The gate counts events at or above threshold - dm/2, because that is the
    half-bin edge b_value_aki itself selects from. Counting from the threshold
    instead makes the gate reject a threshold the estimator could have used, and
    on the reference catalogue at dm 0.1 the two counts differ at every step.
    """
    catalog = pd.read_csv(DATA / "kahramanmaras_180d.csv")
    mags = catalog["mw"].to_numpy()
    curve = b_stability(mags, dm=0.1, min_events=50)
    for threshold, n in zip(curve.thresholds, curve.n, strict=True):
        # Every threshold that survived the gate reports the count the estimator
        # actually used, which is the half-bin one.
        assert n == int((mags >= threshold - 0.1 / 2 - 1e-9).sum())
