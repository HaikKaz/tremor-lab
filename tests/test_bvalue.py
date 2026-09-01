"""The Aki b-value: exact arithmetic, then recovery of a known b from synthetic
draws."""

import numpy as np
import pytest

from tremor_lab import constants
from tremor_lab.bvalue import b_value_aki

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
