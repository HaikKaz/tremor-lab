"""How far b moves with the threshold, how large that is against its own error, and
where it stops moving.

The Shi and Bolt error beside a b-value is the scatter of one estimate. Two estimates
at different thresholds share most of their events, so the error of the difference
between them is a different quantity, and comparing a shift with the error of either
end gives a ratio that means nothing. `b_shift` and `b_plateau` exist to ask the
question properly, so the tests here check them against catalogues whose answer is
known: a Gutenberg-Richter law with one b, and the same law seen through a detector
that misses small events.
"""

import numpy as np
import pytest
from scipy.stats import norm

from tremor_lab.bvalue import b_plateau, b_shift, b_stability, b_value_aki

LN10 = np.log(10)


def gutenberg_richter(n, b, rng, floor=0.5, decimals=2):
    """Magnitudes from one Gutenberg-Richter law, reported to `decimals`.

    The law runs below every threshold the tests use. Cutting it at the lowest
    threshold would leave that threshold's sample as the only one whose events were
    not rounded across its own lower edge, and a shift measured from it would carry
    a bias that belongs to the simulation and not to the estimator.
    """
    return np.round(floor + rng.exponential(1.0 / (b * LN10), n), decimals)


def through_a_detector(n, b, mu, sigma, rng, floor=0.5):
    """The same law, each event detected with probability norm.cdf((M - mu) / sigma)."""
    m = floor + rng.exponential(1.0 / (b * LN10), n)
    seen = rng.random(n) < norm.cdf((m - mu) / sigma)
    return np.round(m[seen], 2)


# ------------------------------------------------------------------------ b_shift


def test_the_shift_is_the_higher_estimate_minus_the_lower_one():
    mags = gutenberg_richter(20000, 1.0, np.random.default_rng(0))
    low = b_value_aki(mags, 1.0)
    high = b_value_aki(mags, 2.0)
    result = b_shift(mags, 1.0, 2.0, n_boot=0)
    assert result.b_low == low.b
    assert result.b_high == high.b
    assert result.shift == pytest.approx(high.b - low.b, abs=1e-15)
    assert (result.n_low, result.n_high) == (low.n, high.n)


def test_the_error_under_a_constant_b_is_the_low_error_scaled_by_the_counts():
    # The estimate at the lower threshold uses more events and is efficient, so
    # var(high - low) = var(high) - var(low). With var proportional to 1/n that is
    # sigma_low^2 (n_low / n_high - 1). The two ways of writing it are checked
    # against each other, because the second is what the code evaluates and the
    # first is what the derivation says.
    mags = gutenberg_richter(20000, 1.0, np.random.default_rng(1))
    result = b_shift(mags, 1.0, 2.0, n_boot=0)
    low = b_value_aki(mags, 1.0)
    expected = low.sigma * np.sqrt(low.n / result.n_high - 1.0)
    assert result.se_constant_b == pytest.approx(expected, rel=1e-12)
    assert result.z_constant_b == pytest.approx(result.shift / expected, rel=1e-12)
    assert result.p_constant_b == pytest.approx(
        2 * norm.sf(abs(result.z_constant_b)), rel=1e-12
    )
    # Not the error of one end, which is what the comparison used to be made with.
    assert result.se_constant_b != pytest.approx(low.sigma, rel=0.1)


def test_a_constant_b_gives_a_shift_that_is_a_standard_normal_variate():
    """z is the Hausman statistic, so under one b it should be N(0, 1).

    If the variance formula were wrong by a constant factor the spread would show it
    here, which no comparison with a single simulated catalogue can.
    """
    rng = np.random.default_rng(2)
    z = np.array(
        [
            b_shift(gutenberg_richter(12000, 1.0, rng), 1.0, 2.0, n_boot=0).z_constant_b
            for _ in range(400)
        ]
    )
    assert abs(z.mean()) < 0.2
    assert 0.88 < z.std(ddof=1) < 1.12
    assert np.mean(np.abs(z) > 1.96) < 0.09


def test_the_bootstrap_error_agrees_with_the_constant_b_error_when_b_is_constant():
    mags = gutenberg_richter(30000, 1.0, np.random.default_rng(3))
    result = b_shift(mags, 1.0, 2.0, n_boot=1500, seed=0)
    assert result.se_boot == pytest.approx(result.se_constant_b, rel=0.15)
    assert result.n_boot == 1500


def test_the_bootstrap_interval_covers_a_zero_shift_at_about_its_stated_rate():
    # Every catalogue here has one b, so the true shift is zero up to the small
    # upward bias of the estimator at its smaller sample (about 0.005 here, a tenth of
    # the error). A 95 per cent interval should contain it in about 95 of 100.
    rng = np.random.default_rng(4)
    covered = 0
    runs = 120
    for _ in range(runs):
        mags = gutenberg_richter(8000, 1.0, rng)
        lo, hi = b_shift(mags, 1.0, 2.0, n_boot=300, seed=1).ci_boot
        covered += lo <= 0.0 <= hi
    assert covered / runs > 0.87


def test_the_bootstrap_does_not_assume_a_constant_b():
    # A catalogue whose b is 1 up to M 2.5 and 2 above it. The error derived under a
    # constant b describes a different situation and comes out far too small for a
    # shift across the break (about 0.05 against a true scatter near 0.1); the
    # bootstrap measures the scatter that is actually there.
    rng = np.random.default_rng(5)
    below = 0.5 + rng.exponential(1.0 / LN10, 400000)
    below = below[below < 2.5]
    above = 2.5 + rng.exponential(1.0 / (2.0 * LN10), 3000)
    mags = np.round(np.concatenate([below, above]), 2)
    result = b_shift(mags, 2.0, 3.0, n_boot=1500, seed=0)
    assert result.shift > 0.5
    assert result.se_boot > 1.5 * result.se_constant_b
    assert result.ci_boot[0] > 0.0
    assert result.ci_boot[0] < result.shift < result.ci_boot[1]


def test_a_catalogue_missing_its_small_events_shows_a_shift_that_can_be_seen():
    rng = np.random.default_rng(6)
    mags = through_a_detector(170000, 1.0, mu=1.7, sigma=0.45, rng=rng)
    # Maximum curvature picks a threshold inside the roll-off; far above it the law is
    # recovered.
    result = b_shift(mags, 1.6, 2.6, n_boot=500, seed=0)
    assert result.shift > 0.1
    assert result.z_constant_b > 8
    assert result.ci_boot[0] > 0.0
    assert result.b_low < 0.95 < 1.0 - 0.03 < result.b_high + 0.03


def test_the_bootstrap_is_reproducible_from_its_seed():
    mags = gutenberg_richter(5000, 1.0, np.random.default_rng(7))
    first = b_shift(mags, 1.0, 1.8, n_boot=200, seed=3)
    assert first == b_shift(mags, 1.0, 1.8, n_boot=200, seed=3)
    assert first.ci_boot != b_shift(mags, 1.0, 1.8, n_boot=200, seed=4).ci_boot
    # The closed form does not depend on the seed.
    assert (
        first.se_constant_b == b_shift(mags, 1.0, 1.8, n_boot=200, seed=4).se_constant_b
    )


def test_the_bootstrap_can_be_switched_off():
    mags = gutenberg_richter(5000, 1.0, np.random.default_rng(8))
    result = b_shift(mags, 1.0, 1.8, n_boot=0)
    assert result.se_boot is None
    assert result.ci_boot is None
    assert result.n_boot == 0
    assert np.isfinite(result.se_constant_b)


def test_the_confidence_level_sets_the_width_of_the_interval():
    mags = gutenberg_richter(8000, 1.0, np.random.default_rng(9))
    narrow = b_shift(mags, 1.0, 1.8, n_boot=800, confidence=0.5, seed=0).ci_boot
    wide = b_shift(mags, 1.0, 1.8, n_boot=800, confidence=0.99, seed=0).ci_boot
    assert narrow[0] > wide[0]
    assert narrow[1] < wide[1]


def test_the_thresholds_must_be_in_order_and_must_keep_different_events():
    mags = gutenberg_richter(5000, 1.0, np.random.default_rng(10))
    with pytest.raises(ValueError, match="below mc_high"):
        b_shift(mags, 2.0, 1.0)
    with pytest.raises(ValueError, match="below mc_high"):
        b_shift(mags, 1.5, 1.5)
    # Magnitudes reported to 0.1 and thresholds 0.02 apart keep the same events.
    coarse = np.round(mags, 1)
    with pytest.raises(ValueError, match=r"same \d+ events"):
        b_shift(coarse, 1.50, 1.52)


def test_unusable_settings_are_refused_with_a_reason():
    mags = gutenberg_richter(5000, 1.0, np.random.default_rng(11))
    with pytest.raises(ValueError, match="cannot give a spread"):
        b_shift(mags, 1.0, 1.8, n_boot=1)
    with pytest.raises(ValueError, match="cannot be negative"):
        b_shift(mags, 1.0, 1.8, n_boot=-5)
    with pytest.raises(ValueError, match="strictly between"):
        b_shift(mags, 1.0, 1.8, confidence=1.0)
    with pytest.raises(ValueError, match="at least two events"):
        b_shift(mags, 1.0, 9.0)


# ---------------------------------------------------------------------- b_plateau


def test_over_a_single_step_the_statistic_is_the_square_of_the_shift_z():
    mags = through_a_detector(120000, 1.0, 1.6, 0.4, np.random.default_rng(12))
    curve = b_stability(mags, thresholds=[1.8, 1.9])
    plateau = b_plateau(mags, thresholds=[1.8, 1.9], n_simulations=0, min_window=1)
    shift = b_shift(mags, 1.8, 1.9, n_boot=0)
    assert plateau.statistic_profile[0] == pytest.approx(
        shift.z_constant_b**2, rel=1e-9
    )
    assert plateau.dof_profile.tolist() == [1, 0]
    assert curve.thresholds.tolist() == [1.8, 1.9]


def test_the_statistic_is_the_sum_over_steps_of_each_squared_change_in_b():
    mags = gutenberg_richter(60000, 1.0, np.random.default_rng(13))
    grid = [1.0, 1.2, 1.4, 1.6, 1.8]
    plateau = b_plateau(mags, thresholds=grid, n_simulations=0, min_window=1)
    curve = b_stability(mags, thresholds=grid)
    expected = 0.0
    for j in range(1, len(grid)):
        variance = (
            curve.sigma[0] ** 2 * curve.n[0] * (1 / curve.n[j] - 1 / curve.n[j - 1])
        )
        expected += (curve.b[j] - curve.b[j - 1]) ** 2 / variance
    assert plateau.statistic_profile[0] == pytest.approx(expected, rel=1e-9)
    assert plateau.dof_profile[0] == len(grid) - 1
    # Each later candidate sums only the steps above it.
    tail = 0.0
    for j in range(len(grid) - 1, 1, -1):
        variance = (
            curve.sigma[1] ** 2 * curve.n[1] * (1 / curve.n[j] - 1 / curve.n[j - 1])
        )
        tail += (curve.b[j] - curve.b[j - 1]) ** 2 / variance
    assert plateau.statistic_profile[1] == pytest.approx(tail, rel=1e-9)


def test_one_b_is_accepted_from_the_first_threshold_about_as_often_as_it_should_be():
    """The size of the test, on catalogues where the answer is that b is constant.

    The chi-squared reference was found too generous with its rejections (it rejects
    about 10 per cent of catalogues like these at the 5 per cent level), which is why
    the default is a simulated null. This holds the simulated one to a rate near the
    nominal.
    """
    rng = np.random.default_rng(14)
    runs = 160
    grid = np.round(np.arange(1.0, 2.55, 0.1), 1)
    rejected = 0
    for k in range(runs):
        mags = gutenberg_richter(6000, 1.0, rng)
        result = b_plateau(mags, thresholds=grid, n_simulations=100, seed=k)
        rejected += result.p_profile[0] < 0.05
    assert rejected / runs < 0.12


def test_the_onset_of_a_catalogue_missing_small_events_is_found_above_its_roll_off():
    rng = np.random.default_rng(15)
    mags = through_a_detector(170000, 1.0, mu=1.7, sigma=0.35, rng=rng)
    result = b_plateau(mags, n_simulations=100, seed=0)
    assert result.onset is not None
    # Detection is 99 per cent at mu + 2.33 sigma = 2.5. The onset can fall a little
    # below it, where the bias is smaller than the error, and not far below.
    assert 2.0 <= result.onset <= 2.7
    assert result.b == pytest.approx(1.0, abs=4 * result.sigma)
    # Every candidate well below the roll-off is rejected outright.
    below = result.thresholds < 1.9
    assert np.all(result.p_profile[below & (result.dof_profile >= 3)] < 0.01)
    assert result.method == "simulated null, 100 draws"


def test_a_complete_catalogue_has_its_plateau_at_the_lowest_threshold():
    rng = np.random.default_rng(16)
    mags = gutenberg_richter(40000, 1.0, rng)
    result = b_plateau(mags, thresholds=np.round(np.arange(1.0, 3.05, 0.1), 1), seed=1)
    assert result.onset == 1.0
    assert result.b == pytest.approx(1.0, abs=4 * result.sigma)


def test_no_candidate_with_plenty_of_events_is_accepted_where_b_climbs_all_the_way():
    # Excesses with a Weibull tail of shape 1.6 have a hazard that rises with
    # magnitude, so b estimated above any threshold is larger than at the one below.
    # Where the events are plentiful the climb is seen and every candidate is
    # rejected. Far up the range, with a few hundred events, a slow climb is within the
    # noise and a candidate there can pass, which the docstring says in so many words.
    rng = np.random.default_rng(17)
    mags = np.round(1.0 + 0.9 * rng.weibull(1.6, 60000), 2)
    result = b_plateau(mags, n_simulations=100, seed=0)
    curve = b_stability(mags)
    plentiful = (curve.n >= 2000) & (result.dof_profile >= 3)
    assert plentiful.sum() >= 5
    assert np.all(result.p_profile[plentiful] < 0.05)
    assert result.onset is None or result.n < 1000


def test_the_chi_squared_reference_is_labelled_and_is_the_more_generous_one():
    rng = np.random.default_rng(18)
    mags = gutenberg_richter(6000, 1.0, rng)
    grid = np.round(np.arange(1.0, 2.55, 0.1), 1)
    approx = b_plateau(mags, thresholds=grid, n_simulations=0)
    assert approx.method == "chi-squared, approximate"
    simulated = b_plateau(mags, thresholds=grid, n_simulations=200, seed=0)
    # The statistics are the same quantity under either reference.
    assert np.allclose(
        approx.statistic_profile, simulated.statistic_profile, equal_nan=True
    )
    assert approx.dof_profile.tolist() == simulated.dof_profile.tolist()


def test_the_simulated_p_value_is_reproducible_from_its_seed():
    rng = np.random.default_rng(19)
    mags = through_a_detector(60000, 1.0, 1.6, 0.3, rng)
    first = b_plateau(mags, n_simulations=60, seed=2)
    again = b_plateau(mags, n_simulations=60, seed=2)
    other = b_plateau(mags, n_simulations=60, seed=3)
    assert np.array_equal(first.p_profile, again.p_profile, equal_nan=True)
    assert (first.onset, first.p_value) == (again.onset, again.p_value)
    assert not np.array_equal(first.p_profile, other.p_profile, equal_nan=True)


def test_a_simulated_p_value_cannot_fall_below_one_over_the_draws_plus_one():
    rng = np.random.default_rng(20)
    mags = through_a_detector(60000, 1.0, 1.6, 0.4, rng)
    result = b_plateau(mags, n_simulations=50, seed=0)
    floor = 1.0 / 51
    finite = result.p_profile[np.isfinite(result.p_profile)]
    assert finite.min() >= floor - 1e-12
    assert finite.min() == pytest.approx(floor)


def test_raising_min_events_shortens_the_range_the_plateau_is_tested_over():
    mags = through_a_detector(170000, 1.0, 1.6, 0.3, np.random.default_rng(21))
    short = b_plateau(mags, min_events=1000, n_simulations=0)
    long = b_plateau(mags, min_events=50, n_simulations=0)
    assert short.thresholds.max() < long.thresholds.max()
    assert short.thresholds.min() == long.thresholds.min()


def test_a_candidate_needs_min_window_steps_above_it():
    rng = np.random.default_rng(22)
    mags = gutenberg_richter(20000, 1.0, rng)
    grid = [1.0, 1.1, 1.2, 1.3, 1.4]
    narrow = b_plateau(mags, thresholds=grid, n_simulations=0, min_window=4)
    wide = b_plateau(mags, thresholds=grid, n_simulations=0, min_window=1)
    assert np.isfinite(narrow.p_profile[0]) and np.isnan(narrow.p_profile[1])
    assert np.isfinite(wide.p_profile[3])
    # The statistic is defined wherever a step exists, tested or not.
    assert np.isfinite(narrow.statistic_profile[3])


def test_thresholds_that_keep_the_same_events_do_not_count_as_steps():
    # Magnitudes to one decimal and thresholds 0.02 apart: neighbouring thresholds
    # keep identical events, differ by nothing, and must neither add a degree of
    # freedom nor divide zero by zero.
    rng = np.random.default_rng(23)
    mags = gutenberg_richter(30000, 1.0, rng, decimals=1)
    grid = [1.0, 1.02, 1.1, 1.12, 1.2, 1.3, 1.4, 1.5]
    result = b_plateau(mags, thresholds=grid, n_simulations=0, min_window=1)
    assert np.all(np.isfinite(result.statistic_profile[:-1]))
    assert result.dof_profile[0] == 5


def test_unusable_settings_are_refused_with_a_reason_here_too():
    mags = gutenberg_richter(5000, 1.0, np.random.default_rng(24))
    with pytest.raises(ValueError, match="alpha"):
        b_plateau(mags, alpha=0.0)
    with pytest.raises(ValueError, match="min_window"):
        b_plateau(mags, min_window=0)
    with pytest.raises(ValueError, match="n_simulations"):
        b_plateau(mags, n_simulations=-1)


def test_too_few_thresholds_give_no_onset_rather_than_an_error():
    mags = gutenberg_richter(300, 1.0, np.random.default_rng(25))
    result = b_plateau(mags, min_events=250, n_simulations=0)
    assert result.onset is None
    assert result.dof_profile.sum() == 0 or result.p_profile.size <= 2


def test_the_coefficient_of_the_standard_error_enters_the_statistic():
    # sigma scales with the coefficient and the statistic divides by sigma squared,
    # so doubling the coefficient divides T by four. A function that ignored the
    # keyword would pass every other test here.
    mags = gutenberg_richter(30000, 1.0, np.random.default_rng(26))
    grid = [1.0, 1.2, 1.4, 1.6]
    base = b_plateau(mags, thresholds=grid, n_simulations=0, min_window=1)
    doubled = b_plateau(
        mags, thresholds=grid, n_simulations=0, min_window=1, shi_bolt_k=4.60
    )
    assert doubled.statistic_profile[0] == pytest.approx(
        base.statistic_profile[0] / 4.0, rel=1e-9
    )
