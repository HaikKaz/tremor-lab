"""The Omori-Utsu fit: the likelihood itself, the k shortcut, and parameter recovery."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize

from tremor_lab import constants
from tremor_lab.cli import fit_verdict
from tremor_lab.omori import (
    bootstrap_omori,
    fit_omori,
    omori_fit_test,
    omori_nll,
)

DATA = Path(__file__).parent / "data"


def omori_sample(n, p, c, t_end, seed):
    """Draw n occurrence times whose density on (0, t_end] is proportional to
    (c + t)^-p, by inverting the cumulative distribution."""
    rng = np.random.default_rng(seed)
    lo, hi = c ** (1 - p), (c + t_end) ** (1 - p)
    return (lo + rng.random(n) * (hi - lo)) ** (1 / (1 - p)) - c


def integrated_rate(c, p, t_end):
    return ((c + t_end) ** (1 - p) - c ** (1 - p)) / (1 - p)


def test_the_likelihood_matches_a_direct_evaluation():
    t = np.array([0.5, 1.0, 2.0])
    k, c, p = 10.0, 0.4, 1.2
    expected = k * integrated_rate(c, p, 3.0) - np.sum(np.log(k) - p * np.log(t + c))
    assert omori_nll((k, c, p), t, 3.0) == pytest.approx(expected)


def test_the_likelihood_rejects_non_positive_parameters():
    t = np.array([0.5, 1.0])
    assert omori_nll((-1.0, 0.4, 1.2), t, 2.0) == np.inf
    assert omori_nll((10.0, 0.0, 1.2), t, 2.0) == np.inf
    assert omori_nll((10.0, 0.4, -0.1), t, 2.0) == np.inf


def test_the_integrated_rate_is_continuous_through_p_equals_one():
    # p = 1 needs a logarithm rather than the power form; the two must agree in
    # the limit or the likelihood surface would have a seam at p = 1.
    t = np.array([0.5, 1.0, 2.0, 10.0])
    at_one = omori_nll((10.0, 0.4, 1.0), t, 100.0)
    just_above = omori_nll((10.0, 0.4, 1.0 + 1e-9), t, 100.0)
    assert at_one == pytest.approx(just_above, abs=1e-6)


@pytest.mark.filterwarnings("ignore:invalid value encountered in subtract")
def test_the_profiled_fit_finds_the_same_optimum_as_a_full_three_parameter_search():
    t = omori_sample(3000, p=1.15, c=0.5, t_end=180.0, seed=0)
    profiled = fit_omori(t, t_end=180.0)
    full = minimize(
        omori_nll,
        (100.0, 0.5, 1.1),
        args=(t, 180.0),
        method="Nelder-Mead",
        options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 20000},
    )
    k, c, p = full.x
    assert profiled.p == pytest.approx(p, rel=1e-3)
    assert profiled.c == pytest.approx(c, rel=1e-3)
    assert profiled.k == pytest.approx(k, rel=1e-3)
    # The decisive comparison: the profiled optimum is at least as likely.
    assert omori_nll((profiled.k, profiled.c, profiled.p), t, 180.0) <= full.fun + 1e-6


def test_the_fit_does_not_depend_on_the_starting_point():
    t = omori_sample(3000, p=1.15, c=0.5, t_end=180.0, seed=1)
    reference = fit_omori(t, t_end=180.0)
    for c0, p0 in [(0.05, 0.7), (2.0, 1.6), (0.5, 1.0)]:
        alternative = fit_omori(t, t_end=180.0, c0=c0, p0=p0)
        assert alternative.p == pytest.approx(reference.p, rel=1e-4)
        assert alternative.c == pytest.approx(reference.c, rel=1e-3)


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize(("p_true", "c_true"), [(1.15, 0.5), (1.05, 0.1), (0.9, 0.2)])
def test_known_decay_parameters_are_recovered(p_true, c_true, seed):
    t = omori_sample(20000, p_true, c_true, t_end=180.0, seed=seed)
    fit = fit_omori(t, t_end=180.0)
    assert fit.p == pytest.approx(p_true, abs=0.02)
    assert fit.c == pytest.approx(c_true, rel=0.12)


def test_the_fitted_rate_integrates_to_the_observed_event_count():
    # k = n / integral is exact, not approximate: the fitted rate integrated over the
    # observation interval returns the number of events that were actually seen.
    t = omori_sample(3000, 1.15, 0.5, t_end=180.0, seed=0)
    fit = fit_omori(t, t_end=180.0)
    assert fit.k * integrated_rate(fit.c, fit.p, 180.0) == pytest.approx(fit.n)


def test_the_productivity_constant_is_recovered():
    # k is the least well-constrained of the three. It is not searched for but derived
    # from the fitted c and p, and the integral is steep in both, so its sampling error
    # stays near a few per cent where p is already recovered to a fraction of one.
    p_true, c_true, k_true = 1.15, 0.5, 2000.0
    n = round(k_true * integrated_rate(c_true, p_true, 180.0))
    t = omori_sample(n, p_true, c_true, t_end=180.0, seed=0)
    assert fit_omori(t, t_end=180.0).k == pytest.approx(k_true, rel=0.10)


def test_the_observation_interval_defaults_to_the_last_event():
    t = omori_sample(2000, 1.15, 0.5, t_end=180.0, seed=0)
    assert fit_omori(t) == fit_omori(t, t_end=float(t.max()))


def test_a_longer_observation_interval_changes_the_answer():
    # t_end is not cosmetic: claiming to have watched for 180 days when the last
    # event was at 179.6 is a different likelihood.
    t = omori_sample(2000, 1.15, 0.5, t_end=180.0, seed=0)
    assert fit_omori(t, t_end=180.0).p != fit_omori(t, t_end=3650.0).p


def test_non_positive_times_are_discarded():
    t = omori_sample(2000, 1.15, 0.5, t_end=180.0, seed=0)
    assert fit_omori(np.concatenate([[-5.0, 0.0], t])).n == 2000


def test_too_few_times_is_an_error():
    with pytest.raises(ValueError, match="at least two times after"):
        fit_omori([1.0])


def test_the_bootstrap_is_reproducible_from_its_seed():
    t = omori_sample(1500, 1.15, 0.5, t_end=180.0, seed=0)
    first = bootstrap_omori(t, n_boot=30, seed=7)
    assert first == bootstrap_omori(t, n_boot=30, seed=7)
    assert first != bootstrap_omori(t, n_boot=30, seed=8)


def test_the_bootstrap_spread_covers_the_fitted_value():
    t = omori_sample(3000, 1.15, 0.5, t_end=180.0, seed=0)
    fit = fit_omori(t, t_end=180.0)
    spread = bootstrap_omori(t, n_boot=40, t_end=180.0, seed=0)
    assert abs(fit.p - 1.15) < 4 * spread.p_std
    assert spread.p_std > 0 and spread.c_std > 0


def test_the_bootstrap_narrows_as_the_catalogue_grows():
    small = bootstrap_omori(
        omori_sample(500, 1.15, 0.5, 180.0, 0), n_boot=40, t_end=180.0, seed=0
    )
    large = bootstrap_omori(
        omori_sample(20000, 1.15, 0.5, 180.0, 0), n_boot=40, t_end=180.0, seed=0
    )
    assert large.p_std < small.p_std


def test_bootstrap_size_follows_the_constant(monkeypatch):
    monkeypatch.setattr(constants, "N_BOOT", 12)
    t = omori_sample(500, 1.15, 0.5, t_end=180.0, seed=0)
    assert bootstrap_omori(t).n_boot == 12


def test_a_negative_bootstrap_count_is_refused():
    t = omori_sample(500, 1.15, 0.5, t_end=180.0, seed=0)
    with pytest.raises(ValueError, match="cannot be negative"):
        bootstrap_omori(t, n_boot=-5)


def test_a_fit_can_start_after_the_origin():
    # The early hours of a real sequence are the least complete; being able to
    # exclude them is how c is diagnosed rather than merely reported.
    t = omori_sample(4000, 1.15, 0.5, t_end=180.0, seed=0)
    late = fit_omori(t, t_end=180.0, t_start=1.0)
    assert late.n == int((t > 1.0).sum())
    assert late.p == pytest.approx(1.15, abs=0.08)


def test_starting_at_zero_is_unchanged():
    t = omori_sample(2000, 1.15, 0.5, t_end=180.0, seed=1)
    assert fit_omori(t, t_end=180.0) == fit_omori(t, t_end=180.0, t_start=0.0)


def test_times_from_the_fitted_model_pass_the_fit_test():
    t = omori_sample(3000, 1.15, 0.5, t_end=180.0, seed=0)
    fit = fit_omori(t, t_end=180.0)
    assert omori_fit_test(t, fit, t_end=180.0, n_simulations=60).p_value > 0.05


def test_times_that_are_not_an_omori_decay_are_rejected():
    # The negative control: without it a test that never fails is evidence of
    # nothing. A uniform sequence will not do, because the model covers it at
    # p = 0, where the rate is constant. A rate that *rises* with time cannot be
    # k / (c + t)^p for any positive p, so it must be rejected.
    rng = np.random.default_rng(0)
    t = np.sort(180.0 * rng.random(3000) ** (1 / 3))
    fit = fit_omori(t, t_end=180.0)
    assert omori_fit_test(t, fit, t_end=180.0, n_simulations=60).p_value < 0.05


def test_a_second_sequence_inside_the_window_is_detected():
    # One Omori decay plus a burst starting on day 30, which is what a large
    # aftershock does to a catalogue. The model should be rejected.
    main = omori_sample(2500, 1.15, 0.5, t_end=180.0, seed=0)
    burst = 30.0 + omori_sample(1200, 1.2, 0.3, t_end=150.0, seed=1)
    t = np.sort(np.concatenate([main, burst]))
    fit = fit_omori(t, t_end=180.0)
    assert omori_fit_test(t, fit, t_end=180.0, n_simulations=60).p_value < 0.05


def test_the_asymptotic_p_value_is_anti_conservative_and_is_labelled_as_such():
    # The parameters are fitted to the times being tested, so the fitted curve
    # hugs the data and the textbook Kolmogorov p-value is far too generous. On
    # sequences drawn from the model it essentially never rejects. The default
    # calibrates by parametric bootstrap; the asymptotic form remains available
    # but must say what it is.
    t = omori_sample(1200, 1.15, 0.5, t_end=180.0, seed=0)
    fit = fit_omori(t, t_end=180.0)
    asymptotic = omori_fit_test(t, fit, t_end=180.0, n_simulations=0)
    calibrated = omori_fit_test(t, fit, t_end=180.0, n_simulations=60)
    assert asymptotic.method == "asymptotic, uncalibrated"
    assert calibrated.method.startswith("parametric bootstrap")
    assert asymptotic.statistic == pytest.approx(calibrated.statistic)
    assert asymptotic.p_value > calibrated.p_value


def test_the_calibrated_p_value_is_reproducible_from_its_seed():
    t = omori_sample(600, 1.15, 0.5, t_end=180.0, seed=0)
    fit = fit_omori(t, t_end=180.0)
    first = omori_fit_test(t, fit, t_end=180.0, n_simulations=40, seed=5)
    assert first == omori_fit_test(t, fit, t_end=180.0, n_simulations=40, seed=5)
    assert first != omori_fit_test(t, fit, t_end=180.0, n_simulations=40, seed=6)


def test_calibration_rejects_the_reference_sequence_the_asymptotic_test_accepts():
    # The Kahramanmaras window contains the M 7.6 Elbistan event and its own
    # aftershocks, so a single Omori decay should not describe it. Only the
    # calibrated test notices.
    frame = pd.read_csv(Path(__file__).parent / "data" / "kahramanmaras_180d.csv")
    times = frame.loc[frame["mw"] >= 3.5, "dt_days"].to_numpy()
    fit = fit_omori(times)
    assert omori_fit_test(times, fit, n_simulations=0).p_value > 0.4
    assert omori_fit_test(times, fit, n_simulations=120, seed=0).p_value < 0.10


@pytest.mark.full_calibration
def test_the_published_default_is_paid_in_full():
    """One call at the real replicate count, so the shipped default is exercised.

    The rest of the suite lowers it for speed (see conftest). If nothing ran at
    600, a default that had drifted to something unusable would pass every test.
    """
    catalog = pd.read_csv(DATA / "kahramanmaras_180d.csv")
    times = catalog.loc[catalog["mw"] >= 3.5, "dt_days"].to_numpy()
    fit = fit_omori(times)
    test = omori_fit_test(times, fit, seed=0)

    assert constants.N_FIT_SIMULATIONS == 600
    assert test.method == "parametric bootstrap, 600 replicates"
    # The reference sequence sits on the threshold: 0.042 with a standard error
    # of 0.008. Both implementations reach it, and neither can call it either
    # way, which is why the verdict says borderline rather than choosing.
    assert test.p_value == pytest.approx(0.042, abs=0.005)
    assert test.p_value_se == pytest.approx(0.008, abs=0.002)
    assert "borderline" in fit_verdict(test)


def test_the_observation_interval_is_held_fixed_across_every_refit():
    """The interval belongs to the catalogue, not to the resample.

    Each bootstrap resample has its own last event, and letting each refit take
    the interval from its own draw makes the answer identical whatever interval
    the caller asked for - which is how this was found: deleting `t_end=t_end`
    from the refit left the whole suite green while understating the uncertainty
    on c by 62 per cent at a 365-day interval.
    """
    catalog = pd.read_csv(DATA / "kahramanmaras_180d.csv")
    times = catalog.loc[catalog["mw"] >= 3.5, "dt_days"].to_numpy()
    spreads = [
        bootstrap_omori(times, n_boot=40, seed=0, t_end=t_end).c_std
        for t_end in (180.0, 365.0, 3650.0)
    ]
    # A longer observation interval over the same events means the decay is less
    # constrained, so the spread has to grow with it.
    assert spreads[0] < spreads[1] < spreads[2]


def test_a_single_resample_is_refused_rather_than_reported_as_nan():
    """A spread over one number has no denominator.

    ddof=1 on a single value is a division by zero, and the report printed
    "+/- nan" beside p and c as though it were a published uncertainty. Zero
    means "do not estimate the spread" and is allowed; one cannot be.
    """
    times = np.linspace(0.01, 180.0, 500)
    with pytest.raises(ValueError, match="cannot give a spread"):
        bootstrap_omori(times, n_boot=1)
    # Zero is a legitimate way to skip it, and two is the smallest real estimate.
    assert bootstrap_omori(times, n_boot=0).n_boot == 0
    assert np.isfinite(bootstrap_omori(times, n_boot=2, seed=0).p_std)
