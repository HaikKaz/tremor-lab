"""The Omori-Utsu fit: the likelihood itself, the k shortcut, and parameter recovery."""

import numpy as np
import pytest
from scipy.optimize import minimize

from tremor_lab import constants
from tremor_lab.omori import bootstrap_omori, fit_omori, omori_nll


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
    with pytest.raises(ValueError, match="at least two positive times"):
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
