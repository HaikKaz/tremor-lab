"""Modified Omori-Utsu aftershock decay, fitted by maximum likelihood."""

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import minimize
from scipy.stats import kstest

from tremor_lab import constants


class Omori(NamedTuple):
    """Modified Omori-Utsu parameters and the number of events fitted."""

    p: float
    c: float
    k: float
    n: int


class OmoriBootstrap(NamedTuple):
    """Bootstrap standard errors of the Omori-Utsu decay and offset parameters."""

    p_std: float
    c_std: float
    n_boot: int


def omori_nll(params: tuple[float, float, float], times_days: ArrayLike, t_end: float):
    """
    Negative log-likelihood of the modified Omori-Utsu process on (0, t_end].

    The rate is n(t) = k / (c + t)^p, so the log-likelihood of unbinned occurrence
    times is sum(log k - p log(t + c)) minus the integrated rate over the observation
    interval.

    Parameters
    ----------
    params : tuple of float
        (k, c, p), ordered as in the reference implementation this package replaces.
    times_days : array_like
        Elapsed times since the mainshock in days, strictly positive.
    t_end : float
        End of the observation interval in days.

    Returns
    -------
    float
        Negative log-likelihood; infinite where a parameter is non-positive.

    References
    ----------
    Ogata, Y. (1983). Utsu, T., Ogata, Y. and Matsu'ura, R. S. (1995).
    """
    k, c, p = params
    if k <= 0 or c <= 0 or p <= 0:
        return np.inf
    t = np.asarray(times_days, float)
    return float(
        k * _integrated_rate(c, p, t_end) - np.sum(np.log(k) - p * np.log(t + c))
    )


def fit_omori(
    times_days: ArrayLike,
    t_end: float | None = None,
    c0: float | None = None,
    p0: float | None = None,
    t_start: float = 0.0,
) -> Omori:
    """
    Maximum-likelihood fit of n(t) = k / (c + t)^p to unbinned post-mainshock times.

    Only c and p are searched. At any (c, p) the likelihood is maximised by
    k = n / integral, so k is recovered in closed form rather than searched for. This
    removes one dimension and any dependence on a starting value for k, and gives the
    same optimum as the three-parameter search, which `test_omori.py` asserts directly.

    Parameters
    ----------
    times_days : array_like
        Elapsed times since the mainshock in days. Non-positive times are discarded,
        since the rate is undefined at the origin.
    t_end : float, optional
        End of the observation interval in days. Defaults to the largest elapsed time
        supplied, which is the convention the reference values were produced under.
        Note that this is the last event above threshold, not the nominal window
        length, and the two differ.
    c0, p0 : float, optional
        Starting point of the search. Default to `constants.OMORI_C0` (0.5) and
        `constants.OMORI_P0` (1.1).

    Returns
    -------
    Omori
        Named tuple of p, c, k and n, unrounded.

    Raises
    ------
    ValueError
        If fewer than two positive times are supplied, or the optimiser fails to
        converge.

    References
    ----------
    Ogata, Y. (1983). Utsu, T., Ogata, Y. and Matsu'ura, R. S. (1995).
    """
    t = np.sort(np.asarray(times_days, float))
    t = t[t > t_start]
    if t.size < 2:
        raise ValueError(
            f"Omori fit needs at least two times after {t_start}; got {t.size}"
        )
    if t_end is None:
        t_end = float(t[-1])
    c0 = constants.OMORI_C0 if c0 is None else c0
    p0 = constants.OMORI_P0 if p0 is None else p0

    search = minimize(
        _profiled_nll,
        (c0, p0),
        args=(t, t_end, t_start),
        method="Nelder-Mead",
        # fatol is absolute, so it must scale with the likelihood, which grows
        # with n; a fixed 1e-10 is unreachable on a large catalogue and the fit
        # is then rejected as unconverged despite sitting on the optimum.
        options={
            "xatol": 1e-8,
            "fatol": max(
                1e-10, 1e-10 * abs(_profiled_nll((c0, p0), t, t_end, t_start))
            ),
            "maxiter": 5000,
        },
    )
    if not search.success:
        raise ValueError(f"Omori fit did not converge: {search.message}")
    c, p = search.x
    return Omori(
        float(p), float(c), t.size / _integrated_rate(c, p, t_end, t_start), t.size
    )


def bootstrap_omori(
    times_days: ArrayLike,
    n_boot: int | None = None,
    t_end: float | None = None,
    seed: int = 0,
) -> OmoriBootstrap:
    """
    Standard errors of p and c from resampling the occurrence times with replacement.

    The observation interval is a property of the catalogue, not of the resample, so
    `t_end` is held at its original value across every refit rather than being taken
    from each resample's own last event.

    Parameters
    ----------
    times_days : array_like
        Elapsed times since the mainshock in days.
    n_boot : int, optional
        Number of resamples. Defaults to `constants.N_BOOT` (200).
    t_end : float, optional
        End of the observation interval. Defaults to the largest elapsed time supplied.
    seed : int, optional
        Seed of the resampling generator, so a reported uncertainty is reproducible.

    Returns
    -------
    OmoriBootstrap
        Named tuple of the standard deviations of p and c and the number of resamples.

    References
    ----------
    Efron, B. (1979).
    """
    n_boot = constants.N_BOOT if n_boot is None else n_boot
    if n_boot < 0:
        raise ValueError(f"n_boot cannot be negative; got {n_boot}")
    if n_boot == 1:
        # A standard deviation over one resample has no denominator, and the
        # report printed "+/- nan" beside p and c as though it were a published
        # uncertainty. Zero means "do not estimate the spread" and is allowed;
        # one means "estimate it from a single number", which cannot be done.
        raise ValueError(
            "n_boot of 1 cannot give a spread: a standard error needs at least "
            "two resamples. Use 0 to skip the bootstrap, or 2 or more to run it"
        )
    t = np.asarray(times_days, float)
    t = t[t > 0]
    if t_end is None:
        t_end = float(t.max())
    rng = np.random.default_rng(seed)
    fits = [
        fit_omori(rng.choice(t, size=t.size, replace=True), t_end=t_end)
        for _ in range(n_boot)
    ]
    # ddof=1, the divisor in Efron's own formula for a bootstrap standard error,
    # and the one the browser page already used. numpy's default of ddof=0 made
    # the two implementations report different numbers for the same quantity -
    # half a per cent apart at 200 resamples, and 41 per cent apart at the two
    # resamples the browser control will accept.
    return OmoriBootstrap(
        float(np.std([f.p for f in fits], ddof=1)),
        float(np.std([f.c for f in fits], ddof=1)),
        n_boot,
    )


def _integrated_rate(c: float, p: float, t_end: float, t_start: float = 0.0) -> float:
    """Integral of (c + t)^-p over (t_start, t_end], the expected count per unit k."""
    if abs(p - 1.0) < 1e-12:
        return float(np.log((t_end + c) / (t_start + c)))
    return float(((t_end + c) ** (1.0 - p) - (t_start + c) ** (1.0 - p)) / (1.0 - p))


# Returned where a parameter leaves the model's domain. A finite penalty rather
# than infinity, so the optimiser's convergence test never subtracts infinities.
_OUT_OF_DOMAIN = 1e12


def _profiled_nll(
    cp: tuple[float, float], t: np.ndarray, t_end: float, t_start: float = 0.0
) -> float:
    """Negative log-likelihood with k replaced by its maximising value n / integral."""
    c, p = cp
    if c <= 0 or p <= 0:
        return _OUT_OF_DOMAIN
    integral = _integrated_rate(c, p, t_end, t_start)
    if integral <= 0 or not np.isfinite(integral):
        return _OUT_OF_DOMAIN
    n = t.size
    # The terms can still overflow for extreme but in-domain parameters, and a
    # non-finite value here makes the optimiser's convergence test subtract
    # infinities. Anything unusable is reported as out of domain instead.
    value = p * np.sum(np.log(t + c)) - n * np.log(n / integral) + n
    return float(value) if np.isfinite(value) else _OUT_OF_DOMAIN


class FitTest(NamedTuple):
    """Kolmogorov-Smirnov test of a fitted decay against the times it was fitted to."""

    statistic: float
    p_value: float
    n: int
    method: str
    p_value_se: float


def omori_sample(
    n: int, p: float, c: float, t_end: float, t_start: float, rng
) -> np.ndarray:
    """Draw n occurrence times from k / (c + t)^p on (t_start, t_end].

    The density is inverted directly; k cancels, because conditioning on n
    removes the scale.
    """
    u = rng.random(n)
    if abs(p - 1.0) < 1e-12:
        lo = np.log(t_start + c)
        hi = np.log(t_end + c)
        return np.exp(lo + u * (hi - lo)) - c
    lo = (t_start + c) ** (1.0 - p)
    hi = (t_end + c) ** (1.0 - p)
    return (lo + u * (hi - lo)) ** (1.0 / (1.0 - p)) - c


def omori_fit_test(
    times_days: ArrayLike,
    fit: Omori,
    t_end: float | None = None,
    t_start: float = 0.0,
    n_simulations: int | None = None,
    seed: int = 0,
) -> FitTest:
    """
    Test whether the fitted decay actually describes the occurrence times.

    Under the modified Omori-Utsu model the transformed times, each event's
    expected count since the start of the interval, are the arrival times of a
    Poisson process of unit rate. Rescaled by the total expected count they are
    therefore uniform on (0, 1), and a Kolmogorov-Smirnov test against that
    uniform is the standard residual analysis for a point process. A small
    p-value means the sequence is not a single Omori decay, most often because a
    large aftershock has started a sequence of its own inside the window.

    The test is on times, which are continuous, so no discreteness correction is
    needed; the same test on binned magnitudes would reject on the binning alone.

    Parameters
    ----------
    times_days : array_like
        The elapsed times the fit was made on.
    fit : Omori
        The fitted parameters.
    t_end : float, optional
        End of the observation interval. Defaults to the largest time supplied.
    t_start : float, optional
        Start of the observation interval, 0 by default.

    Because c and p were estimated from the very times being tested, the fitted
    curve hugs the data and the statistic is systematically smaller than the
    standard Kolmogorov distribution assumes. Taking the textbook p-value here
    would be badly anti-conservative: on sequences drawn from the model it
    rejects at 5 per cent in 0 per cent of cases, with a mean p-value near 0.87
    rather than 0.5. The null distribution is therefore obtained by parametric
    bootstrap, simulating from the fitted model, refitting each replicate, and
    comparing statistics. Set `n_simulations=0` to fall back to the uncalibrated
    asymptotic form, which is reported as such.

    Parameters
    ----------
    times_days : array_like
        The elapsed times the fit was made on.
    fit : Omori
        The fitted parameters.
    t_end : float, optional
        End of the observation interval. Defaults to the largest time supplied.
    t_start : float, optional
        Start of the observation interval, 0 by default.
    n_simulations : int, optional
        Replicates used to calibrate the null distribution. Defaults to
        `constants.N_FIT_SIMULATIONS` (600). Zero selects the asymptotic form,
        whose p-value is far too generous here and must carry no verdict.
    seed : int, optional
        Seed of the simulation, so a reported p-value is reproducible.

    Returns
    -------
    FitTest
        The KS statistic, its p-value, the number of times tested, how the
        p-value was obtained, and the Monte Carlo standard error of the p-value
        itself. A p-value below 0.05 is the conventional signal that the model
        is inadequate, but read it against that standard error: a simulated
        p-value within two of them of the threshold decides nothing, and the
        remedy is more replicates rather than a firmer verdict.

    References
    ----------
    Ogata, Y. (1988). Lilliefors, H. W. (1967).
    """
    t = np.sort(np.asarray(times_days, float))
    t = t[t > t_start]
    if t.size == 0:
        return FitTest(
            float("nan"), float("nan"), 0, "no times after t_start", float("nan")
        )
    if t_end is None:
        t_end = float(t[-1])
    if not (t_end > t_start) or not np.isfinite(fit.c) or not np.isfinite(fit.p):
        return FitTest(
            float("nan"),
            float("nan"),
            int(t.size),
            "interval or fit unusable",
            float("nan"),
        )
    total = _integrated_rate(fit.c, fit.p, t_end, t_start)
    if not np.isfinite(total) or total <= 0:
        return FitTest(
            float("nan"),
            float("nan"),
            int(t.size),
            "interval or fit unusable",
            float("nan"),
        )
    if abs(fit.p - 1.0) < 1e-12:
        transformed = np.log((t + fit.c) / (t_start + fit.c))
    else:
        transformed = (
            (t + fit.c) ** (1.0 - fit.p) - (t_start + fit.c) ** (1.0 - fit.p)
        ) / (1.0 - fit.p)
    scaled = np.clip(transformed / total, 0.0, 1.0)
    observed = float(kstest(scaled, "uniform").statistic)

    if n_simulations is None:
        n_simulations = constants.N_FIT_SIMULATIONS
    if n_simulations < 1:
        asymptotic = float(kstest(scaled, "uniform").pvalue)
        return FitTest(
            observed, asymptotic, int(t.size), "asymptotic, uncalibrated", 0.0
        )

    rng = np.random.default_rng(seed)
    exceeded = 0
    used = 0
    for _ in range(n_simulations):
        drawn = np.sort(omori_sample(t.size, fit.p, fit.c, t_end, t_start, rng))
        try:
            refit = fit_omori(drawn, t_end=t_end, t_start=t_start)
        except ValueError:
            continue
        simulated = omori_fit_test(
            drawn, refit, t_end=t_end, t_start=t_start, n_simulations=0
        ).statistic
        if not np.isfinite(simulated):
            continue
        used += 1
        if simulated >= observed:
            exceeded += 1
    if used == 0:
        return FitTest(
            observed, float("nan"), int(t.size), "simulation failed", float("nan")
        )
    p_value = (1.0 + exceeded) / (1.0 + used)
    # The p-value is itself an estimate. Its binomial standard error decides
    # whether a verdict against a threshold means anything at this many
    # replicates, so it travels with the p-value rather than being implied.
    return FitTest(
        observed,
        p_value,
        int(t.size),
        f"parametric bootstrap, {used} replicates",
        float(np.sqrt(p_value * (1.0 - p_value) / used)),
    )
