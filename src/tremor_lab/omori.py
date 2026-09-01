"""Modified Omori-Utsu aftershock decay, fitted by maximum likelihood."""

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import minimize

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
    t = t[t > 0]
    if t.size < 2:
        raise ValueError(f"Omori fit needs at least two positive times; got {t.size}")
    if t_end is None:
        t_end = float(t[-1])
    c0 = constants.OMORI_C0 if c0 is None else c0
    p0 = constants.OMORI_P0 if p0 is None else p0

    search = minimize(
        _profiled_nll,
        (c0, p0),
        args=(t, t_end),
        method="Nelder-Mead",
        # fatol is absolute, so it must scale with the likelihood, which grows
        # with n; a fixed 1e-10 is unreachable on a large catalogue and the fit
        # is then rejected as unconverged despite sitting on the optimum.
        options={
            "xatol": 1e-8,
            "fatol": max(1e-10, 1e-10 * abs(_profiled_nll((c0, p0), t, t_end))),
            "maxiter": 5000,
        },
    )
    if not search.success:
        raise ValueError(f"Omori fit did not converge: {search.message}")
    c, p = search.x
    return Omori(float(p), float(c), t.size / _integrated_rate(c, p, t_end), t.size)


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
    t = np.asarray(times_days, float)
    t = t[t > 0]
    if t_end is None:
        t_end = float(t.max())
    rng = np.random.default_rng(seed)
    fits = [
        fit_omori(rng.choice(t, size=t.size, replace=True), t_end=t_end)
        for _ in range(n_boot)
    ]
    return OmoriBootstrap(
        float(np.std([f.p for f in fits])),
        float(np.std([f.c for f in fits])),
        n_boot,
    )


def _integrated_rate(c: float, p: float, t_end: float) -> float:
    """Integral of (c + t)^-p over (0, t_end], the expected count per unit k."""
    if abs(p - 1.0) < 1e-12:
        return float(np.log((t_end + c) / c))
    return float(((t_end + c) ** (1.0 - p) - c ** (1.0 - p)) / (1.0 - p))


def _profiled_nll(cp: tuple[float, float], t: np.ndarray, t_end: float) -> float:
    """Negative log-likelihood with k replaced by its maximising value n / integral."""
    c, p = cp
    if c <= 0 or p <= 0:
        return np.inf
    integral = _integrated_rate(c, p, t_end)
    if integral <= 0:
        return np.inf
    n = t.size
    return float(p * np.sum(np.log(t + c)) - n * np.log(n / integral) + n)
