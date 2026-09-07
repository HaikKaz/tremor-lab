"""Gutenberg-Richter b-value and its uncertainty."""

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from tremor_lab import constants
from tremor_lab.completeness import fmd
from tremor_lab.grid import at_or_above


class BValue(NamedTuple):
    """b-value, its standard error, and the number of events it was estimated from."""

    b: float
    sigma: float
    n: int


def b_value_aki(
    mags: ArrayLike,
    mc: float,
    dm: float | None = None,
    shi_bolt_k: float | None = None,
) -> BValue:
    """
    b-value by the Aki (1965) maximum-likelihood estimator.

    b = 1 / (ln10 (mean(M) - (Mc - dm/2)))

    with the Shi and Bolt (1982) standard error

    sigma = k b^2 sqrt(sum (M - mean(M))^2 / (n (n - 1))).

    Events below Mc - dm/2 are discarded. The half-bin offset is the lower edge of the
    completeness bin: a magnitude reported as Mc stands for the interval Mc +/- dm/2,
    so the smallest complete magnitude is Mc - dm/2 rather than Mc.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    mc : float
        Completeness magnitude.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM` (0.1). Pass 0.0 for unbinned
        magnitudes, where no half-bin offset applies.
    shi_bolt_k : float, optional
        Coefficient of the standard error. Defaults to `constants.SHI_BOLT_K` (2.30),
        the published rounding of ln(10).

    Returns
    -------
    BValue
        Named tuple of b, sigma and n.

    Raises
    ------
    ValueError
        If fewer than two events lie at or above the completeness threshold, where the
        standard error is undefined.

    References
    ----------
    Aki, K. (1965). Shi, Y. and Bolt, B. A. (1982).
    """
    dm = constants.DM if dm is None else dm
    shi_bolt_k = constants.SHI_BOLT_K if shi_bolt_k is None else shi_bolt_k
    threshold = mc - dm / 2
    m = np.asarray(mags, float)
    m = m[at_or_above(m, threshold)]
    n = m.size
    if n < 2:
        raise ValueError(
            f"b-value needs at least two events at or above {threshold}; got {n}"
        )
    mean_m = m.mean()
    b = 1.0 / (np.log(10) * (mean_m - threshold))
    sigma = shi_bolt_k * b**2 * np.sqrt(((m - mean_m) ** 2).sum() / (n * (n - 1)))
    return BValue(float(b), float(sigma), n)


class BStability(NamedTuple):
    """b and its standard error as a function of the threshold applied."""

    thresholds: NDArray[np.float64]
    b: NDArray[np.float64]
    sigma: NDArray[np.float64]
    n: NDArray[np.int64]


def b_stability(
    mags: ArrayLike,
    thresholds: ArrayLike | None = None,
    dm: float | None = None,
    min_events: int = 50,
) -> BStability:
    """
    b against the threshold it was estimated at.

    Above a correctly estimated completeness magnitude the b-value should not
    depend on where the threshold is put; a curve that keeps climbing is the
    standard sign that completeness has been placed too low, or that the sample
    is not a single Gutenberg-Richter population. The Shi and Bolt error is the
    scatter at one threshold and says nothing about this, so the curve carries
    information the headline figure cannot.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    thresholds : array_like, optional
        Thresholds to evaluate. Defaults to a grid of 2.5 magnitude units
        starting at the mode of the incremental distribution.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM`.
    min_events : int, optional
        Thresholds retaining fewer events than this are dropped, since b is not
        meaningful there. 50 by default.

    Returns
    -------
    BStability
        Arrays of threshold, b, sigma and n, one entry per usable threshold.

    References
    ----------
    Cao, A. and Gao, S. S. (2002). Woessner, J. and Wiemer, S. (2005).
    """
    dm = constants.DM if dm is None else dm
    m = np.asarray(mags, float)
    if thresholds is None:
        # Anchored at the mode of the incremental distribution, not at the
        # smallest magnitude: completeness cannot lie below the mode, which is
        # the premise of the maximum-curvature method, and one anomalously small
        # event moves the minimum while leaving the mode where it was.
        edges, inc, _ = fmd(m, dm)
        start = float(edges[int(np.argmax(inc))])
        thresholds = np.round(start + dm * np.arange(round(2.5 / dm) + 1), 10)
    thresholds = np.atleast_1d(np.asarray(thresholds, float))

    kept, bs, sigmas, ns = [], [], [], []
    for threshold in thresholds:
        if int(at_or_above(m, threshold - dm / 2).sum()) < min_events:
            continue
        estimate = b_value_aki(m, threshold, dm=dm)
        kept.append(threshold)
        bs.append(estimate.b)
        sigmas.append(estimate.sigma)
        ns.append(estimate.n)
    return BStability(
        np.array(kept), np.array(bs), np.array(sigmas), np.array(ns, dtype=np.int64)
    )


def mc_b_stability(
    mags: ArrayLike,
    dm: float | None = None,
    span: float = 0.5,
    min_events: int = 50,
) -> float | None:
    """
    Completeness magnitude by the b-value stability method.

    The lowest threshold at which b has stopped changing: b there differs from
    the mean of b over the next `span` magnitude units by no more than its own
    standard error. Where maximum curvature answers "where is the peak of the
    incremental distribution", this answers "from where onward does the slope
    stop moving", and the two disagreeing is itself worth reporting.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM`.
    span : float, optional
        Width of the averaging window above each candidate, 0.5 by default.
    min_events : int, optional
        Minimum events above a threshold for it to be a candidate.

    Returns
    -------
    float or None
        The completeness magnitude, or None if b never stabilises over the
        range available.

    References
    ----------
    Cao, A. and Gao, S. S. (2002). Woessner, J. and Wiemer, S. (2005).
    """
    dm = constants.DM if dm is None else dm
    # How many bins the averaging window spans, rounded the same way magnitudes
    # are binned: halves up, on the grid. Python's round() is banker's rounding
    # and sends 2.5 down to 2 where JavaScript's Math.round sends it up to 3, so
    # at dm 0.2 the package and the browser page averaged over different windows
    # and reported different completeness magnitudes, 4.0 against 4.2, from
    # identical b-stability curves. Where span is not a whole number of bins the
    # window is the next bin up: 0.5 at dm 0.2 averages over 0.6.
    steps = int(np.floor(span / dm + 0.5 + constants.GRID_TOLERANCE))
    if steps < 1:
        raise ValueError(
            f"a bin width of {dm} rounds the averaging window of {span} down "
            f"to no bins at all, so there is nothing to average over and "
            f"stability cannot be tested; use a narrower dm or a wider span"
        )
    curve = b_stability(mags, dm=dm, min_events=min_events)
    for i in range(len(curve.thresholds) - steps):
        window = curve.b[i : i + steps + 1]
        if abs(window.mean() - curve.b[i]) <= curve.sigma[i]:
            return float(curve.thresholds[i])
    return None


def b_value_tinti(
    mags: ArrayLike,
    mc: float,
    dm: float | None = None,
    shi_bolt_k: float | None = None,
) -> BValue:
    """
    b-value by the exact maximum-likelihood estimator for binned magnitudes.

    b = ln(1 + dm / mean(M - Mc)) / (dm ln10)

    Where `b_value_aki` applies Utsu's half-bin offset to the Aki estimator,
    which is a first-order approximation, this is the exact solution for
    magnitudes reported on a grid of width dm. The two converge as dm shrinks:
    on the reference catalogue they differ by 0.32 per cent at dm 0.1 and by
    0.004 per cent at dm 0.01.

    This is the estimator the independent package `seismostats` uses, and
    `examples/compare_with_seismostats.py` checks that this implementation
    reproduces it. The package's own reported values use `b_value_aki`, which is
    the convention of the thesis and of the spreadsheet implementation; this
    function exists so the choice can be tested rather than assumed.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    mc : float
        Completeness magnitude.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM` (0.1).
    shi_bolt_k : float, optional
        Coefficient of the standard error. Defaults to `constants.SHI_BOLT_K`.

    Returns
    -------
    BValue
        Named tuple of b, sigma and n.

    Raises
    ------
    ValueError
        If fewer than two events lie at or above the completeness threshold.

    References
    ----------
    Tinti, S. and Mulargia, F. (1987). Shi, Y. and Bolt, B. A. (1982).
    """
    dm = constants.DM if dm is None else dm
    shi_bolt_k = constants.SHI_BOLT_K if shi_bolt_k is None else shi_bolt_k
    threshold = mc - dm / 2
    m = np.asarray(mags, float)
    m = m[at_or_above(m, threshold)]
    n = m.size
    if n < 2:
        raise ValueError(
            f"b-value needs at least two events at or above {threshold}; got {n}"
        )
    mean_above = (m - mc).mean()
    if mean_above <= 0:
        raise ValueError("mean magnitude does not exceed the completeness magnitude")
    b = np.log1p(dm / mean_above) / (dm * np.log(10))
    mean_m = m.mean()
    sigma = shi_bolt_k * b**2 * np.sqrt(((m - mean_m) ** 2).sum() / (n * (n - 1)))
    return BValue(float(b), float(sigma), n)
