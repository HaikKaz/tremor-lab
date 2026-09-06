"""Frequency-magnitude distribution and the completeness magnitude."""

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from tremor_lab import constants


class FMD(NamedTuple):
    """Binned frequency-magnitude distribution."""

    edges: NDArray[np.float64]
    inc: NDArray[np.int64]
    cum: NDArray[np.int64]


def fmd(mags: ArrayLike, dm: float | None = None) -> FMD:
    """
    Incremental and cumulative event counts per magnitude bin.

    Each magnitude is assigned to the nearest bin centre rather than to a half-open
    interval. Reported magnitudes arrive on an exact decimal grid that binary floating
    point cannot represent, so a value sitting on a bin boundary is otherwise liable to
    fall on the wrong side of it. The cumulative count at a bin is the number of events
    at or above it, which is the quantity plotted against magnitude on a
    Gutenberg-Richter diagram.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    dm : float, optional
        Bin width. Defaults to `constants.DM` (0.1).

    Returns
    -------
    FMD
        Named tuple of bin centres, incremental counts and cumulative counts.

    References
    ----------
    Gutenberg, B. and Richter, C. F. (1944).
    """
    dm = constants.DM if dm is None else dm
    m = np.asarray(mags, float)
    lo = np.floor(m.min() / dm) * dm
    hi = np.ceil(m.max() / dm) * dm
    n_bins = int(np.rint((hi - lo) / dm)) + 1
    edges = np.round(lo + dm * np.arange(n_bins), _edge_decimals(dm))
    inc = np.bincount(np.rint((m - lo) / dm).astype(int), minlength=n_bins)
    return FMD(edges, inc, np.cumsum(inc[::-1])[::-1])


def mc_maxcurvature(
    mags: ArrayLike, dm: float | None = None, correction: float | None = None
) -> float:
    """
    Completeness magnitude by the maximum-curvature method.

    The point of maximum curvature of the cumulative frequency-magnitude distribution
    is the mode of the incremental distribution. The method is known to underestimate
    completeness, so a positive correction is added; +0.2 is the value used throughout
    this work and by the spreadsheet implementation.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    dm : float, optional
        Bin width. Defaults to `constants.DM` (0.1).
    correction : float, optional
        Added to the modal magnitude. Defaults to `constants.MC_CORRECTION` (0.2).
        Pass 0.0 for the uncorrected maximum-curvature estimate.

    Returns
    -------
    float
        Completeness magnitude, rounded to the resolution of the bin width.

    References
    ----------
    Wiemer, S. and Wyss, M. (2000).
    """
    dm = constants.DM if dm is None else dm
    correction = constants.MC_CORRECTION if correction is None else correction
    edges, inc, _ = fmd(mags, dm)
    return float(np.round(edges[np.argmax(inc)] + correction, _edge_decimals(dm)))


def _edge_decimals(dm: float) -> int:
    """Decimals that render a bin label as 3.2 rather than 3.2000000000000004."""
    return int(np.ceil(-np.log10(dm))) + 1


class GoodnessOfFit(NamedTuple):
    """Result of the frequency-magnitude goodness-of-fit test."""

    mc: float | None
    r_value: float | None
    candidates: NDArray[np.float64]
    r_values: NDArray[np.float64]


def mc_goodness_of_fit(
    mags: ArrayLike,
    dm: float | None = None,
    confidence: float = 90.0,
    max_candidates: int = 30,
) -> GoodnessOfFit:
    """
    Completeness magnitude by goodness of fit to a Gutenberg-Richter law.

    For each candidate completeness magnitude a b-value is estimated, the
    cumulative distribution it implies is synthesised, and the residual between
    observed and synthetic counts is expressed as a percentage of the observed
    total. The lowest candidate whose fit explains at least `confidence` per
    cent of the observation is returned. Unlike maximum curvature this asks
    whether the data above the threshold actually look like a power law, rather
    than only where the incremental distribution peaks.

    The test is applied to binned counts rather than by a distributional test on
    the magnitudes themselves, because reported magnitudes are discrete and a
    continuous test rejects on the discreteness alone once the sample is large.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM`.
    confidence : float, optional
        Percentage of the observed distribution the fit must explain, 90 by
        default; 95 is the stricter level also in common use.
    max_candidates : int, optional
        How many bins upward to try before giving up.

    Returns
    -------
    GoodnessOfFit
        The chosen completeness magnitude and its R value, or None for both
        where no candidate reaches the level, together with every candidate
        tried and the R value it achieved.

    References
    ----------
    Wiemer, S. and Wyss, M. (2000).
    """
    dm = constants.DM if dm is None else dm
    m = np.asarray(mags, float)
    edges, _, cum = fmd(m, dm)
    chosen_mc, chosen_r = None, None
    candidates, r_values = [], []

    for i, mc in enumerate(edges[:max_candidates]):
        above = m[m >= mc - dm / 2]
        if above.size < 2:
            continue
        mean_m = above.mean()
        if mean_m <= mc - dm / 2:
            continue
        b = 1.0 / (np.log(10) * (mean_m - (mc - dm / 2)))
        observed = cum[i:]
        if observed[0] <= 0:
            continue
        # Synthetic cumulative counts for the fitted law, anchored at the
        # candidate so both distributions start from the same total.
        synthetic = observed[0] * 10.0 ** (-b * (edges[i:] - mc))
        r = 100.0 * (1.0 - np.abs(observed - synthetic).sum() / observed.sum())
        candidates.append(float(mc))
        r_values.append(float(r))
        if chosen_mc is None and r >= confidence:
            chosen_mc, chosen_r = float(mc), float(r)

    return GoodnessOfFit(chosen_mc, chosen_r, np.array(candidates), np.array(r_values))
