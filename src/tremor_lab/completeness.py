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
