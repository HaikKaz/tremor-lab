"""Comparing magnitudes against the decimal grid they are reported on.

Reported magnitudes are decimal: 3.4, 3.5, 4.1. None of those is exactly
representable in binary floating point, and neither is a bin edge computed from
them. Asking `m >= mc - dm / 2` directly therefore decides membership by rounding
error. At the default bin width of 0.1 the errors fell the harmless way on the
reference catalogue, which is why this went unnoticed; at 0.2, a documented and
supported setting, they did not. `4.2 - 0.2 / 2` evaluates to 4.1000000000000005,
which every event reported at exactly M 4.1 fails - dropping a whole magnitude
class out of its own completeness sample and moving b by six times its own
standard error.

Both routines here take the same tolerance, `constants.GRID_TOLERANCE`, so a
magnitude a reader would call equal to an edge is treated as equal everywhere in
the package. It is a published default like any other and can be changed by
keyword, by reassignment, or from a settings file.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from tremor_lab import constants

__all__ = ["at_or_above", "bin_index"]


def at_or_above(
    mags: ArrayLike, threshold: float, tolerance: float | None = None
) -> NDArray[np.bool_]:
    """
    Which magnitudes sit at or above a threshold, treating the grid as exact.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    threshold : float
        The magnitude to compare against, typically a bin edge such as mc - dm/2.
    tolerance : float, optional
        How close counts as equal. Defaults to `constants.GRID_TOLERANCE`.

    Returns
    -------
    ndarray of bool
        True where the magnitude is at or above the threshold.
    """
    tolerance = constants.GRID_TOLERANCE if tolerance is None else tolerance
    return np.asarray(mags, float) >= threshold - tolerance


def bin_index(
    mags: ArrayLike, lo: float, dm: float, tolerance: float | None = None
) -> NDArray[np.int64]:
    """
    The bin each magnitude falls in, counting from the bin centred on `lo`.

    Magnitudes go to the nearest bin centre, and a magnitude exactly half way
    between two centres goes up. That is the spreadsheet's convention and the
    browser page's, so all three implementations agree on which bin an event
    belongs to.

    The tolerance is what makes that convention true. Without it, a magnitude
    mathematically half way between two centres lands on whichever side the
    division happened to round to: at dm 0.2, M 3.5 went up and M 3.3 went down,
    so a catalogue holding one event at every magnitude produced bins of 1, 3, 1,
    2, 3 rather than 2 throughout - and handed the maximum-curvature method a
    mode that was an artefact of arithmetic.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    lo : float
        Centre of the lowest bin.
    dm : float
        Bin width.
    tolerance : float, optional
        How close to the half-way point counts as being on it. Defaults to
        `constants.GRID_TOLERANCE`.

    Returns
    -------
    ndarray of int
        Zero-based bin index for each magnitude.
    """
    tolerance = constants.GRID_TOLERANCE if tolerance is None else tolerance
    offset = (np.asarray(mags, float) - lo) / dm
    return np.floor(offset + 0.5 + tolerance).astype(np.int64)
