"""Gutenberg-Richter b-value and its uncertainty."""

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike

from tremor_lab import constants


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
    m = m[m >= threshold]
    n = m.size
    if n < 2:
        raise ValueError(
            f"b-value needs at least two events at or above {threshold}; got {n}"
        )
    mean_m = m.mean()
    b = 1.0 / (np.log(10) * (mean_m - threshold))
    sigma = shi_bolt_k * b**2 * np.sqrt(((m - mean_m) ** 2).sum() / (n * (n - 1)))
    return BValue(float(b), float(sigma), n)
