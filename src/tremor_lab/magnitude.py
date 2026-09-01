"""Magnitude-energy relations, the Bath energy screen, and homogenisation to Mw.

Every function accepts scalars or arrays and returns a NumPy value of the same shape.
Coefficients default to the published values in `tremor_lab.constants` and can be
overridden per call or, by reassigning them there, for a whole session.
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from tremor_lab import constants


def energy_joules(
    mw: ArrayLike, a: float | None = None, b: float | None = None
) -> NDArray[np.float64]:
    """
    Radiated seismic energy from moment magnitude.

    log10(E) = a Mw + b, E in joules.

    Parameters
    ----------
    mw : array_like
        Moment magnitude.
    a, b : float, optional
        Slope and intercept of the energy-magnitude relation. Default to
        `constants.ENERGY_A` (1.5) and `constants.ENERGY_B` (4.8).

    Returns
    -------
    ndarray
        Radiated energy in joules.

    References
    ----------
    Gutenberg, B. and Richter, C. F. (1956). Kanamori, H. (1977).
    """
    a = constants.ENERGY_A if a is None else a
    b = constants.ENERGY_B if b is None else b
    return 10.0 ** (a * np.asarray(mw, float) + b)


def bath_mag(mw_main: ArrayLike, delta_mb: float | None = None) -> NDArray[np.float64]:
    """
    Largest-aftershock magnitude expected under Bath's law.

    Parameters
    ----------
    mw_main : array_like
        Mainshock moment magnitude.
    delta_mb : float, optional
        Magnitude deficit. Defaults to `constants.DELTA_MB` (1.15); the reported
        plausible range is 1.1 to 1.2.

    Returns
    -------
    ndarray
        Expected magnitude of the largest aftershock.

    References
    ----------
    Bath, M. (1965).
    """
    delta_mb = constants.DELTA_MB if delta_mb is None else delta_mb
    return np.asarray(mw_main, float) - delta_mb


def bath_ratio(
    mw_sec: ArrayLike,
    mw_main: ArrayLike,
    delta_mb: float | None = None,
    a: float | None = None,
) -> NDArray[np.float64]:
    """
    Energy of a secondary event as a multiple of the Bath-expected largest aftershock.

    10 ** (a (Mw_sec - (Mw_main - delta_mb))). A ratio near unity is the ordinary
    aftershock expectation; an order of magnitude or more marks a candidate under the
    energy screen.

    Parameters
    ----------
    mw_sec : array_like
        Moment magnitude of the secondary event.
    mw_main : array_like
        Mainshock moment magnitude.
    delta_mb : float, optional
        Bath magnitude deficit. Defaults to `constants.DELTA_MB`.
    a : float, optional
        Slope of the energy-magnitude relation. Defaults to `constants.ENERGY_A`.

    Returns
    -------
    ndarray
        Energy ratio, dimensionless.

    References
    ----------
    Bath, M. (1965). Kanamori, H. (1977).
    """
    a = constants.ENERGY_A if a is None else a
    return 10.0 ** (a * (np.asarray(mw_sec, float) - bath_mag(mw_main, delta_mb)))


def ms_to_mw(
    ms: ArrayLike,
    branch: float | None = None,
    low_slope: float | None = None,
    low_intercept: float | None = None,
    high_slope: float | None = None,
    high_intercept: float | None = None,
) -> NDArray[np.float64]:
    """
    Surface-wave magnitude placed on the moment scale.

    Mw = 0.67 Ms + 2.07 at or below Ms 6.1, and Mw = 0.99 Ms + 0.08 above it. The
    branches were calibrated on 3.0 <= Ms <= 6.1 and 6.2 <= Ms <= 8.2 and are applied
    here without a range check, matching the spreadsheet implementation this package is
    cross-validated against. They are mildly discontinuous across the uncalibrated
    6.1-6.2 gap: Ms 6.1 gives 6.157 on the lower branch and 6.119 on the upper. Supply
    the coefficients to substitute a regional relation.

    Parameters
    ----------
    ms : array_like
        Surface-wave magnitude.
    branch : float, optional
        Magnitude at which the relation switches branch. Defaults to
        `constants.MS_MW_BRANCH`.
    low_slope, low_intercept : float, optional
        Coefficients at or below `branch`. Default to `constants.MS_MW_LOW_SLOPE` and
        `constants.MS_MW_LOW_INTERCEPT`.
    high_slope, high_intercept : float, optional
        Coefficients above `branch`. Default to `constants.MS_MW_HIGH_SLOPE` and
        `constants.MS_MW_HIGH_INTERCEPT`.

    Returns
    -------
    ndarray
        Moment magnitude.

    References
    ----------
    Scordilis, E. M. (2006).
    """
    branch = constants.MS_MW_BRANCH if branch is None else branch
    low_slope = constants.MS_MW_LOW_SLOPE if low_slope is None else low_slope
    high_slope = constants.MS_MW_HIGH_SLOPE if high_slope is None else high_slope
    if low_intercept is None:
        low_intercept = constants.MS_MW_LOW_INTERCEPT
    if high_intercept is None:
        high_intercept = constants.MS_MW_HIGH_INTERCEPT
    ms = np.asarray(ms, float)
    return np.where(
        ms <= branch, low_slope * ms + low_intercept, high_slope * ms + high_intercept
    )


def mb_to_mw(
    mb: ArrayLike, slope: float | None = None, intercept: float | None = None
) -> NDArray[np.float64]:
    """
    Body-wave magnitude placed on the moment scale.

    Mw = 0.85 mb + 1.03, calibrated on 3.5 <= mb <= 6.2 and applied here without a
    range check, for the reason given in `ms_to_mw`.

    Parameters
    ----------
    mb : array_like
        Body-wave magnitude.
    slope, intercept : float, optional
        Coefficients of the relation. Default to `constants.MB_MW_SLOPE` and
        `constants.MB_MW_INTERCEPT`.

    Returns
    -------
    ndarray
        Moment magnitude.

    References
    ----------
    Scordilis, E. M. (2006).
    """
    slope = constants.MB_MW_SLOPE if slope is None else slope
    intercept = constants.MB_MW_INTERCEPT if intercept is None else intercept
    return slope * np.asarray(mb, float) + intercept


def to_mw(mag: ArrayLike, mtype: ArrayLike) -> NDArray[np.float64]:
    """
    Homogenise reported magnitudes to Mw according to their scale label.

    Labels beginning "mw" (mw, mwb, mww) pass through. "ms" and "mb" are converted by
    `ms_to_mw` and `mb_to_mw`, which read their coefficients from `tremor_lab.constants`
    at call time, so changing a coefficient there changes this function too. Any other
    label, including ml, md and blanks, is taken as reported: no global relation is
    published for those scales, so substituting one would be less faithful than leaving
    the catalogue value alone.

    Parameters
    ----------
    mag : array_like
        Reported magnitude.
    mtype : array_like of str
        Scale label, matched case-insensitively after stripping whitespace.

    Returns
    -------
    ndarray
        Magnitude on the moment scale.

    References
    ----------
    Scordilis, E. M. (2006).
    """
    mag = np.asarray(mag, float)
    label = np.strings.lower(np.strings.strip(np.asarray(mtype, dtype=np.str_)))
    return np.select(
        [np.strings.startswith(label, "mw"), label == "ms", label == "mb"],
        [mag, ms_to_mw(mag), mb_to_mw(mag)],
        default=mag,
    )
