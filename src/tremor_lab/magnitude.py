"""Magnitude-energy relations, the Bath energy screen, and homogenisation to Mw.

Every function accepts scalars or arrays and returns a NumPy value of the same shape.
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray

# log10(E) = ENERGY_A * M + ENERGY_B, with E in joules.
ENERGY_A = 1.5
ENERGY_B = 4.8

# Central value of the 1.1-1.2 range reported for the mainshock-to-largest-aftershock
# magnitude deficit.
DELTA_MB = 1.15


def energy_joules(mw: ArrayLike) -> NDArray[np.float64]:
    """
    Radiated seismic energy from moment magnitude.

    log10(E) = 1.5 Mw + 4.8, E in joules.

    Parameters
    ----------
    mw : array_like
        Moment magnitude.

    Returns
    -------
    ndarray
        Radiated energy in joules.

    References
    ----------
    Gutenberg, B. and Richter, C. F. (1956). Kanamori, H. (1977).
    """
    return 10.0 ** (ENERGY_A * np.asarray(mw, float) + ENERGY_B)


def bath_mag(mw_main: ArrayLike, delta_mb: float = DELTA_MB) -> NDArray[np.float64]:
    """
    Largest-aftershock magnitude expected under Bath's law.

    Parameters
    ----------
    mw_main : array_like
        Mainshock moment magnitude.
    delta_mb : float, optional
        Magnitude deficit, 1.15 by default.

    Returns
    -------
    ndarray
        Expected magnitude of the largest aftershock.

    References
    ----------
    Bath, M. (1965).
    """
    return np.asarray(mw_main, float) - delta_mb


def bath_ratio(
    mw_sec: ArrayLike, mw_main: ArrayLike, delta_mb: float = DELTA_MB
) -> NDArray[np.float64]:
    """
    Energy of a secondary event as a multiple of the Bath-expected largest aftershock.

    10 ** (1.5 (Mw_sec - (Mw_main - delta_mb))). A ratio near unity is the ordinary
    aftershock expectation; an order of magnitude or more marks a candidate under the
    energy screen.

    Parameters
    ----------
    mw_sec : array_like
        Moment magnitude of the secondary event.
    mw_main : array_like
        Mainshock moment magnitude.
    delta_mb : float, optional
        Bath magnitude deficit, 1.15 by default.

    Returns
    -------
    ndarray
        Energy ratio, dimensionless.

    References
    ----------
    Bath, M. (1965). Kanamori, H. (1977).
    """
    return 10.0 ** (
        ENERGY_A * (np.asarray(mw_sec, float) - bath_mag(mw_main, delta_mb))
    )


def ms_to_mw(ms: ArrayLike) -> NDArray[np.float64]:
    """
    Surface-wave magnitude placed on the moment scale.

    Mw = 0.67 Ms + 2.07 for Ms <= 6.1, and Mw = 0.99 Ms + 0.08 above it. The branches
    were calibrated on 3.0 <= Ms <= 6.1 and 6.2 <= Ms <= 8.2 and are applied here
    without a range check, matching the spreadsheet implementation this package is
    cross-validated against. They are mildly discontinuous across the uncalibrated
    6.1-6.2 gap: Ms 6.1 gives 6.157 on the lower branch and 6.119 on the upper.

    Parameters
    ----------
    ms : array_like
        Surface-wave magnitude.

    Returns
    -------
    ndarray
        Moment magnitude.

    References
    ----------
    Scordilis, E. M. (2006).
    """
    ms = np.asarray(ms, float)
    return np.where(ms <= 6.1, 0.67 * ms + 2.07, 0.99 * ms + 0.08)


def mb_to_mw(mb: ArrayLike) -> NDArray[np.float64]:
    """
    Body-wave magnitude placed on the moment scale.

    Mw = 0.85 mb + 1.03, calibrated on 3.5 <= mb <= 6.2 and applied here without a
    range check, for the reason given in `ms_to_mw`.

    Parameters
    ----------
    mb : array_like
        Body-wave magnitude.

    Returns
    -------
    ndarray
        Moment magnitude.

    References
    ----------
    Scordilis, E. M. (2006).
    """
    return 0.85 * np.asarray(mb, float) + 1.03


def to_mw(mag: ArrayLike, mtype: ArrayLike) -> NDArray[np.float64]:
    """
    Homogenise reported magnitudes to Mw according to their scale label.

    Labels beginning "mw" (mw, mwb, mww) pass through. "ms" and "mb" are converted by
    the Scordilis relations. Any other label, including ml, md and blanks, is taken as
    reported: no global conversion is published for those scales, so substituting one
    would be less faithful than leaving the catalogue value alone.

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
