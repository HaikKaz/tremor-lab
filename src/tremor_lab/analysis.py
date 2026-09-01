"""The standard per-case analysis: the estimators combined on one windowed
catalogue."""

from collections.abc import Mapping
from typing import Any

import pandas as pd

from tremor_lab import constants
from tremor_lab.bvalue import b_value_aki
from tremor_lab.completeness import fmd, mc_maxcurvature
from tremor_lab.magnitude import bath_mag, energy_joules
from tremor_lab.omori import bootstrap_omori, fit_omori

REQUIRED_COLUMNS = ("dt_days", "mw")


def analyze_case(
    catalog_df: pd.DataFrame,
    mainshock: Mapping,
    mc_threshold: float | None = None,
    window_days: float | None = None,
    dm: float | None = None,
    mc_correction: float | None = None,
    n_boot: int | None = None,
    seed: int = 0,
    min_events_for_mc: int | None = None,
    min_events_for_fit: int | None = None,
) -> dict[str, Any]:
    """
    Completeness, b-value, decay and energy for one aftershock sequence.

    The b-value and the decay fit are attempted only where the catalogue above the
    threshold is dense enough to support them. Where it is not, they are returned as
    None with the reason under "note", rather than a fit being forced on a sparse
    sequence and reported as though it were stable.

    Parameters
    ----------
    catalog_df : DataFrame
        Must carry dt_days (days since the mainshock) and mw. This is the shape
        `tremor_lab.catalog.read_catalog` produces.
    mainshock : mapping
        Needs mw; lat, lon and t are not used here, having already been applied when
        the catalogue was windowed.
    mc_threshold : float, optional
        Fixed magnitude threshold for the b-value and the decay fit. When omitted the
        estimated completeness magnitude is used.
    window_days : float, optional
        Aftershock window in days. Defaults to `constants.WINDOW_DAYS` (180).
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM`.
    mc_correction : float, optional
        Maximum-curvature correction. Defaults to `constants.MC_CORRECTION`.
    n_boot : int, optional
        Bootstrap resamples for the decay uncertainty. Defaults to `constants.N_BOOT`.
        Pass 0 to skip the bootstrap.
    seed : int, optional
        Seed of the bootstrap, so a reported uncertainty is reproducible.
    min_events_for_mc, min_events_for_fit : int, optional
        Density rules described above. Default to `constants.MIN_EVENTS_FOR_MC` (50)
        and `constants.MIN_EVENTS_FOR_FIT` (100).

    Returns
    -------
    dict
        Keys n_events, mc, threshold, n_above, b_value, omori, omori_bootstrap,
        mainshock_energy_j, bath_mag, fmd and note.

    Raises
    ------
    KeyError
        If the catalogue lacks the columns named in `REQUIRED_COLUMNS`.
    """
    missing = [c for c in REQUIRED_COLUMNS if c not in catalog_df.columns]
    if missing:
        raise KeyError(
            f"catalogue is missing {missing}; expected the columns produced by "
            f"tremor_lab.catalog.read_catalog"
        )
    window_days = constants.WINDOW_DAYS if window_days is None else window_days
    n_boot = constants.N_BOOT if n_boot is None else n_boot
    if min_events_for_mc is None:
        min_events_for_mc = constants.MIN_EVENTS_FOR_MC
    if min_events_for_fit is None:
        min_events_for_fit = constants.MIN_EVENTS_FOR_FIT

    post = catalog_df[
        (catalog_df["dt_days"] > 0) & (catalog_df["dt_days"] <= window_days)
    ]
    mags = post["mw"].to_numpy()

    mc = (
        mc_maxcurvature(mags, dm=dm, correction=mc_correction)
        if len(post) > min_events_for_mc
        else None
    )
    threshold = mc if mc_threshold is None else mc_threshold
    above = post if threshold is None else post[post["mw"] >= threshold]

    result: dict[str, Any] = {
        "n_events": len(post),
        "mc": mc,
        "threshold": threshold,
        "n_above": len(above) if threshold is not None else None,
        "b_value": None,
        "omori": None,
        "omori_bootstrap": None,
        "mainshock_energy_j": float(energy_joules(mainshock["mw"])),
        "bath_mag": float(bath_mag(mainshock["mw"])),
        "fmd": fmd(mags, dm=dm) if len(post) else None,
        "note": None,
    }

    if threshold is None:
        result["note"] = (
            f"fewer than {min_events_for_mc} events in the window; "
            f"no completeness magnitude estimated"
        )
        return result
    if len(above) <= min_events_for_fit:
        result["note"] = (
            f"only {len(above)} events at or above M {threshold}; "
            f"fewer than the {min_events_for_fit} required for a stable fit"
        )
        return result

    result["b_value"] = b_value_aki(above["mw"].to_numpy(), threshold, dm=dm)
    times = above["dt_days"].to_numpy()
    result["omori"] = fit_omori(times)
    if n_boot:
        result["omori_bootstrap"] = bootstrap_omori(times, n_boot=n_boot, seed=seed)
    return result
