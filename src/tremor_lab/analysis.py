"""The standard per-case analysis: the estimators combined on one windowed
catalogue."""

from collections.abc import Mapping
from typing import Any

import pandas as pd

from tremor_lab import constants
from tremor_lab.bvalue import b_stability, b_value_aki, mc_b_stability
from tremor_lab.completeness import fmd, mc_goodness_of_fit, mc_maxcurvature
from tremor_lab.grid import at_or_above
from tremor_lab.magnitude import bath_mag, energy_joules
from tremor_lab.omori import bootstrap_omori, fit_omori, omori_fit_test

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
    n_fit_simulations: int | None = None,
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
    n_fit_simulations : int, optional
        Replicates used to calibrate the decay fit test, the slowest thing here.
        Defaults to `constants.N_FIT_SIMULATIONS` (600). Zero falls back to the
        uncalibrated asymptotic p-value, which is reported without a verdict.

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

    usable = catalog_df.dropna(subset=list(REQUIRED_COLUMNS))
    n_unusable = len(catalog_df) - len(usable)
    post = usable[(usable["dt_days"] > 0) & (usable["dt_days"] <= window_days)]
    mags = post["mw"].to_numpy()

    dense_enough = len(post) > min_events_for_mc
    mc = (
        mc_maxcurvature(mags, dm=dm, correction=mc_correction) if dense_enough else None
    )
    # Maximum curvature answers only "where does the incremental distribution
    # peak". Two further methods ask whether the data above a threshold actually
    # look like a Gutenberg-Richter law, and where the slope stops moving. Where
    # the three disagree, the spread is the honest uncertainty in completeness,
    # and it is usually far larger than the b-value's own standard error.
    mc_methods = {"maximum_curvature": mc}
    stability = None
    if dense_enough:
        mc_methods["goodness_of_fit"] = mc_goodness_of_fit(mags, dm=dm).mc
        try:
            mc_methods["b_stability"] = mc_b_stability(mags, dm=dm)
        except ValueError as reason:
            # b-stability is one of three completeness estimates and the only one
            # that can refuse: at a bin width wider than twice the averaging
            # window there is nothing to average over. That is a reason to report
            # no answer from this method, not to abandon the whole analysis - the
            # browser page has always reported the rest, and the two disagreed
            # about whether such a run is possible at all.
            mc_methods["b_stability"] = None
            mc_methods["b_stability_note"] = str(reason)
        stability = b_stability(mags, dm=dm)
    threshold = mc if mc_threshold is None else mc_threshold
    bin_width = constants.DM if dm is None else dm
    # At or above the LOWER EDGE of the threshold's bin, not the threshold itself.
    # A magnitude reported as Mc stands for the interval Mc +/- dm/2, which is the
    # whole premise of the half-bin offset that b_value_aki then puts in its
    # denominator. Selecting at Mc while computing with Mc - dm/2 drops the lower
    # half of the completeness class from a sample the formula assumes contains it.
    # On a catalogue reported on the dm grid the two select the same events, which
    # is why this went unseen; on a homogenised one they do not, because Ms -> Mw
    # maps a 0.1 grid onto a 0.067 grid, and b came out several standard errors low.
    above = (
        post
        if threshold is None
        else post[at_or_above(post["mw"].to_numpy(), threshold - bin_width / 2)]
    )

    result: dict[str, Any] = {
        "n_events": len(post),
        "mc": mc,
        "threshold": threshold,
        "n_above": len(above) if threshold is not None else None,
        # The bound the sample was actually taken at, so the report can say it
        # rather than leave a reader to infer it from the threshold.
        "sample_floor": (
            None if threshold is None else round(threshold - bin_width / 2, 10)
        ),
        "b_value": None,
        "omori": None,
        "omori_bootstrap": None,
        "mainshock_energy_j": float(energy_joules(mainshock["mw"])),
        "bath_mag": float(bath_mag(mainshock["mw"])),
        "fmd": fmd(mags, dm=dm) if len(post) else None,
        "mc_methods": mc_methods,
        "b_stability": stability,
        "omori_fit_test": None,
        "n_unusable": n_unusable,
        "omori_warning": None,
        "note": None,
    }

    if threshold is None:
        result["note"] = (
            f"{len(post)} events in the window; estimating a completeness "
            f"magnitude needs more than {min_events_for_mc}"
        )
        return result
    if len(above) <= min_events_for_fit:
        result["note"] = (
            f"only {len(above)} events at or above M {threshold}; a stable fit "
            f"needs more than {min_events_for_fit}"
        )
        return result

    result["b_value"] = b_value_aki(above["mw"].to_numpy(), threshold, dm=dm)
    times = above["dt_days"].to_numpy()
    result["omori"] = fit_omori(times)
    result["omori_warning"] = _omori_caution(result["omori"])
    result["omori_fit_test"] = omori_fit_test(
        times, result["omori"], n_simulations=n_fit_simulations, seed=seed
    )
    if n_boot:
        result["omori_bootstrap"] = bootstrap_omori(times, n_boot=n_boot, seed=seed)
    return result


def _omori_caution(fit) -> str | None:
    """
    Say when a converged Omori fit should not be believed.

    Two cases. The offset c can collapse onto zero, the boundary of the model,
    where the likelihood has no interior maximum: that happens when the earliest
    events sit essentially on the origin time, and c is then not identified at
    all. And p can land far outside the range compiled from real sequences,
    which usually means the sample is too small or too sparse to constrain it.
    Neither suppresses the fit; both mark it.

    References
    ----------
    Utsu, T., Ogata, Y. and Matsu'ura, R. S. (1995).
    """
    cautions = []
    if fit.c < constants.OMORI_C_FLOOR:
        cautions.append(
            f"c has collapsed to {fit.c:.2g}, at the boundary of the model, so it "
            f"is not identified by this catalogue"
        )
    if not constants.OMORI_P_MIN <= fit.p <= constants.OMORI_P_MAX:
        cautions.append(
            f"p = {fit.p:.3f} lies outside {constants.OMORI_P_MIN} to "
            f"{constants.OMORI_P_MAX}, the range reported for real sequences"
        )
    return "; ".join(cautions) if cautions else None
