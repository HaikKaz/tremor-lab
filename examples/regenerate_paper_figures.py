"""Regenerate every number in Section 5 of the paper, and the supplementary tables.

    python examples/regenerate_paper_figures.py

Reads the bundled 2023 Kahramanmaras catalogue (tests/data) and the 1999 Hector Mine
and 2019 Ridgecrest catalogues in this directory, recomputes each value from the
package's own public functions, and writes docs/tremor_lab_regeneration_data.json.
`examples/make_paper_figures.py` draws Figures 1 to 3 from that file.

What the file holds, and where each part appears in the paper:

    b_sweep       b, its standard error and n at every threshold (Table S1,
                  Figure 2)
    b_at_mc       b at each of the three completeness estimates (Table 4)
    omori_sweep   the decay refitted at each threshold (Figure 3, Section 5.4)
    time_bands    completeness magnitude in time bands (Table S5, Figure S1)
    start_cut     b after removing the first days of a sequence (Tables S6 and S7)
    band_b        b estimated in the time bands (Table S8, Supplement S7)
    fit_test      the calibrated residual test of a single decay (Table 5)
    refit_tests   the same test repeated at higher thresholds (Section 5.4)
    shifts        the change in b between thresholds, with a bootstrap interval
                  (Table 4, Supplement Table S4)
    plateau       the threshold above which b is statistically constant, and the
                  test's whole profile (Table 4, Supplement Table S2, Figure 2)
    plateau_grid_top
                  the plateau onset for grids that end at different magnitudes
                  (Supplement Table S3)
    floor         the lowest bins of the Kahramanmaras catalogue (Section 5.2)
    magnitude_types
                  the magnitude types behind the two USGS catalogues, and b with
                  the type Mh left out (Section 5.3, Supplement S5)

Run it from the project's own environment. It takes about five minutes, most of it the
calibrated decay-fit tests, which refit the decay once per replicate. Pass
--skip-fit-test to leave those out and keep the rest.

The synthetic experiment and the detection-model fits are in
`synthetic_incompleteness.py`, and the calibration of the decay-fit test in
`calibration_study.py`; each writes its own file in docs/.

The Hector Mine and Ridgecrest catalogues are the USGS ComCat downloads described in
the header of hectormine.toml and ridgecrest.toml. Fetching them again later will
return slightly different files, because ComCat revises magnitudes and adds events.
"""

import argparse
import json
from itertools import pairwise
from pathlib import Path

import numpy as np
import pandas as pd

from tremor_lab import (
    __version__,
    b_plateau,
    b_shift,
    b_value_aki,
    constants,
    fit_omori,
    mc_maxcurvature,
)
from tremor_lab.bvalue import b_stability, mc_b_stability
from tremor_lab.catalog import read_catalog
from tremor_lab.completeness import fmd, mc_goodness_of_fit
from tremor_lab.grid import at_or_above
from tremor_lab.omori import omori_fit_test, omori_sample

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT_PATH = REPO / "docs" / "tremor_lab_regeneration_data.json"

DM = 0.1
WINDOW_DAYS = 180
SEED = 0

# Time bands of Table S5, in days after the mainshock, edges (lo, hi].
BANDS = [
    (0.0, 0.20),
    (0.20, 0.38),
    (0.38, 0.60),
    (0.60, 1.0),
    (1.0, 2.0),
    (2.0, 5.0),
    (30.0, 180.0),
]

# Tables S6 and S7: start times (days removed from the front of the sequence) and
# the thresholds at which b is estimated after each cut.
START_DAYS = [0, 1, 2, 7]
START_THRESHOLDS = {
    "kahramanmaras": [3.5, 3.8, 4.0, 4.4],
    "ridgecrest": [1.3, 1.5, 2.0, 2.5, 3.0],
}

# Table 5: the threshold each catalogue is tested at, and the replicates behind its
# p-value. 5,000 replicates for Catalogue A give a p-value with a standard error of
# about 0.003; 600 is the package default and is enough where p is far below 0.05.
FIT_TEST_THRESHOLD = {"kahramanmaras": 3.5, "hectormine": 1.7, "ridgecrest": 1.3}
FIT_TEST_REPLICATES = {"kahramanmaras": 5000, "hectormine": 600, "ridgecrest": 600}

# Section 5.4: the same test repeated at higher thresholds, where the decay is refitted
# above the early-time completeness magnitude. 2,000 replicates give a p-value with a
# standard error of about 0.005 near 0.05, enough to tell a verdict from a borderline.
REFIT_THRESHOLDS = {
    "kahramanmaras": [3.5, 3.8, 4.0, 4.4],
    "hectormine": [1.7, 2.0, 2.2, 2.6],
}
REFIT_REPLICATES = 2000

# The change in b between thresholds (Section 5.2): resamples behind each bootstrap
# interval, and draws behind the plateau test's simulated null. 5,000 resamples put the
# ends of a 95 per cent interval at about the 125th value from each side; 1,000 draws
# give a p-value near 0.05 a standard error of 0.007.
SHIFT_BOOTSTRAPS = 5000
PLATEAU_DRAWS = 1000
# Where the grid of candidate thresholds ends, for the check on how far the plateau
# onset depends on it (Supplement S3). Each catalogue's default top is the last entry.
PLATEAU_TOPS = {
    "kahramanmaras": (4.0, 4.2, 4.4, 4.6, 4.9),
    "hectormine": (3.0, 3.2, 3.4, 3.6, 4.0),
    "ridgecrest": (2.6, 2.8, 3.0, 3.3, 3.6),
}
# Contrasts every catalogue gets, in magnitude units above its maximum-curvature
# completeness magnitude. They are the same for all three, so no catalogue's contrast
# is chosen to suit it.
SHIFT_STEPS = (0.5, 1.0)
# Thresholds at which b is compared across subsets of the magnitude types (Section 5.3).
TYPE_THRESHOLDS = np.round(np.arange(1.5, 3.6001, DM), 2)
# Magnitude bands of the composition table, lower edges.
TYPE_BANDS = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5)

# The origin times are the mainshocks' own ComCat timestamps, fractions of a second
# included. The window excludes events at or before the origin time, so a time
# rounded down to the whole second leaves the mainshock's own record inside its
# aftershock sample as the first event. Both USGS catalogues list it a fraction of a
# second after the rounded time (0.46 s and 0.04 s), and the single M 7.1 event moves
# b at high thresholds by up to 0.1. `load_usgs` refuses a sample that contains it.
MAINSHOCKS = {
    "hectormine": {
        "t": "1999-10-16 09:46:44.46",
        "lat": 34.6033,
        "lon": -116.265,
        "mw": 7.1,
    },
    "ridgecrest": {
        "t": "2019-07-06 03:19:53.04",
        "lat": 35.7695,
        "lon": -117.5993,
        "mw": 7.1,
    },
}

USGS_COLUMNS = {
    "datetime": "time",
    "lat": "latitude",
    "lon": "longitude",
    "mag": "mag",
    "mag_type": "magType",
}


def b_row(mags, threshold):
    """b, its standard error and n at one threshold, rounded as the paper quotes."""
    b = b_value_aki(mags, threshold, dm=DM)
    return {
        "threshold": round(float(threshold), 2),
        "b": round(b.b, 3),
        "sigma": round(b.sigma, 3),
        "n": int(b.n),
    }


def mc_estimates(mags):
    return {
        "maxcurvature": mc_maxcurvature(mags, dm=DM),
        "goodness_of_fit": mc_goodness_of_fit(mags, dm=DM).mc,
        "b_stability": mc_b_stability(mags, dm=DM),
    }


def b_at_each_mc(mags, estimates):
    return {
        name: None if mc is None else b_row(mags, mc) for name, mc in estimates.items()
    }


def b_sweep(mags, thresholds=None):
    curve = b_stability(mags, thresholds=thresholds, dm=DM)
    return [
        {
            "threshold": round(float(t), 2),
            "b": round(float(b), 3),
            "sigma": round(float(s), 3),
            "n": int(n),
        }
        for t, b, s, n in zip(
            curve.thresholds, curve.b, curve.sigma, curve.n, strict=True
        )
    ]


def omori_sweep(mags, times, thresholds):
    out = []
    for thr in thresholds:
        sel = at_or_above(mags, thr - DM / 2)
        t = times[sel]
        t = t[t > 0]
        if t.size < 2:
            continue
        fit = fit_omori(t)
        out.append(
            {
                "threshold": round(float(thr), 2),
                "p": round(fit.p, 4),
                "c": round(fit.c, 4),
                "k": round(fit.k, 2),
                "n": int(fit.n),
            }
        )
    return out


def time_bands(mags, times):
    rows = []
    for lo, hi in BANDS:
        band = mags[(times > lo) & (times <= hi)]
        n = int(band.size)
        mc = mc_maxcurvature(band, dm=DM) if n >= 10 else None
        rows.append({"band_lo": lo, "band_hi": hi, "n": n, "mc": mc})
    return rows


def band_b(mags, times, fixed_threshold=None):
    """b in each time band, at the band's own Mc or at one fixed threshold.

    Each row carries the threshold it was estimated at, which is the band's own
    maximum-curvature completeness magnitude unless one was fixed.
    """
    rows = []
    for lo, hi in BANDS:
        band = mags[(times > lo) & (times <= hi)]
        if band.size < 10:
            continue
        if fixed_threshold is None:
            threshold = mc_maxcurvature(band, dm=DM)
        else:
            threshold = fixed_threshold
        try:
            row = b_row(band, threshold)
        except ValueError:
            continue
        rows.append({"band_lo": lo, "band_hi": hi, **row})
    return rows


def start_cut(mags, times, thresholds):
    rows = []
    for start in START_DAYS:
        keep = times > start
        cells = [b_row(mags[keep], thr) for thr in thresholds]
        rows.append({"start_days": start, "cells": cells})
    return rows


def null_statistics(times, fit, n_simulations, seed):
    """Simulated Kolmogorov-Smirnov statistics under the fitted decay.

    The same loop as `tremor_lab.omori.omori_fit_test` (draw from the fitted model,
    refit, take the statistic), kept here so the 95th percentile of the null, which
    the package does not return, can be reported beside the p-value.
    """
    t_end = float(np.max(times))
    rng = np.random.default_rng(seed)
    stats = []
    for _ in range(n_simulations):
        drawn = np.sort(omori_sample(times.size, fit.p, fit.c, t_end, 0.0, rng))
        try:
            refit = fit_omori(drawn, t_end=t_end, t_start=0.0)
        except ValueError:
            continue
        stat = omori_fit_test(
            drawn, refit, t_end=t_end, t_start=0.0, n_simulations=0
        ).statistic
        if np.isfinite(stat):
            stats.append(stat)
    return np.array(stats)


def fit_test(mags, times, threshold, n_simulations):
    """Table 5 row: the residual test of one decay, with its simulated null."""
    sel = at_or_above(mags, threshold - DM / 2) & (times > 0)
    t = np.sort(times[sel])
    fit = fit_omori(t)
    package = omori_fit_test(t, fit, n_simulations=n_simulations, seed=SEED)
    null = null_statistics(t, fit, n_simulations, SEED)
    p_value = (1.0 + np.sum(null >= package.statistic)) / (1.0 + null.size)
    if abs(p_value - package.p_value) > 1e-12:
        raise RuntimeError(
            f"the simulated null no longer matches the package's own test: "
            f"{p_value} against {package.p_value}"
        )
    return {
        "threshold": threshold,
        "n": int(package.n),
        "ks": round(package.statistic, 4),
        "null_95th_percentile": round(float(np.quantile(null, 0.95)), 4),
        "p": round(package.p_value, 4),
        "p_se": round(package.p_value_se, 4),
        "replicates": int(null.size),
        "method": package.method,
    }


def refit_tests(mags, times, thresholds):
    """The residual test of one decay, repeated at each of several thresholds."""
    rows = []
    for threshold in thresholds:
        sel = at_or_above(mags, threshold - DM / 2) & (times > 0)
        t = np.sort(times[sel])
        result = omori_fit_test(
            t, fit_omori(t), n_simulations=REFIT_REPLICATES, seed=SEED
        )
        rows.append(
            {
                "threshold": threshold,
                "n": int(result.n),
                "ks": round(result.statistic, 4),
                "p": round(result.p_value, 4),
                "p_se": round(result.p_value_se, 4),
                "replicates": REFIT_REPLICATES,
            }
        )
    return rows


def shift_row(mags, contrast, mc_low, mc_high):
    """b at two nested thresholds and what the difference is worth.

    Two errors are kept because they answer different questions: the closed-form
    error holds if b is the same at both thresholds and tests that, the bootstrap
    holds either way and sizes the shift. `ratio_to_low_error` is the shift over the
    Shi and Bolt error at the lower threshold, the comparison the revised paper no
    longer makes; it is here so the difference between the two readings stays visible.
    """
    shift = b_shift(mags, mc_low, mc_high, dm=DM, n_boot=SHIFT_BOOTSTRAPS, seed=SEED)
    low_error = b_value_aki(mags, mc_low, dm=DM).sigma
    lo, hi = shift.ci_boot
    return {
        "contrast": contrast,
        "mc_low": round(shift.mc_low, 2),
        "mc_high": round(shift.mc_high, 2),
        "n_low": shift.n_low,
        "n_high": shift.n_high,
        "b_low": round(shift.b_low, 3),
        "b_high": round(shift.b_high, 3),
        "shift": round(shift.shift, 3),
        "ci95": [round(lo, 3), round(hi, 3)],
        "se_boot": round(shift.se_boot, 4),
        "z_boot": round(shift.shift / shift.se_boot, 1),
        "se_constant_b": round(shift.se_constant_b, 4),
        "z_constant_b": round(shift.z_constant_b, 1),
        "ratio_to_low_error": round(shift.shift / low_error, 1),
        "resamples": SHIFT_BOOTSTRAPS,
    }


def plateau_block(mags, thresholds):
    """The plateau test, with the simulated null, and its chi-squared counterpart."""
    plateau = b_plateau(
        mags,
        thresholds=thresholds,
        dm=DM,
        n_simulations=PLATEAU_DRAWS,
        seed=SEED,
    )
    approximate = b_plateau(mags, thresholds=thresholds, dm=DM, n_simulations=0)

    def cell(value, digits):
        return None if not np.isfinite(value) else round(float(value), digits)

    profile = [
        {
            "threshold": round(float(t), 2),
            "statistic": cell(stat, 2),
            "steps": int(steps),
            "p": cell(p, 4),
        }
        for t, stat, steps, p in zip(
            plateau.thresholds,
            plateau.statistic_profile,
            plateau.dof_profile,
            plateau.p_profile,
            strict=True,
        )
    ]
    return {
        "onset": plateau.onset if plateau.onset is None else round(plateau.onset, 2),
        "b": None if plateau.b is None else round(plateau.b, 3),
        "sigma": None if plateau.sigma is None else round(plateau.sigma, 3),
        "n": plateau.n,
        "statistic": None if plateau.statistic is None else round(plateau.statistic, 2),
        "steps": plateau.dof,
        "p": None if plateau.p_value is None else round(plateau.p_value, 4),
        "alpha": plateau.alpha,
        "method": plateau.method,
        "onset_chi_squared": approximate.onset,
        "profile": profile,
    }


def plateau_by_grid_top(mags, tops, first_candidate):
    """The plateau onset for grids that start where the default one does and end at
    each of `tops`.

    The onset is the lowest candidate that passes with at least three thresholds
    above it, so it depends on where the grid ends: a short grid has fewer steps to
    reject on, and a long one reaches thin samples whose errors are wide. The table
    shows how far the answer moves, so that the onset is read as a property of the
    grid as well as of the catalogue.
    """
    rows = []
    for top in tops:
        grid = np.round(np.arange(first_candidate, top + 1e-9, DM), 2)
        plateau = b_plateau(
            mags, thresholds=grid, dm=DM, n_simulations=PLATEAU_DRAWS, seed=SEED
        )
        rows.append(
            {
                "top": round(float(top), 2),
                "onset": None if plateau.onset is None else round(plateau.onset, 2),
                "b": None if plateau.b is None else round(plateau.b, 3),
                "n": plateau.n,
            }
        )
    return rows


def shifts_for(mags, estimates, plateau):
    """The contrasts of Table 4, from the maximum-curvature magnitude upward."""
    low = estimates["maxcurvature"]
    contrasts = [(f"{step:+.1f} above Mc", low + step) for step in SHIFT_STEPS]
    if estimates["b_stability"] is not None:
        contrasts.append(("b-stability", estimates["b_stability"]))
    if plateau["onset"] is not None:
        contrasts.append(("plateau onset", plateau["onset"]))
    rows, seen = [], set()
    for label, high in contrasts:
        high = round(float(high), 2)
        if high <= low or (round(low, 2), high) in seen:
            continue
        try:
            rows.append(shift_row(mags, label, low, high))
        except ValueError:
            continue
        seen.add((round(low, 2), high))
    return rows


def floor_table(mags, bins=8):
    """The lowest incremental bins, for a catalogue truncated at a floor."""
    edges, incremental, _ = fmd(mags, DM)
    start = int(np.argmax(incremental > 0))
    return {
        "floor": round(float(mags.min()), 2),
        "bins": [
            {"magnitude": round(float(e), 2), "n": int(n)}
            for e, n in zip(
                edges[start : start + bins],
                incremental[start : start + bins],
                strict=True,
            )
        ],
        "mode": round(float(edges[int(np.argmax(incremental))]), 2),
    }


def type_groups(labels):
    """Magnitude-type labels reduced to the groups the paper discusses."""
    clean = np.char.lower(np.char.strip(np.asarray(labels, dtype=str)))
    out = np.full(clean.shape, "other", dtype=object)
    for group in ("ml", "mlr", "mh", "mc", "md", "mb", "ms"):
        out[clean == group] = group
    out[np.char.startswith(clean, "mw")] = "mw"
    return out


def magnitude_types(mags, groups, scale_counts, subsets):
    """Which types make up each magnitude band, and b with some types left out.

    `scale_counts` is what `read_catalog` converted: only surface-wave and body-wave
    magnitudes are converted to Mw, and the local, duration, coda and hand-read types
    are taken as the network reports them. `subsets` maps a name to the types to keep,
    or to None for all of them. A subset that is one type is not a thinned copy of the
    complete catalogue: wherever that type is a minority its distribution is
    incomplete by construction, so a curve of b for it is a check on the plateau and
    not a second estimate of b.
    """
    present = sorted(set(groups))
    bands = []
    for lo, hi in pairwise([*TYPE_BANDS, np.inf]):
        sel = (mags >= lo - 1e-9) & (mags < hi - 1e-9)
        counts = {g: int(np.sum(groups[sel] == g)) for g in present}
        bands.append(
            {
                "band_lo": lo,
                "band_hi": None if np.isinf(hi) else hi,
                "n": int(sel.sum()),
                "counts": counts,
            }
        )
    sweeps = {}
    for name, kept_types in subsets.items():
        keep = (
            np.ones(mags.size, bool)
            if kept_types is None
            else np.isin(groups, kept_types)
        )
        sub = mags[keep]
        curve = b_stability(sub, thresholds=TYPE_THRESHOLDS, dm=DM)
        rows = [
            {"threshold": round(float(t), 2), "b": round(float(b), 3), "n": int(n)}
            for t, b, n in zip(curve.thresholds, curve.b, curve.n, strict=True)
        ]
        # The default grid, as for the whole catalogue, so that the all-types row of
        # this table is the plateau of Table 4 and not a second answer to it.
        onset = b_plateau(sub, dm=DM, n_simulations=PLATEAU_DRAWS, seed=SEED).onset
        sweeps[name] = {
            "n": int(sub.size),
            "b_sweep": rows,
            "plateau_onset": onset if onset is None else round(onset, 2),
        }
    return {
        "scale_counts_in_window": scale_counts,
        "types": present,
        "bands": bands,
        "subsets": sweeps,
    }


def event_type_check(name, mags, frame):
    """Rows that ComCat types as something other than an earthquake, and what they do.

    The two USGS files keep the rows ComCat labels `quarry blast`. They are a few dozen
    low-magnitude events inside the window. This leaves them out and compares the b
    curve, the maximum-curvature magnitude and the plateau onset with the full sample,
    so that the paper can say how much they matter instead of assuming it.
    """
    kinds = pd.read_csv(HERE / f"{name}.csv", usecols=["type"])["type"]
    kinds = kinds.to_numpy(dtype=str)[np.asarray(frame.attrs["kept_rows"])]
    other = kinds != "earthquake"
    kept = mags[~other]
    full = b_stability(mags, dm=DM)
    cut = b_stability(kept, dm=DM)
    b_full = {
        round(float(t), 2): float(b)
        for t, b in zip(full.thresholds, full.b, strict=True)
    }
    b_cut = {
        round(float(t), 2): float(b) for t, b in zip(cut.thresholds, cut.b, strict=True)
    }
    shared = sorted(set(b_full) & set(b_cut))
    change = {t: abs(b_full[t] - b_cut[t]) for t in shared}
    worst = max(change, key=change.get)

    def onset(sample):
        value = b_plateau(sample, dm=DM, n_simulations=PLATEAU_DRAWS, seed=SEED).onset
        return value if value is None else round(value, 2)

    return {
        "kinds": {str(k): int(np.sum(kinds == k)) for k in sorted(set(kinds[other]))},
        "n_other": int(other.sum()),
        "thresholds_compared": len(shared),
        "max_abs_change_in_b": round(change[worst], 4),
        "at_threshold": worst,
        "mc_maxcurvature": [
            mc_maxcurvature(mags, dm=DM),
            mc_maxcurvature(kept, dm=DM),
        ],
        "plateau_onset": [onset(mags), onset(kept)],
    }


def grid_reading_check(mags, sweep):
    """b from the magnitudes as published against b after rounding them to the 0.1 grid.

    ComCat gives most magnitudes to two decimals, and the half-bin floor of the
    estimator assumes they sit on the 0.1 grid. The paper uses the published values.
    This rounds the same events half up to the grid and repeats the sweep, so that the
    paper can say how far the two readings differ at each threshold, and whether the
    difference could produce the climb of b with threshold.
    """
    on_grid = np.floor(mags / DM + 0.5 + 1e-9) * DM
    rows = []
    for row in sweep:
        published = b_value_aki(mags, row["threshold"], dm=DM)
        rounded = b_value_aki(on_grid, row["threshold"], dm=DM)
        rows.append(
            {
                "threshold": row["threshold"],
                "b_published": round(float(published.b), 4),
                "b_grid": round(float(rounded.b), 4),
                "n_published": int(published.n),
                "n_grid": int(rounded.n),
            }
        )
    return rows


def analyse(name, mags, times, thresholds, omori_thresholds, skip_fit_test):
    estimates = mc_estimates(mags)
    plateau = plateau_block(mags, thresholds)
    result = {
        "n_events": int(mags.size),
        "mc": estimates,
        "b_at_mc": b_at_each_mc(mags, estimates),
        "b_sweep": b_sweep(mags, thresholds),
        "plateau": plateau,
        "plateau_grid_top": plateau_by_grid_top(
            mags, PLATEAU_TOPS[name], plateau["profile"][0]["threshold"]
        ),
        "shifts": shifts_for(mags, estimates, plateau),
    }
    if omori_thresholds is not None:
        result["omori_sweep"] = omori_sweep(mags, times, omori_thresholds)
    if name in START_THRESHOLDS:
        result["start_cut"] = start_cut(mags, times, START_THRESHOLDS[name])
        result["band_b"] = band_b(mags, times)
    if name == "kahramanmaras":
        result["time_bands"] = time_bands(mags, times)
        result["band_b_fixed_3p5"] = band_b(mags, times, fixed_threshold=3.5)
    if not skip_fit_test:
        result["fit_test"] = fit_test(
            mags, times, FIT_TEST_THRESHOLD[name], FIT_TEST_REPLICATES[name]
        )
        if name in REFIT_THRESHOLDS:
            result["refit_tests"] = refit_tests(mags, times, REFIT_THRESHOLDS[name])
    return result


def load_kahramanmaras():
    df = pd.read_csv(REPO / "tests" / "data" / "kahramanmaras_180d.csv")
    mags = df["mw"].to_numpy(float)
    times = df["dt_days"].to_numpy(float)
    keep = (times > 0) & (times <= WINDOW_DAYS)
    return mags[keep], times[keep]


def load_usgs_frame(name):
    """The analysed window as a table, with each event's magnitude-type label."""
    mainshock = MAINSHOCKS[name]
    df = read_catalog(
        HERE / f"{name}.csv",
        columns=USGS_COLUMNS,
        mainshock=mainshock,
        window_days=WINDOW_DAYS,
        radius_km=100,
    )
    # No aftershock is as large as the mainshock it follows. One that is, in a
    # sample that is meant to hold only aftershocks, is the mainshock's own record,
    # counted because the origin time given here was earlier than the catalogue's.
    if (df["mw"] >= mainshock["mw"]).any():
        raise ValueError(
            f"{name}: an event of M {mainshock['mw']} or larger is inside the "
            f"aftershock window. That is the mainshock itself; give its origin time "
            f"to the precision of the catalogue's own timestamp"
        )
    raw = pd.read_csv(HERE / f"{name}.csv", usecols=[USGS_COLUMNS["mag_type"]])
    labels = raw[USGS_COLUMNS["mag_type"]].to_numpy(dtype=str)[
        np.asarray(df.attrs["kept_rows"])
    ]
    return df, labels


def load_usgs(name):
    df, _ = load_usgs_frame(name)
    return df["mw"].to_numpy(float), df["dt_days"].to_numpy(float)


def summarise(results):
    for name, r in results.items():
        if name == "meta":
            continue
        print(f"\n{name}: n={r['n_events']}  Mc {r['mc']}")
        for label, row in r["b_at_mc"].items():
            if row is not None:
                print(f"  b at {label:16s} M {row['threshold']}: {row['b']}")
        pl = r["plateau"]
        print(
            f"  plateau onset M {pl['onset']} "
            f"(chi-squared reference: {pl['onset_chi_squared']})"
        )
        for row in r["shifts"]:
            print(
                f"  shift {row['contrast']:>14s} M {row['mc_low']} to "
                f"{row['mc_high']}: {row['shift']:+.3f} {row['ci95']} "
                f"z {row['z_boot']}"
            )
        if "fit_test" in r:
            ft = r["fit_test"]
            print(
                f"  fit test at M {ft['threshold']}: n={ft['n']} KS {ft['ks']} "
                f"null95 {ft['null_95th_percentile']} p {ft['p']} +/- {ft['p_se']}"
            )
        for row in r.get("refit_tests", []):
            print(
                f"  refit at M {row['threshold']}: n={row['n']} KS {row['ks']} "
                f"p {row['p']} +/- {row['p_se']}"
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--skip-fit-test",
        action="store_true",
        help="leave out the calibrated decay-fit test of Table 5 (the slow part)",
    )
    args = parser.parse_args()

    mags_a, times_a = load_kahramanmaras()
    frame_b, labels_b = load_usgs_frame("hectormine")
    frame_c, labels_c = load_usgs_frame("ridgecrest")
    mags_b, times_b = frame_b["mw"].to_numpy(float), frame_b["dt_days"].to_numpy(float)
    mags_c, times_c = frame_c["mw"].to_numpy(float), frame_c["dt_days"].to_numpy(float)

    # The Kahramanmaras catalogue is truncated at M 3.0, so its sweep is taken from
    # there. The default grid starts at the mode of the incremental distribution.
    grid_a = np.round(np.arange(3.0, 4.9001, DM), 2)

    results = {
        "meta": {
            "tremor_lab_version": __version__,
            "dm": DM,
            "window_days": WINDOW_DAYS,
            "seed": SEED,
            "n_fit_simulations_default": constants.N_FIT_SIMULATIONS,
            "shift_bootstraps": SHIFT_BOOTSTRAPS,
            "plateau_draws": PLATEAU_DRAWS,
        },
        "kahramanmaras": analyse(
            "kahramanmaras",
            mags_a,
            times_a,
            grid_a,
            np.round(np.arange(3.5, 4.4001, DM), 2),
            args.skip_fit_test,
        ),
        "hectormine": analyse(
            "hectormine",
            mags_b,
            times_b,
            None,
            np.round(np.arange(1.7, 2.6001, DM), 2),
            args.skip_fit_test,
        ),
        "ridgecrest": analyse(
            "ridgecrest", mags_c, times_c, None, None, args.skip_fit_test
        ),
    }

    results["kahramanmaras"]["floor"] = floor_table(mags_a)
    results["hectormine"]["event_types"] = event_type_check(
        "hectormine", mags_b, frame_b
    )
    results["ridgecrest"]["event_types"] = event_type_check(
        "ridgecrest", mags_c, frame_c
    )
    results["hectormine"]["grid_reading"] = grid_reading_check(
        mags_b, results["hectormine"]["b_sweep"]
    )
    results["ridgecrest"]["grid_reading"] = grid_reading_check(
        mags_c, results["ridgecrest"]["b_sweep"]
    )
    results["hectormine"]["magnitude_types"] = magnitude_types(
        mags_b,
        type_groups(labels_b),
        frame_b.attrs["scale_counts_in_window"],
        {
            "all types": None,
            "without mh": ["ml", "mc", "md", "mb", "mw"],
            "ml only": ["ml"],
        },
    )
    results["ridgecrest"]["magnitude_types"] = magnitude_types(
        mags_c,
        type_groups(labels_c),
        frame_c.attrs["scale_counts_in_window"],
        {"all types": None, "local types only (ml, mlr)": ["ml", "mlr"]},
    )

    OUT_PATH.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_PATH}")
    summarise(results)


if __name__ == "__main__":
    main()
