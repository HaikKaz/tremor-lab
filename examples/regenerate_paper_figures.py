"""Regenerate every number in Section 4 of the paper, and the supplementary tables.

    python examples/regenerate_paper_figures.py

Reads the bundled 2023 Kahramanmaras catalogue (tests/data) and the 1999 Hector Mine
and 2019 Ridgecrest catalogues in this directory, recomputes each value from the
package's own public functions, and writes docs/tremor_lab_regeneration_data.json.
`examples/make_paper_figures.py` draws Figures 1 to 3 from that file.

What the file holds, and where each part appears in the paper:

    b_sweep       b, its standard error and n at every threshold (Tables 2 to 4,
                  Figure 1)
    b_at_mc       b at each of the three completeness estimates (Section 4.2)
    omori_sweep   the decay refitted at each threshold (Figure 3, Section 4.4)
    time_bands    completeness magnitude in time bands (Table 5, Figure 2)
    start_cut     b after removing the first days of a sequence (Tables S1 and S2)
    band_b        b estimated in the time bands (Table S3, Section 4.5)
    fit_test      the calibrated residual test of a single decay (Table 6)
    refit_tests   the same test repeated at higher thresholds (Section 4.4)

Run it from the project's own environment. It takes about three minutes, nearly all of
it the calibrated decay-fit tests, which refit the decay once per replicate. Pass
--skip-fit-test to leave those out and keep the rest.

The Hector Mine and Ridgecrest catalogues are the USGS ComCat downloads described in
the header of hectormine.toml and ridgecrest.toml. Fetching them again later will
return slightly different files, because ComCat revises magnitudes and adds events.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tremor_lab import __version__, b_value_aki, constants, fit_omori, mc_maxcurvature
from tremor_lab.bvalue import b_stability, mc_b_stability
from tremor_lab.catalog import read_catalog
from tremor_lab.completeness import mc_goodness_of_fit
from tremor_lab.grid import at_or_above
from tremor_lab.omori import omori_fit_test, omori_sample

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT_PATH = REPO / "docs" / "tremor_lab_regeneration_data.json"

DM = 0.1
WINDOW_DAYS = 180
SEED = 0

# Time bands of Table 5, in days after the mainshock, edges (lo, hi].
BANDS = [
    (0.0, 0.20),
    (0.20, 0.38),
    (0.38, 0.60),
    (0.60, 1.0),
    (1.0, 2.0),
    (2.0, 5.0),
    (30.0, 180.0),
]

# Tables S1 and S2: start times (days removed from the front of the sequence) and
# the thresholds at which b is estimated after each cut.
START_DAYS = [0, 1, 2, 7]
START_THRESHOLDS = {
    "kahramanmaras": [3.5, 3.8, 4.0, 4.4],
    "ridgecrest": [1.3, 1.5, 2.0, 2.5, 3.0],
}

# Table 6: the threshold each catalogue is tested at, and the replicates behind its
# p-value. 5,000 replicates for Catalogue A give a p-value with a standard error of
# about 0.003; 600 is the package default and is enough where p is far below 0.05.
FIT_TEST_THRESHOLD = {"kahramanmaras": 3.5, "hectormine": 1.7, "ridgecrest": 1.3}
FIT_TEST_REPLICATES = {"kahramanmaras": 5000, "hectormine": 600, "ridgecrest": 600}

# Section 4.4: the same test repeated at higher thresholds, where the decay is refitted
# above the early-time completeness magnitude. 2,000 replicates give a p-value with a
# standard error of about 0.005 near 0.05, enough to tell a verdict from a borderline.
REFIT_THRESHOLDS = {
    "kahramanmaras": [3.5, 3.8, 4.0, 4.4],
    "hectormine": [1.7, 2.0, 2.2, 2.6],
}
REFIT_REPLICATES = 2000

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
    """Table 6 row: the residual test of one decay, with its simulated null."""
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


def analyse(name, mags, times, thresholds, omori_thresholds, skip_fit_test):
    estimates = mc_estimates(mags)
    result = {
        "n_events": int(mags.size),
        "mc": estimates,
        "b_at_mc": b_at_each_mc(mags, estimates),
        "b_sweep": b_sweep(mags, thresholds),
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


def load_usgs(name):
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
    return df["mw"].to_numpy(float), df["dt_days"].to_numpy(float)


def summarise(results):
    for name, r in results.items():
        if name == "meta":
            continue
        print(f"\n{name}: n={r['n_events']}  Mc {r['mc']}")
        for label, row in r["b_at_mc"].items():
            if row is not None:
                print(f"  b at {label:16s} M {row['threshold']}: {row['b']}")
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
        help="leave out the calibrated decay-fit test of Table 6 (the slow part)",
    )
    args = parser.parse_args()

    mags_a, times_a = load_kahramanmaras()
    mags_b, times_b = load_usgs("hectormine")
    mags_c, times_c = load_usgs("ridgecrest")

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

    OUT_PATH.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_PATH}")
    summarise(results)


if __name__ == "__main__":
    main()
