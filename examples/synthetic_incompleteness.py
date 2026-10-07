"""A catalogue with a known b and a known loss of detection, put through the package.

    python examples/synthetic_incompleteness.py            # about half an hour
    python examples/synthetic_incompleteness.py --quick    # a few minutes, for a check
    python examples/synthetic_incompleteness.py --workers 4

On a real catalogue the true b is not known, so a b that depends on the threshold can
be called unstable but not wrong. Here the truth is set. Magnitudes are drawn from a
Gutenberg-Richter law with b = 1, and each event is recorded with a probability that
rises smoothly with its magnitude, the cumulative normal of Ogata and Katsura (1993):
half the events at the magnitude mu are recorded, and the rise from 16 to 84 per cent
takes two sigma. That is a catalogue incomplete in the way real ones are, gradually,
with no magnitude below which nothing is seen and above which everything is.

Each simulated catalogue goes through the steps a user would take: the three
completeness estimates, b at each of them, and the plateau test. What is measured is how
far b at the estimated completeness magnitude lies from the true value, how often the
interval the package quotes contains it, and where the plateau test places the onset.
The first scenario has no incompleteness at all and shows what the estimators do when
nothing is wrong. The next four vary the width of the detection function around the
Hector Mine value, and the last repeats the Hector Mine value with magnitudes rounded to
0.01, as ComCat reports them, instead of to the 0.1 grid the estimator assumes. The
complete catalogues are also used to measure how often the approximate chi-squared
reference of the plateau test, the alternative to its simulated null, rejects a b that
is constant.

The detection function that goes with Hector Mine, mu 1.68 and sigma 0.46 at b = 1.0,
is not chosen by hand. `fit_detection` finds the b, mu and sigma that best reproduce
the observed curve of b against threshold, and the same fit is run on the other two
catalogues. It is a descriptive fit: the points on the curve share events, so the
residuals are correlated and their sum of squares is not a chi-squared statistic. What
it shows is whether a smooth detection loss with believable parameters can account
for the climb at all.

Every replicate has its own seed, derived from its scenario and index, so the result
does not depend on the number of workers and an interrupted run can be resumed with the
same command: finished scenarios are read back from the output file.

Writes docs/tremor_lab_synthetic_incompleteness.json.
"""

import argparse
import json
import math
import platform
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.stats import norm

from tremor_lab import (
    __version__,
    b_plateau,
    b_shift,
    b_stability,
    b_value_aki,
    mc_maxcurvature,
)
from tremor_lab.bvalue import mc_b_stability
from tremor_lab.completeness import mc_goodness_of_fit

HERE = Path(__file__).resolve().parent
OUT_PATH = HERE.parent / "docs" / "tremor_lab_synthetic_incompleteness.json"

DM = 0.1
B_TRUE = 1.0
LN10 = math.log(10.0)
Z95 = 1.959963984540054
# The thresholds at which the mean b curve is kept, and the plateau test's null draws.
CURVE = np.round(np.arange(1.0, 3.5001, DM), 2)
SIMULATIONS = 199
CHI_SQUARED_REPLICATES = 2000

# Hector Mine has 4,344 events at M 2.1 and above. At sigma 0.46, 170,000 events drawn
# from M 0.5 give about 4,390 there. They give about 17,200 recorded events from M 1.0,
# more than the 13,524 of the real file, so the simulated detection function is not a
# fit of the real network below M 1.7.
GRADUAL = {
    "b": B_TRUE,
    "mu": 1.68,
    "n_generated": 170_000,
    "floor": 0.5,
    "record_from": 1.0,
    "decimals": 1,
    "replicates": 500,
}
SCENARIOS = [
    {
        "name": "complete",
        "b": B_TRUE,
        "mu": None,
        "sigma": 0.0,
        "n_generated": 8_500,
        "floor": 1.65,
        "record_from": 1.7,
        "decimals": 1,
        "replicates": 500,
    },
    *(
        {**GRADUAL, "name": f"sigma {s:.2f}", "sigma": s}
        for s in (0.15, 0.30, 0.46, 0.70)
    ),
    {**GRADUAL, "name": "sigma 0.46, rounded to 0.01", "sigma": 0.46, "decimals": 2},
]


def wilson(successes: int, trials: int) -> tuple[float, float]:
    """Wilson score interval for a proportion, 95 per cent."""
    if trials == 0:
        return (float("nan"), float("nan"))
    p = successes / trials
    z2 = Z95**2
    denominator = 1.0 + z2 / trials
    centre = (p + z2 / (2.0 * trials)) / denominator
    half = (
        Z95 * math.sqrt(p * (1.0 - p) / trials + z2 / (4.0 * trials**2)) / denominator
    )
    return (max(0.0, centre - half), min(1.0, centre + half))


def catalogue(scenario: dict, rng: np.random.Generator) -> np.ndarray:
    """Recorded magnitudes: Gutenberg-Richter, thinned by the detection function."""
    m = scenario["floor"] + rng.exponential(
        1.0 / (scenario["b"] * LN10), scenario["n_generated"]
    )
    if scenario["mu"] is not None:
        seen = rng.random(m.size) < norm.cdf((m - scenario["mu"]) / scenario["sigma"])
        m = m[seen]
    m = np.round(m, scenario["decimals"])
    return m[m >= scenario["record_from"] - 1e-9]


def oracle_threshold(scenario: dict) -> float:
    """The lowest threshold at which 99 per cent of events are recorded.

    It is what an analyst who knew the detection function would choose, and gives the
    benchmark the estimated completeness magnitudes are compared against.
    """
    if scenario["mu"] is None:
        return scenario["record_from"]
    level = scenario["mu"] + 2.326 * scenario["sigma"]
    return float(np.round(np.ceil(level / DM - 1e-9) * DM, 2))


def b_at(mags: np.ndarray, threshold: float) -> dict | None:
    try:
        estimate = b_value_aki(mags, threshold, dm=DM)
    except ValueError:
        return None
    return {
        "mc": round(float(threshold), 2),
        "b": estimate.b,
        "sigma": estimate.sigma,
        "n": estimate.n,
    }


def one_replicate(job: tuple[int, int, int]) -> dict:
    """Simulate one catalogue and run the package's estimators on it."""
    scenario_index, index, simulations = job
    scenario = SCENARIOS[scenario_index]
    seed = 2_000_003 * (scenario_index + 1) + index
    rng = np.random.default_rng(seed)
    mags = catalogue(scenario, rng)

    row: dict = {"n": int(mags.size)}
    estimates = {
        "maxcurvature": mc_maxcurvature(mags, dm=DM),
        "goodness_of_fit": mc_goodness_of_fit(mags, dm=DM).mc,
        "b_stability": mc_b_stability(mags, dm=DM),
    }
    for label, mc in estimates.items():
        row[label] = None if mc is None else b_at(mags, mc)
    row["oracle"] = b_at(mags, oracle_threshold(scenario))

    plateau = b_plateau(mags, dm=DM, n_simulations=simulations, seed=seed)
    tested = np.flatnonzero(np.isfinite(plateau.p_profile))
    row["first_candidate"] = (
        None
        if tested.size == 0
        else {
            "threshold": round(float(plateau.thresholds[tested[0]]), 2),
            "p": float(plateau.p_profile[tested[0]]),
        }
    )
    row["plateau"] = (
        None
        if plateau.onset is None
        else {
            "mc": round(plateau.onset, 2),
            "b": plateau.b,
            "sigma": plateau.sigma,
            "n": plateau.n,
            "p": plateau.p_value,
        }
    )
    row["shift"] = None
    maxcurvature = estimates["maxcurvature"]
    if plateau.onset is not None and maxcurvature < plateau.onset:
        shift = b_shift(mags, maxcurvature, plateau.onset, dm=DM, n_boot=0)
        row["shift"] = {"shift": shift.shift, "z": shift.z_constant_b}

    curve = b_stability(mags, thresholds=CURVE, dm=DM)
    values: list[float | None] = [None] * CURVE.size
    for t, b in zip(curve.thresholds, curve.b, strict=True):
        values[round((t - CURVE[0]) / DM)] = float(b)
    row["curve"] = values
    return row


def expected_b(thresholds, b: float, mu: float, sigma: float) -> np.ndarray:
    """The b the Aki estimator converges to at each threshold under the detection loss.

    Above a floor f the recorded magnitudes have density proportional to
    10^(-b x) times the cumulative normal at x, so the estimator's mean excess has a
    closed form up to one integral, done here on a fine grid. Magnitudes are treated as
    continuous; rounding to the dm grid shifts the result by well under 0.01.
    """
    out = []
    for t in np.atleast_1d(thresholds):
        floor = t - DM / 2
        x = np.linspace(floor, 8.0, 20001)
        w = 10.0 ** (-b * x) * norm.cdf((x - mu) / sigma)
        mean = np.trapezoid(x * w, x) / np.trapezoid(w, x)
        out.append(1.0 / (LN10 * (mean - floor)))
    return np.array(out)


def fit_detection(thresholds, b, sigma) -> dict:
    """Best b, mu and detection width for an observed curve of b against threshold."""
    thresholds = np.asarray(thresholds, float)
    b = np.asarray(b, float)
    sigma = np.asarray(sigma, float)

    def residual(p):
        return (expected_b(thresholds, p[0], p[1], p[2]) - b) / sigma

    best = None
    lower = [0.3, thresholds[0] - 2.0, 0.02]
    upper = [2.0, thresholds[-1] + 1.0, 3.0]
    for mu0 in (thresholds[0] - 0.5, thresholds[0], thresholds[0] + 0.5):
        for width in (0.2, 0.5, 0.8):
            try:
                fit = least_squares(residual, [1.0, mu0, width], bounds=(lower, upper))
            except ValueError:
                continue
            if best is None or fit.cost < best.cost:
                best = fit
    b_true, mu, width = best.x
    return {
        "b_true": float(b_true),
        "mu": float(mu),
        "sigma": float(width),
        "weighted_ss": float(2.0 * best.cost),
        "points": int(thresholds.size),
        "at_bound": bool(
            np.any(np.isclose(best.x, lower, atol=1e-6))
            or np.any(np.isclose(best.x, upper, atol=1e-6))
        ),
        "thresholds": [round(float(t), 2) for t in thresholds],
        "observed": [round(float(v), 4) for v in b],
        "fitted": [
            round(float(v), 4) for v in expected_b(thresholds, b_true, mu, width)
        ],
    }


def chi_squared_size(replicates: int = CHI_SQUARED_REPLICATES) -> dict:
    """How often the chi-squared reference rejects a catalogue whose b is constant.

    The complete scenario again, with the seeds of its first replicates, so that the
    catalogues are the ones the simulated null was tried on, but with the plateau
    statistic referred to the chi-squared distribution in place of the simulated
    null. The simulated null's own rate is `plateau_first_candidate` of that
    scenario. The first candidate is the lowest threshold with at least three steps
    above it, and every rejection there is a false one.
    """
    scenario = SCENARIOS[0]
    rejected = 0
    used = 0
    for index in range(replicates):
        rng = np.random.default_rng(2_000_003 + index)
        plateau = b_plateau(catalogue(scenario, rng), dm=DM, n_simulations=0)
        tested = np.flatnonzero(np.isfinite(plateau.p_profile))
        if tested.size == 0:
            continue
        used += 1
        rejected += int(plateau.p_profile[tested[0]] < 0.05)
    lo, hi = wilson(rejected, used)
    return {
        "scenario": scenario["name"],
        "replicates": used,
        "rejected_at_0.05": rejected / used if used else float("nan"),
        "rejected_ci95": [lo, hi],
    }


def real_catalogues(min_events: int = 300) -> dict:
    """The detection model fitted to the three catalogues of the paper.

    The curve for each is taken over the thresholds that keep at least `min_events`
    events, so that the highest points are not a few dozen events each, and starts
    where the paper's own tables do: at the mode of the incremental distribution for
    the two USGS catalogues, at the file's floor of M 3.0 for Kahramanmaras.
    """
    import regenerate_paper_figures as paper

    sources = {
        "kahramanmaras": (
            paper.load_kahramanmaras()[0],
            np.round(np.arange(3.0, 4.9001, DM), 2),
        ),
        "hectormine": (paper.load_usgs("hectormine")[0], None),
        "ridgecrest": (paper.load_usgs("ridgecrest")[0], None),
    }
    out = {}
    for name, (mags, grid) in sources.items():
        curve = b_stability(mags, thresholds=grid, dm=DM, min_events=min_events)
        fit = fit_detection(curve.thresholds, curve.b, curve.sigma)
        fit["min_events"] = min_events
        out[name] = fit
    return out


def block(found: list[dict], replicates: int, b_true: float) -> dict:
    """Summary of one estimator over the replicates in which it returned a value."""
    result: dict = {"estimated": len(found) / replicates}
    if not found:
        return result
    mc = np.array([f["mc"] for f in found])
    b = np.array([f["b"] for f in found])
    sigma = np.array([f["sigma"] for f in found])
    inside = int(np.sum(np.abs(b - b_true) <= Z95 * sigma))
    lo, hi = wilson(inside, len(found))
    result.update(
        {
            "mc_median": float(np.median(mc)),
            "mc_p10": float(np.percentile(mc, 10)),
            "mc_p90": float(np.percentile(mc, 90)),
            "n_median": float(np.median([f["n"] for f in found])),
            "b_mean": float(b.mean()),
            "bias": float((b - b_true).mean()),
            "bias_se": float(b.std(ddof=1) / math.sqrt(b.size)) if b.size > 1 else None,
            "sigma_mean": float(sigma.mean()),
            "sd_of_b": float(b.std(ddof=1)) if b.size > 1 else None,
            "coverage": inside / len(found),
            "coverage_ci95": [lo, hi],
        }
    )
    return result


def summarise(scenario: dict, rows: list[dict], seconds: float) -> dict:
    n = len(rows)
    out = {
        "scenario": {k: v for k, v in scenario.items() if k != "replicates"},
        "replicates": n,
        "seconds": round(seconds, 1),
        "events_recorded_mean": float(np.mean([r["n"] for r in rows])),
        "oracle_threshold": oracle_threshold(scenario),
    }
    for label in (
        "maxcurvature",
        "goodness_of_fit",
        "b_stability",
        "plateau",
        "oracle",
    ):
        out[label] = block(
            [r[label] for r in rows if r[label] is not None], n, scenario["b"]
        )

    first = [r["first_candidate"] for r in rows if r["first_candidate"] is not None]
    rejected = sum(1 for f in first if f["p"] < 0.05)
    lo, hi = wilson(rejected, len(first))
    out["plateau_first_candidate"] = {
        "threshold_median": float(np.median([f["threshold"] for f in first])),
        "rejected_at_0.05": rejected / len(first),
        "rejected_ci95": [lo, hi],
    }

    # Only where the onset lies above the maximum-curvature magnitude is there a shift
    # to measure. Its significance is deliberately not summarised: the onset was
    # chosen by a test of those very differences, so a test of the shift between
    # that onset and the starting point would be circular.
    shifts = [r["shift"] for r in rows if r["shift"] is not None]
    lo, hi = wilson(len(shifts), n)
    out["shift_maxcurvature_to_onset"] = {
        "replicates_with_a_shift": len(shifts) / n,
        "replicates_ci95": [lo, hi],
        "median_shift": float(np.median([s["shift"] for s in shifts]))
        if shifts
        else None,
    }

    table = np.array([[np.nan if v is None else v for v in r["curve"]] for r in rows])
    with np.errstate(invalid="ignore"):
        mean_curve = np.nanmean(table, axis=0)
    counts = np.sum(np.isfinite(table), axis=0)
    expected = (
        None
        if scenario["mu"] is None
        else expected_b(CURVE, scenario["b"], scenario["mu"], scenario["sigma"])
    )
    out["curve"] = {
        "thresholds": [float(t) for t in CURVE],
        "mean_b": [
            None if not np.isfinite(v) else round(float(v), 4) for v in mean_curve
        ],
        "replicates": [int(c) for c in counts],
        "expected_b": None
        if expected is None
        else [round(float(v), 4) for v in expected],
    }
    return out


def show(result: dict) -> None:
    s = result["scenario"]
    where = "complete" if s["mu"] is None else f"mu {s['mu']}, sigma {s['sigma']}"
    print(
        f"\n{s['name']}: {where}, b {s['b']}, recorded to {s['decimals']} decimal(s), "
        f"{result['events_recorded_mean']:.0f} events, {result['replicates']} "
        f"replicates, {result['seconds']:.0f} s"
    )
    print(f"  threshold with 99 per cent detection: M {result['oracle_threshold']}")
    for label in (
        "maxcurvature",
        "goodness_of_fit",
        "b_stability",
        "plateau",
        "oracle",
    ):
        r = result[label]
        if "mc_median" not in r:
            print(f"  {label:16s} never estimated")
            continue
        lo, hi = r["coverage_ci95"]
        print(
            f"  {label:16s} estimated {r['estimated']:.0%}  Mc {r['mc_median']:.1f} "
            f"[{r['mc_p10']:.1f}, {r['mc_p90']:.1f}]  b {r['b_mean']:.3f} "
            f"(bias {r['bias']:+.3f})  95% interval holds b: {r['coverage']:.1%} "
            f"[{lo:.1%}, {hi:.1%}]"
        )
    f = result["plateau_first_candidate"]
    print(
        f"  plateau test rejects at its first candidate, M "
        f"{f['threshold_median']:.1f}, in {f['rejected_at_0.05']:.1%} "
        f"[{f['rejected_ci95'][0]:.1%}, {f['rejected_ci95'][1]:.1%}]"
    )
    sh = result["shift_maxcurvature_to_onset"]
    if sh["median_shift"] is not None:
        print(
            f"  the plateau onset lies above the maximum-curvature Mc in "
            f"{sh['replicates_with_a_shift']:.1%} of replicates, where b rises by a "
            f"median {sh['median_shift']:+.3f} between them"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--quick", action="store_true", help="8 replicates, 99 draws")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--out", type=Path, default=OUT_PATH)
    parser.add_argument(
        "--skip-real",
        action="store_true",
        help="leave out the fits to the real catalogues",
    )
    args = parser.parse_args()

    simulations = 99 if args.quick else SIMULATIONS
    previous: dict = {}
    if args.out.exists() and not args.quick:
        saved = json.loads(args.out.read_text(encoding="utf-8"))
        previous = {r["scenario"]["name"]: r for r in saved["scenarios"]}
        previous["__real__"] = saved.get("real_catalogues")
        previous["__reference__"] = saved.get("chi_squared_reference")

    results = []
    with Pool(args.workers) as pool:
        for si, scenario in enumerate(SCENARIOS):
            replicates = 8 if args.quick else scenario["replicates"]
            cached = previous.get(scenario["name"])
            if cached and cached["replicates"] == replicates:
                results.append(cached)
                show(cached)
                continue
            started = time.time()
            rows = pool.map(
                one_replicate,
                [(si, i, simulations) for i in range(replicates)],
                chunksize=2,
            )
            result = summarise(scenario, rows, time.time() - started)
            results.append(result)
            show(result)
            if not args.quick:
                write(
                    args.out,
                    results,
                    previous.get("__real__"),
                    simulations,
                    previous.get("__reference__"),
                )

    reference = previous.get("__reference__")
    if reference is None and not args.quick:
        reference = chi_squared_size()
        print(
            f"\nchi-squared reference on a constant b: rejects at the first candidate "
            f"in {reference['rejected_at_0.05']:.1%} of {reference['replicates']} "
            f"catalogues [{reference['rejected_ci95'][0]:.1%}, "
            f"{reference['rejected_ci95'][1]:.1%}]"
        )
    real = previous.get("__real__")
    if real is None and not args.skip_real and not args.quick:
        real = real_catalogues()
        for name, fit in real.items():
            print(
                f"\n{name}: b {fit['b_true']:.3f}, mu {fit['mu']:.2f}, sigma "
                f"{fit['sigma']:.2f} over M {fit['thresholds'][0]} to "
                f"{fit['thresholds'][-1]} ({fit['points']} points, weighted SS "
                f"{fit['weighted_ss']:.1f})"
            )
    if not args.quick:
        write(args.out, results, real, simulations, reference)
        print(f"\nwrote {args.out}")


def write(
    path: Path,
    results: list[dict],
    real: dict | None,
    simulations: int,
    reference: dict | None = None,
) -> None:
    payload = {
        "meta": {
            "tremor_lab_version": __version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "python": platform.python_version(),
            "dm": DM,
            "b_true": B_TRUE,
            "plateau_simulations": simulations,
        },
        "scenarios": results,
        "chi_squared_reference": reference,
        "real_catalogues": real,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
