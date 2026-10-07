"""How often the decay-fit test rejects a decay that is correct, and one that is not.

    python examples/calibration_study.py                 # the full study, about an hour
    python examples/calibration_study.py --quick         # a few minutes, for a check
    python examples/calibration_study.py --workers 4

`omori_fit_test` asks whether a fitted modified Omori-Utsu decay describes the times it
was fitted to. The decay's two parameters are estimated from those same times, so the
textbook Kolmogorov-Smirnov p-value does not apply and the package replaces it with one
obtained by parametric bootstrap, refitting on every replicate. This script measures
what each of the two p-values does, on sequences whose answer is known.

For every configuration it simulates sequences from a modified Omori-Utsu process, so
the decay is correct by construction, fits each, and runs both versions of the test on
it exactly as `analyze_case` would: the asymptotic p-value (``n_simulations=0``) and
the calibrated one. The share of sequences rejected at a level is the realised size of
the test at that level; a test with the right size rejects 5 per cent of correct decays
at the 5 per cent level. Each share is reported with a Wilson 95 per cent interval,
because a rate measured on a few hundred sequences carries a sampling error of its own
and a figure quoted without one cannot be compared with its target.

The last configuration is not a null. Twenty per cent of its events are a second
Omori sequence beginning one day after the first, which is what a large aftershock does
to a window, and it shows what the test is for: the share it rejects is the power.

A calibrated p-value from B replicates has exactly the form (1 + k) / (1 + B), which is
uniform on a grid when the null is exactly right; B is 199 here so that the 5 per cent
level falls on the grid. The shipped default of 600 sharpens individual p-values and
does not change the size.

Every sequence has its own seed, derived from its configuration and index, so the
result does not depend on the number of workers and an interrupted run can be resumed
with the same command: finished configurations are read back from the output file.

Writes docs/tremor_lab_calibration_study.json.
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

from tremor_lab import __version__, fit_omori
from tremor_lab.omori import omori_fit_test, omori_sample

HERE = Path(__file__).resolve().parent
OUT_PATH = HERE.parent / "docs" / "tremor_lab_calibration_study.json"

WINDOW_DAYS = 180.0
LEVELS = (0.10, 0.05, 0.01)
REPLICATES = 199
Z95 = 1.959963984540054

# Each configuration is a family of sequences. The first is the decay fitted to the
# 2023 Kahramanmaras catalogue above M 3.5 (1,529 events, p 1.16, c 0.50 d); the
# others vary one thing at a time from it. `second` is the share of events that belong
# to a second Omori sequence beginning at `second_start` days, zero for a null.
CONFIGS = [
    {
        "name": "reference",
        "n": 1500,
        "p": 1.16,
        "c": 0.50,
        "second": 0.0,
        "sequences": 2000,
    },
    {
        "name": "few events",
        "n": 200,
        "p": 1.16,
        "c": 0.50,
        "second": 0.0,
        "sequences": 1000,
    },
    {
        "name": "many events",
        "n": 5000,
        "p": 1.16,
        "c": 0.50,
        "second": 0.0,
        "sequences": 500,
    },
    {
        "name": "slow decay, small c",
        "n": 1500,
        "p": 0.90,
        "c": 0.02,
        "second": 0.0,
        "sequences": 1000,
    },
    {
        "name": "fast decay, large c",
        "n": 1500,
        "p": 1.40,
        "c": 2.00,
        "second": 0.0,
        "sequences": 1000,
    },
    {
        "name": "second sequence",
        "n": 1500,
        "p": 1.16,
        "c": 0.50,
        "second": 0.20,
        "second_start": 1.0,
        "second_p": 1.10,
        "second_c": 0.05,
        "sequences": 500,
    },
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


def draw_sequence(config: dict, rng: np.random.Generator) -> np.ndarray:
    """Occurrence times in days after the mainshock for one simulated sequence."""
    n_second = round(config["n"] * config["second"])
    first = omori_sample(
        config["n"] - n_second, config["p"], config["c"], WINDOW_DAYS, 0.0, rng
    )
    if not n_second:
        return np.sort(first)
    start = config["second_start"]
    second = start + omori_sample(
        n_second,
        config["second_p"],
        config["second_c"],
        WINDOW_DAYS - start,
        0.0,
        rng,
    )
    return np.sort(np.concatenate([first, second]))


def one_sequence(job: tuple[int, int, int]) -> dict:
    """Simulate, fit and test one sequence. Returns the numbers worth keeping."""
    config_index, index, replicates = job
    config = CONFIGS[config_index]
    seed = 1_000_003 * (config_index + 1) + index
    rng = np.random.default_rng(seed)
    times = draw_sequence(config, rng)
    try:
        fit = fit_omori(times)
    except ValueError:
        return {"failed": True}
    calibrated = omori_fit_test(times, fit, n_simulations=replicates, seed=seed)
    asymptotic = omori_fit_test(times, fit, n_simulations=0)
    return {
        "failed": False,
        "ks": calibrated.statistic,
        "p_calibrated": calibrated.p_value,
        "p_asymptotic": asymptotic.p_value,
        "fit_p": fit.p,
        "fit_c": fit.c,
    }


def summarise(config: dict, rows: list[dict], replicates: int, seconds: float) -> dict:
    done = [r for r in rows if not r["failed"]]
    n = len(done)
    p_cal = np.array([r["p_calibrated"] for r in done])
    p_asy = np.array([r["p_asymptotic"] for r in done])
    out = {
        "config": {k: v for k, v in config.items() if k != "sequences"},
        "sequences": n,
        "fit_failures": len(rows) - n,
        "replicates": replicates,
        "mean_fit_p": float(np.mean([r["fit_p"] for r in done])),
        "median_fit_c": float(np.median([r["fit_c"] for r in done])),
        "seconds": round(seconds, 1),
        "calibrated": {"mean_p": float(p_cal.mean()), "rejection": {}},
        "asymptotic": {"mean_p": float(p_asy.mean()), "rejection": {}},
    }
    for label, p in (("calibrated", p_cal), ("asymptotic", p_asy)):
        for level in LEVELS:
            k = int(np.sum(p <= level + 1e-12))
            lo, hi = wilson(k, n)
            out[label]["rejection"][f"{level:.2f}"] = {
                "rate": k / n,
                "rejected": k,
                "ci95": [lo, hi],
            }
    return out


def show(result: dict) -> None:
    c = result["config"]
    print(
        f"\n{c['name']}: n={c['n']} p={c['p']} c={c['c']}"
        + (f", {c['second']:.0%} in a second sequence" if c["second"] else "")
        + f"  ({result['sequences']} sequences, {result['replicates']} replicates, "
        f"{result['seconds']:.0f} s)"
    )
    for label in ("asymptotic", "calibrated"):
        block = result[label]
        cells = []
        for level, cell in block["rejection"].items():
            lo, hi = cell["ci95"]
            cells.append(f"{level}: {cell['rate']:.1%} [{lo:.1%}, {hi:.1%}]")
        print(f"  {label:10s} mean p {block['mean_p']:.3f}   " + "   ".join(cells))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--quick", action="store_true", help="40 sequences, 49 replicates"
    )
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--out", type=Path, default=OUT_PATH)
    args = parser.parse_args()

    replicates = 49 if args.quick else REPLICATES
    previous = {}
    if args.out.exists() and not args.quick:
        previous = {
            r["config"]["name"]: r
            for r in json.loads(args.out.read_text(encoding="utf-8"))["results"]
        }

    results = []
    with Pool(args.workers) as pool:
        for ci, config in enumerate(CONFIGS):
            sequences = 40 if args.quick else config["sequences"]
            cached = previous.get(config["name"])
            if (
                cached
                and cached["sequences"] + cached["fit_failures"] == sequences
                and cached["replicates"] == replicates
            ):
                results.append(cached)
                show(cached)
                continue
            started = time.time()
            rows = pool.map(
                one_sequence,
                [(ci, i, replicates) for i in range(sequences)],
                chunksize=4,
            )
            result = summarise(config, rows, replicates, time.time() - started)
            results.append(result)
            show(result)
            if not args.quick:
                write(args.out, results, replicates)
    if not args.quick:
        write(args.out, results, replicates)
        print(f"\nwrote {args.out}")


def write(path: Path, results: list[dict], replicates: int) -> None:
    payload = {
        "meta": {
            "tremor_lab_version": __version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "python": platform.python_version(),
            "replicates": replicates,
            "window_days": WINDOW_DAYS,
        },
        "results": results,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
