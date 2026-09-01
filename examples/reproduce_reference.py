"""Reproduce the reference values for the 2023 Kahramanmaras sequence.

    python examples/reproduce_reference.py

Prints each reported value beside the one this installation computes, and exits
non-zero if any of them disagree, so the check can be run unattended.

The reported values come from the validated Python/SciPy reference implementation and
were reproduced independently by the Google Apps Script implementation, whose optimiser
is a hand-written Nelder-Mead simplex.
"""

import sys
from pathlib import Path

import pandas as pd

from tremor_lab import __version__, analyze_case

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "tests" / "data" / "kahramanmaras_180d.csv"

MAINSHOCK = {"t": "2023-02-06 01:17:32", "lat": 37.1578, "lon": 36.8092, "mw": 7.8}
WINDOW_DAYS = 180
THRESHOLD = 3.5

# quantity, reported value, tolerance, decimals to display
REPORTED = [
    ("events in window", 3469, 0, 0),
    ("Mc (maximum curvature)", 3.4, 0, 1),
    ("b-value", 0.844, 0.0005, 3),
    ("b-value sigma", 0.019, 0.0005, 3),
    ("events above M 3.5", 1529, 0, 0),
    ("Omori p", 1.161, 0.001, 3),
    ("Omori c", 0.497, 0.001, 3),
    ("Omori k", 359, 1, 1),
]


def main() -> int:
    result = analyze_case(
        pd.read_csv(CATALOG),
        MAINSHOCK,
        mc_threshold=THRESHOLD,
        window_days=WINDOW_DAYS,
    )
    computed = [
        result["n_events"],
        result["mc"],
        result["b_value"].b,
        result["b_value"].sigma,
        result["b_value"].n,
        result["omori"].p,
        result["omori"].c,
        result["omori"].k,
    ]

    print(f"Tremor Lab {__version__}")
    print(
        f"catalogue: {CATALOG.name}, {WINDOW_DAYS}-day window, threshold M {THRESHOLD}"
    )
    print(f"\n{'quantity':<24}{'reported':>10}{'this run':>12}   agreement")
    agree = True
    for (name, reported, tolerance, decimals), value in zip(
        REPORTED, computed, strict=True
    ):
        matches = abs(value - reported) <= tolerance
        agree &= matches
        verdict = "ok" if matches else "DIFFERS"
        print(f"{name:<24}{reported:>10}{value:>12.{decimals}f}   {verdict}")

    spread = result["omori_bootstrap"]
    print(
        f"\nbootstrap standard errors over {spread.n_boot} resamples: "
        f"p +/- {spread.p_std:.3f}, c +/- {spread.c_std:.3f}"
    )
    print(
        f"mainshock Mw {MAINSHOCK['mw']}: energy {result['mainshock_energy_j']:.3e} J, "
        f"Bath expectation M {result['bath_mag']:.2f}"
    )
    print(
        "\nAll reported values reproduced."
        if agree
        else "\nSome values differ from those reported."
    )
    return 0 if agree else 1


if __name__ == "__main__":
    sys.exit(main())
