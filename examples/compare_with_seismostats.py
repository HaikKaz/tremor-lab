"""Compare Tremor Lab against an independent third-party implementation.

    python examples/compare_with_seismostats.py

`seismostats` is an independent package from the Swiss Seismological Service at
ETH Zurich, implementing the same estimators from the same literature. It is not
a dependency of Tremor Lab; install it into a separate environment to run this:

    py -m venv .venv-compare
    .venv-compare\\Scripts\\python.exe -m pip install seismostats
    .venv-compare\\Scripts\\python.exe -m pip install -e .
    .venv-compare\\Scripts\\python.exe examples/compare_with_seismostats.py

Every difference this prints is a documented difference of convention between two
published estimators, not a disagreement about the data. The script names the
convention behind each one and exits non-zero if any gap exceeds what that
convention predicts.
"""

import sys
from pathlib import Path

import pandas as pd

from tremor_lab.bvalue import b_value_aki, b_value_tinti, mc_b_stability
from tremor_lab.completeness import mc_maxcurvature

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "tests" / "data" / "kahramanmaras_180d.csv"
THRESHOLD = 3.5
DM = 0.1


def main() -> int:
    try:
        from seismostats.analysis import (
            estimate_b,
            estimate_mc_b_stability,
            estimate_mc_maxc,
        )
    except ImportError:
        print(__doc__)
        print("seismostats is not installed here; nothing to compare.")
        # Skipping is right on a machine that simply does not have the
        # comparator, and wrong in CI, where it is installed on purpose: a job
        # whose whole point is the third-party comparison reported success
        # having compared nothing at all. --required turns the skip into a
        # failure, and the workflow passes it.
        if "--required" in sys.argv:
            print(
                "but --required was given, so this counts as a failure: the "
                "comparison this run exists to perform did not happen"
            )
            return 1
        return 0

    mags = pd.read_csv(CATALOG)["mw"].to_numpy()
    print(f"catalogue: {CATALOG.name}, {mags.size} magnitudes\n")
    ok = True

    # ---------------------------------------------- completeness magnitude --
    mine_mc = mc_maxcurvature(mags, dm=DM)
    theirs_mc, _ = estimate_mc_maxc(mags, fmd_bin=DM, correction_factor=0.2)
    agree = abs(mine_mc - float(theirs_mc)) < 1e-9
    ok &= agree
    print("Maximum-curvature completeness (Wiemer and Wyss 2000)")
    print(f"  Tremor Lab   {mine_mc}")
    print(f"  seismostats  {round(float(theirs_mc), 10)}")
    print(f"  {'identical' if agree else 'DIFFER'}\n")

    # ------------------------------------------------------------ b-value --
    mine_aki = b_value_aki(mags, THRESHOLD, dm=DM)
    mine_tinti = b_value_tinti(mags, THRESHOLD, dm=DM)
    theirs_b, _std, theirs_n = estimate_b(
        mags, mc=THRESHOLD, delta_m=DM, return_std=True, return_n=True
    )
    same_n = mine_aki.n == theirs_n
    reproduces = abs(mine_tinti.b - theirs_b) < 1e-6
    ok &= same_n and reproduces
    apart = abs(mine_aki.b / mine_tinti.b - 1) * 100

    print(f"b-value at Mc {THRESHOLD}")
    print(f"  Aki 1965 + Utsu half-bin, ours    b {mine_aki.b:.6f}  n {mine_aki.n}")
    print(f"  Tinti and Mulargia 1987, ours     b {mine_tinti.b:.6f}  n {mine_tinti.n}")
    print(f"  Tinti and Mulargia, seismostats   b {theirs_b:.6f}  n {theirs_n}")
    print(f"  same sample                       {'yes' if same_n else 'NO'}")
    print(f"  exact MLE reproduced to 1e-6      {'yes' if reproduces else 'NO'}")
    print(f"  half-bin form against exact MLE   {apart:.3f} per cent apart\n")

    print("The half-bin form is the first-order approximation of the exact")
    print("binned MLE, so the gap must fall with the bin width:")
    for dm in (0.4, 0.2, 0.1, 0.05, 0.01):
        a = b_value_aki(mags, THRESHOLD, dm=dm).b
        t = b_value_tinti(mags, THRESHOLD, dm=dm).b
        print(f"    dM {dm:<5} Aki {a:.6f}   Tinti {t:.6f}   gap {abs(a - t):.2e}")
    print()

    # -------------------------------------------- completeness by stability --
    mine_stab = mc_b_stability(mags, dm=DM)
    theirs_stab, _ = estimate_mc_b_stability(mags, delta_m=DM, stability_range=0.5)
    close = theirs_stab is not None and abs(mine_stab - float(theirs_stab)) <= DM + 1e-9
    ok &= close
    verdict = "within one bin" if close else "DIFFER BY MORE THAN ONE BIN"
    print("Completeness by b-value stability (Cao and Gao 2002)")
    print(f"  Tremor Lab   {mine_stab}")
    print(f"  seismostats  {theirs_stab}")
    print(f"  {verdict}. The criterion walks up a curve of b-values that the")
    print("  two estimators place differently, so this follows from the")
    print("  b-value convention above rather than being independent of it.\n")

    print("Standard errors are not compared directly. Both use Shi and Bolt")
    print("(1982), but seismostats takes ln(10) as the coefficient where Tremor")
    print("Lab takes the published rounding of 2.30, which is SHI_BOLT_K and")
    print("adjustable.")

    print()
    print(
        "Every comparison behaved as its convention predicts."
        if ok
        else "A comparison fell outside its expected tolerance."
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
