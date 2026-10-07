"""Break the threshold-dependence code seventeen ways; a test must fail each time.

    python examples/mutation_check_threshold.py

A passing suite is evidence of correctness only once it has been shown that it can
fail. For `b_shift` and `b_plateau` that is done here: each mutation below changes one
line of `src/tremor_lab/bvalue.py` in a scratch copy of the package, the tests in
`tests/test_threshold_dependence.py` run against the copy, and the mutation counts as
caught if at least one of them fails. The repository itself is never modified.

The mutations are the plausible mistakes: the wrong error for a difference of nested
estimates (one end's error, a missing subtraction, a sum where a difference belongs), a
one-sided p-value, a bootstrap that resamples half the sample or ignores its seed or its
confidence level, a flipped sign, a statistic that is not squared, a plateau taken from
the top instead of the bottom, a simulated p-value without its plus one, a null drawn
at the wrong scale or with too few events, and off-by-one errors in the window rules.
Each is applied by exact string replacement, and the script stops with an error if the
line it expects has changed, so it cannot report a pass on code it no longer matches.

It takes a few minutes. Exit status is 0 only if every mutation was caught.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / "src"
TARGET = Path("tremor_lab") / "bvalue.py"
TESTS = "tests/test_threshold_dependence.py"

MUTATIONS = {
    "shift: error taken from one end only": (
        "se = low.sigma * float(np.sqrt(low.n / high.n - 1.0))",
        "se = low.sigma",
    ),
    "shift: the subtraction of one forgotten": (
        "se = low.sigma * float(np.sqrt(low.n / high.n - 1.0))",
        "se = low.sigma * float(np.sqrt(low.n / high.n))",
    ),
    "shift: p-value one-sided": (
        "p_value = float(2.0 * norm.sf(abs(z)))",
        "p_value = float(norm.sf(abs(z)))",
    ),
    "shift: bootstrap resamples half the sample": (
        "resample = sample[rng.integers(0, sample.size, sample.size)]",
        "resample = sample[rng.integers(0, sample.size, sample.size // 2)]",
    ),
    "shift: confidence level ignored": (
        "tail = 100.0 * (1.0 - confidence) / 2.0",
        "tail = 2.5",
    ),
    "shift: bootstrap shift has the wrong sign": (
        "shifts[i] = (1.0 / gap_high - 1.0 / gap_low) / np.log(10)",
        "shifts[i] = (1.0 / gap_low - 1.0 / gap_high) / np.log(10)",
    ),
    "shift: bootstrap seed ignored": (
        "rng = np.random.default_rng(seed)\n        shifts",
        "rng = np.random.default_rng(0)\n        shifts",
    ),
    "plateau: step not squared": (
        "step_b**2 / np.where(informative, step_v, 1.0)",
        "np.abs(step_b) / np.where(informative, step_v, 1.0)",
    ),
    "plateau: variance of a step summed, not differenced": (
        "step_v = 1.0 / n[1:] - 1.0 / n[:-1]\n"
        "        informative = step_v > 0.0\n"
        "        terms = np.where(\n"
        "            informative, step_b**2",
        "step_v = 1.0 / n[1:] + 1.0 / n[:-1]\n"
        "        informative = step_v > 0.0\n"
        "        terms = np.where(\n"
        "            informative, step_b**2",
    ),
    "plateau: onset taken from the top": (
        "onset = next((i for i in tested if p[i] >= alpha), None)",
        "onset = next((i for i in reversed(tested) if p[i] >= alpha), None)",
    ),
    "plateau: simulated p-value without the plus one": (
        "return (1.0 + exceeded) / (1.0 + used)",
        "return exceeded / used if exceeded else 0.0",
    ),
    "plateau: obvious-rejection skip swallows candidates": (
        "_OBVIOUS_REJECTION = 1e-10",
        "_OBVIOUS_REJECTION = 0.5",
    ),
    "plateau: coefficient ignored": (
        "sigma = curve.sigma * coefficient / constants.SHI_BOLT_K",
        "sigma = curve.sigma",
    ),
    "plateau: min_window off by one": (
        "tested = [i for i in range(count) if dof[i] >= min_window]",
        "tested = [i for i in range(count) if dof[i] > min_window]",
    ),
    "plateau: null drawn with a quarter of the events": (
        "x = floors[0] + mean_excess * rng.exponential(size=n0)",
        "x = floors[0] + mean_excess * rng.exponential(size=max(n0 // 4, 2))",
    ),
    "plateau: null drawn at the wrong scale": (
        "x = floors[0] + mean_excess * rng.exponential(size=n0)",
        "x = floors[0] + 0.5 * mean_excess * rng.exponential(size=n0)",
    ),
    "plateau: steps above a candidate counted from the wrong end": (
        "tail_steps = np.cumsum(informative[::-1].astype(np.int64))[::-1]",
        "tail_steps = np.cumsum(informative.astype(np.int64))",
    ),
}


def main() -> int:
    original = (SOURCE / TARGET).read_text(encoding="utf-8")
    caught = 0
    survivors = []
    with tempfile.TemporaryDirectory() as scratch:
        copy = Path(scratch) / "src"
        for name, (old, new) in MUTATIONS.items():
            if original.count(old) != 1:
                print(
                    f"cannot apply {name!r}: the line it expects appears "
                    f"{original.count(old)} times in {TARGET}",
                    file=sys.stderr,
                )
                return 2
            shutil.rmtree(copy, ignore_errors=True)
            shutil.copytree(SOURCE, copy)
            (copy / TARGET).write_text(original.replace(old, new), encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "pytest",
                    TESTS,
                    "-x",
                    "-q",
                    "-p",
                    "no:cacheprovider",
                ],
                cwd=REPO,
                env={**os.environ, "PYTHONPATH": str(copy)},
                capture_output=True,
                text=True,
            )
            failed = [
                ln for ln in result.stdout.splitlines() if ln.startswith("FAILED")
            ]
            if result.returncode != 0:
                caught += 1
                print(f"caught    {name}  [{failed[0][7:80] if failed else 'error'}]")
            else:
                survivors.append(name)
                print(f"SURVIVED  {name}")
    print(f"\n{caught} of {len(MUTATIONS)} mutations caught")
    return 0 if not survivors else 1


if __name__ == "__main__":
    raise SystemExit(main())
