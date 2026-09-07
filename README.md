# Tremor Lab

Aftershock-sequence statistics in Python: the magnitude of completeness, the
Gutenberg-Richter b-value, the modified Omori-Utsu decay, radiated energy and the Bath
energy screen, and magnitude homogenisation to Mw.

Tremor Lab is the scripting counterpart of `SeismoSheet.gs`, a spreadsheet
implementation of the same estimators. The two were written independently and use
different optimisers, so their agreement on a real sequence is a cross-check of both.
This package reproduces the values reported for the 2023 Kahramanmaras sequence with one
command.

## Installation

Python 3.11 or later. From the project folder:

```
py -m venv .venv
.venv\Scripts\python.exe -m pip install -e ".[dev]"
```

Runtime dependencies are NumPy, SciPy and pandas. The development extra adds pytest and
ruff.

## Reproducing the reference values

```
.venv\Scripts\python.exe examples\reproduce_reference.py
```

```
quantity                  reported    this run   agreement
events in window              3469        3469   ok
Mc (maximum curvature)         3.4         3.4   ok
b-value                      0.844       0.844   ok
b-value sigma                0.019       0.019   ok
events above M 3.5            1529        1529   ok
Omori p                      1.161       1.161   ok
Omori c                      0.497       0.497   ok
Omori k                        359       358.9   ok
```

The script exits non-zero if any value disagrees.

That command reads a file that already carries elapsed days, so it exercises the
estimators but not the catalogue pipeline. `examples/rebuild_koeri_fixture.py` reaches
the same values the long way: it reconstructs KOERI-shaped timestamps at whole-second
precision, adds ten events outside the window that must be excluded, and runs the whole
path — parsing, homogenisation, windowing, estimation — to land on 3,469 events, Mc 3.4,
b 0.844 and p 1.161. `tests/test_regression.py` locks that path too.

It is a pipeline test, not provenance: the timestamps are derived from the elapsed days,
so the file cannot corroborate the elapsed days themselves. Only the original KOERI
download can, and it should be archived with the release.

## Use from Python

```python
import pandas as pd
from tremor_lab import analyze_case, read_catalog

mainshock = {"t": "2023-02-06 01:17:32", "lat": 37.1578, "lon": 36.8092, "mw": 7.8}

catalog = read_catalog(
    "koeri_export.csv",
    columns={"date": "Tarih", "time": "Saat", "lat": "Enlem",
             "lon": "Boylam", "mag": "Mag"},
    mainshock=mainshock,
    window_days=180,
    mag_type_col="Tip",
)

result = analyze_case(catalog, mainshock, mc_threshold=3.5)
print(result["mc"], result["b_value"], result["omori"])
```

Individual estimators can be used on their own:

```python
from tremor_lab import b_value_aki, energy_joules, fit_omori, mc_maxcurvature

mc = mc_maxcurvature(mags)                  # 3.4
b, sigma, n = b_value_aki(mags, mc=3.5)     # 0.844, 0.019, 1529
p, c, k, n = fit_omori(times_days)          # 1.161, 0.497, 358.9, 1529
energy_joules(7.8)                          # 3.16e16 joules
```

## Use from the command line

```
.venv\Scripts\tremor-lab.exe run examples\kahramanmaras.toml
```

The settings file describes the catalogue, its columns, the mainshock, the analysis
options and any constant to override, so a run is a file that can be archived beside the
paper rather than a sequence of arguments someone has to remember. An unknown key or an
unknown constant name is refused rather than silently ignored.

```toml
[catalog]
path = "catalogue.csv"          # forward slashes on Windows, or single quotes
dayfirst = false                # true for European dates such as 06.02.2023
# radius_km = 100               # optional distance limit from the epicentre;
                                # omitted means no spatial cut, and the report
                                # says so either way

[catalog.columns]
date = "Tarih"                  # KOERI: separate date and time columns
time = "Saat"
lat  = "Enlem"
lon  = "Boylam"
mag  = "Mag"
mag_type = "Tip"                # optional; enables homogenisation to Mw

[mainshock]
t = "2023-02-06 01:17:32"       # on the same clock as the catalogue
lat = 37.1578
lon = 36.8092
mw = 7.8

[analysis]
window_days = 180
threshold = 3.5                 # omit to use the estimated Mc
n_boot = 200
n_fit_simulations = 600         # replicates behind the decay fit test's p-value;
                                # raise it when the verdict says borderline
seed = 0

[constants]                     # optional; any name from tremor_lab.constants
DELTA_MB = 1.15
```

For a USGS export, replace the `date` and `time` entries with a single
`datetime = "time"`. For a catalogue that already carries elapsed days rather than
timestamps, use `dt_days = "dt_days"` and `mag = "mw"`.

Add `--fmd-out table.csv` to write the frequency-magnitude table for a Gutenberg-Richter
plot.

## Use from a browser

`web/tremor-lab.html` is the whole tool in one file. Open it by double-clicking:
it runs offline, with nothing installed and no server. Drop a catalogue in, or press
the button to load the bundled 2023 Kahramanmaras sequence.

It reports the same quantities as the package, draws the sequence, the
frequency-magnitude distribution and the decay curve, and shows the data behind each
chart. Eight settings are editable on the page - the window, the threshold, the bin
width, the Mc correction, the Bath deficit, the bootstrap count, the fit-test replicate
count and the seed - and the page states which columns it read, how magnitudes were treated, and what it assumed, so
a number never travels without the choices that produced it. The remaining constants,
including the Scordilis coefficients and the Earth radius, are fixed in the page; change
those in the Python package, which is the authority for published values in any case.

A self-test runs on load and reproduces the eight reference values below in the browser.

The page runs the JavaScript implementation transplanted from `SeismoSheet.gs`. This
package remains the authority for published values; the two agree on the reference
catalogue and on the real KOERI and USGS exports they have been checked against.

## What it computes

**Magnitude of completeness**, maximum-curvature method: the mode of the incremental
frequency-magnitude distribution plus a correction, +0.2 by default (Wiemer and Wyss
2000). Magnitudes are assigned to the nearest bin centre rather than to a half-open
interval, because reported magnitudes sit on a decimal grid that binary floating point
cannot represent exactly. On the reference catalogue this agrees with half-open binning.

**Gutenberg-Richter b-value**, Aki (1965) maximum-likelihood estimator on events at or
above Mc, with the Shi and Bolt (1982) standard error:

    b = 1 / (ln10 (mean(M) - (Mc - dm/2)))
    sigma = 2.30 b^2 sqrt(sum (M - mean(M))^2 / (n (n - 1)))

The half-bin offset is the lower edge of the completeness bin: a magnitude reported as
Mc stands for the interval Mc +/- dm/2.

**Modified Omori-Utsu decay**, n(t) = k / (c + t)^p, by maximum likelihood on unbinned
post-mainshock times (Ogata 1983; Utsu, Ogata and Matsuura 1995), minimised with
`scipy.optimize.minimize`. Only c and p are searched: at any (c, p) the likelihood is
maximised by k = n / integral, so k follows in closed form. This removes one dimension
and any dependence on a starting value for k. The test suite asserts the result is at
least as likely as a full three-parameter search, and independent of the starting point.
Uncertainty on p and c comes from a bootstrap that holds the observation interval fixed.

By default the observation interval ends at the last event supplied, not at the nominal
window length. On the reference catalogue those differ, 179.596 days against 180, and
the choice moves p in the third decimal.

**Radiated energy**, log10 E = 1.5 M + 4.8 with E in joules (Gutenberg and Richter 1956;
Kanamori 1977), and the Bath energy ratio 10^(1.5 (M_sec - (M_main - delta_mb))) with
delta_mb defaulting to 1.15 (Bath 1965).

**Magnitude homogenisation to Mw** by the Scordilis (2006) global relations, Ms in two
branches and mb linearly. Labels beginning "mw" pass through; ml, md and blank labels
are taken as reported, since no global relation is published for them.

## Changing the constants

Every published constant is a default, not a fixed law of the package. All 24 live in
`tremor_lab.constants`, including `N_FIT_SIMULATIONS`, the number of replicates used to
simulate the null distribution of the decay fit test.

For a single call, pass the keyword argument:

```python
mc_maxcurvature(mags, dm=0.05, correction=0.0)
bath_ratio(7.5, 7.8, delta_mb=1.2)
```

For a whole session, reassign the name. Functions read these at call time, so the change
reaches everything afterwards, including calls made inside `analyze_case`:

```python
from tremor_lab import constants
constants.DELTA_MB = 1.2
```

Without writing Python, set it in the `[constants]` section of a settings file.

Three of the 24 are the exception to the keyword route: `OMORI_C_FLOOR`, `OMORI_P_MIN`
and `OMORI_P_MAX` judge whether a fitted decay is worth believing rather than entering
the fit, so they have no keyword argument and are changed by reassignment or from a
settings file. The other 21 can be set all three ways.

## Validation

Correctness is demonstrated by tests that run, not by reading the code:

```
.venv\Scripts\python.exe -m pytest
```

Three layers. **Analytic**, where the expected value is exact arithmetic on a published
relation: the energy exponent, the Bath ratio at a one-magnitude anomaly, both Scordilis
conversions, a great-circle distance. **Synthetic recovery**, where magnitudes are drawn
from a Gutenberg-Richter law with a known b, and occurrence times from an Omori process
with known p, c and k by inverse-CDF sampling, and the estimators must recover them
across several seeds. **Regression**, which reproduces the reference values below on the
bundled catalogue.

### Checked against an independent implementation

`examples/compare_with_seismostats.py` runs the same catalogue through
[seismostats](https://pypi.org/project/seismostats/), an independent package from the
Swiss Seismological Service at ETH Zurich. It is not a dependency; install it in a
separate environment to reproduce the comparison.

| quantity | Tremor Lab | seismostats | |
| --- | --- | --- | --- |
| maximum-curvature Mc | 3.4 | 3.4 | identical |
| events above M 3.5 | 1529 | 1529 | identical |
| b-value | 0.844132 | 0.846804 | 0.32% apart |
| completeness by b-stability | 4.1 | 4.0 | one bin apart |

Every difference is a documented difference of convention, not a disagreement about
the data, and the script states which one accounts for each.

The b-value gap is the important one. Tremor Lab reports the **Aki (1965) estimator
with Utsu's half-bin offset**, which is what the thesis and the spreadsheet use.
seismostats reports the **Tinti and Mulargia (1987)** estimator, the exact
maximum-likelihood solution for magnitudes reported on a grid. The half-bin form is the
first-order approximation of it, and the two converge as the bin narrows: they differ by
0.32% at dM 0.1 and by 0.004% at dM 0.01. `b_value_tinti` implements the exact estimator
and reproduces seismostats to six decimals, so the choice can be tested rather than
assumed. Standard errors are not compared directly: both use Shi and Bolt (1982), but
seismostats takes ln(10) as the coefficient where this package takes the published
rounding of 2.30, which is `constants.SHI_BOLT_K` and adjustable.

### Reference values

2023 Kahramanmaras sequence, KOERI catalogue, 180-day window, threshold M >= 3.5:

| quantity | value |
| --- | --- |
| events in window | 3469 |
| Mc (maximum curvature) | 3.4 |
| b-value (Aki) | 0.844 +/- 0.019 on 1529 events |
| Omori | p = 1.161, c = 0.497, k = 359 |

Produced with Python 3.13.7, NumPy 2.5.2, SciPy 1.18.1, pandas 2.3.3 on Windows 11, 5 September 2026.

These were produced by the validated Python/SciPy reference implementation and
reproduced independently by the spreadsheet implementation, whose optimiser is a
hand-written Nelder-Mead simplex rather than SciPy's. Two independent optimisers
agreeing to three decimals on a real sequence is what makes them reference values. If a
change makes `tests/test_regression.py` fail, the change is wrong, not the numbers.

## Limitations

The estimators are standard and inherit their known behaviour. Maximum-curvature Mc is
sensitive to binning and to short-term aftershock incompleteness in the first hours
after a mainshock. The Aki b-value assumes completeness above the chosen threshold. The
Omori fit needs a dense catalogue and is not forced on sparse historical sequences:
where the catalogue above threshold is too thin, the absence is reported rather than a
fit produced. Of the three decay parameters, k is the least well constrained; it is
derived from the fitted c and p rather than searched for, and its sampling error stays
at a few per cent where p is already recovered to a fraction of one.

The Scordilis relations are applied without range checks, matching the spreadsheet
implementation, and the two Ms branches are mildly discontinuous across the uncalibrated
6.1 to 6.2 gap.

This package computes the estimators and the energy diagnostics. The probabilistic
largest-secondary anomaly statistic is a separate layer, added once its specification is
frozen, and is deliberately not part of this foundation.

## When you are ready to publish

Two steps remain before the software is formally citable.

**Continuous integration.** Create a GitHub repository, push, and add
`.github/workflows/tests.yml` running `ruff check`, `ruff format --check` and `pytest`
on every push, so a green badge in this README is evidence rather than assertion.

**Archive for a DOI.** Enable the repository in Zenodo, tag a release and push the tag,
and Zenodo mints a versioned DOI. Put that DOI in `CITATION.cff` and in the paper's
software availability statement.

## Licence and citation

MIT, see `LICENSE`. Citation metadata is in `CITATION.cff`.

## References

Aki, K. (1965). Bath, M. (1965). Efron, B. (1979). Gutenberg, B. and Richter, C. F.
(1944, 1956). Kanamori, H. (1977). Nelder, J. A. and Mead, R. (1965). Ogata, Y. (1983).
Scordilis, E. M. (2006). Shi, Y. and Bolt, B. A. (1982). Utsu, T., Ogata, Y. and
Matsuura, R. S. (1995). Wiemer, S. and Wyss, M. (2000).

Verify each primary source before final citation.
