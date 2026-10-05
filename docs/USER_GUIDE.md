# Tremor Lab user guide

This guide is the reference for running Tremor Lab on a catalogue of your own: what goes
in, what comes out, what each setting does, and what the package does when the input is
awkward. The [README](../README.md) is the overview and carries the validation. The
[design record](PROVENANCE_AND_VALIDATION.md) explains why each estimator is written the
way it is.

Every Python example in this guide is executed by `tests/test_user_guide.py`, so the code
on this page runs against the version of the package it ships with.

- [1. Requirements and installation](#1-requirements-and-installation)
- [2. Choosing how to run it](#2-choosing-how-to-run-it)
- [3. Inputs](#3-inputs)
- [4. The settings file](#4-the-settings-file)
- [5. Outputs](#5-outputs)
- [6. What to expect](#6-what-to-expect)
- [7. How-to recipes](#7-how-to-recipes)
- [8. Messages and what they mean](#8-messages-and-what-they-mean)
- [Appendix: the published constants](#appendix-the-published-constants)

## 1. Requirements and installation

Python 3.11 or later, with NumPy 2, SciPy 1.13 or later and pandas 2.2 or later; the exact
ranges are in `pyproject.toml`. The package is pure Python, has no compiled code of its
own, needs no network access when it runs, and uses one core.

Cost is set by two things: the number of events above the threshold, and the number of
replicates behind the decay-fit test and the bootstrap. With the default settings (200
bootstrap resamples, 600 test replicates) the three catalogues in `examples/` and
`tests/data/` ran in about 5 seconds (3,469 events), 8 seconds (13,524 events) and 11
seconds (28,962 events) on one core of a cloud container, and the largest used about
150 MB of memory. Time grows linearly with the replicate count.

Clone the repository, create an environment and install the package in it.

Windows:

```
git clone https://github.com/HaikKaz/tremor-lab.git
cd tremor-lab
py -m venv .venv
.venv\Scripts\python.exe -m pip install -e .
```

Linux and macOS:

```
git clone https://github.com/HaikKaz/tremor-lab.git
cd tremor-lab
python3 -m venv .venv
.venv/bin/python -m pip install -e .
```

Install `-e ".[dev]"` instead to add pytest and ruff for development. The two scripts
that draw the paper's figures also need matplotlib, which is not a dependency of the
package: `pip install matplotlib`.

Check the installation with three commands, run from the repository folder with the
environment active (or with the full path to `.venv\Scripts\tremor-lab.exe` and
`.venv\Scripts\python.exe` on Windows):

```
tremor-lab --version
python examples/reproduce_reference.py
python -m pytest -q
```

The first prints `Tremor Lab 1.1.2`. The second recomputes the eight reference values for
the 2023 Kahramanmaras sequence and ends with `All reported values reproduced.`; it exits
non-zero if any value disagrees. The third runs the test suite.

## 2. Choosing how to run it

| You want | Use |
| --- | --- |
| A run you can archive beside a paper and repeat exactly | the command line: `tremor-lab run settings.toml` |
| The estimators inside your own script or notebook | the Python API: `read_catalog`, `analyze_case` and the individual estimators |
| A quick look at one file with nothing installed | `web/tremor-lab.html`, opened in a browser |

The browser page runs the same estimators written independently in JavaScript. The Python
package is the authority for published values; where the two differ, the package is
right and the difference is a bug in the page.

A first run needs no catalogue of your own. The repository carries three:

```
tremor-lab run examples/hectormine.toml
```

prints the report described in section 5 in about eight seconds.

## 3. Inputs

### 3.1 The catalogue file

A CSV file with a header row. Which of three shapes it is depends on which column roles
you name.

| Shape | Roles to name | Typical source |
| --- | --- | --- |
| Separate date and clock columns | `date`, `time`, `lat`, `lon`, `mag` | KOERI exports |
| One timestamp column | `datetime`, `lat`, `lon`, `mag` | USGS ComCat CSV exports |
| Elapsed days already computed | `dt_days`, `mag` | a file windowed elsewhere |

`mag_type` is an optional role on the first two shapes. The roles mean this:

| Role | Meaning |
| --- | --- |
| `date`, `time` | Calendar date and clock time in two columns, joined before parsing. A compact numeric date such as `20230206` is read as written. |
| `datetime` | One timestamp column. ISO strings with a trailing `Z` or a numeric UTC offset are accepted when every row carries the same offset; the offset is dropped, not converted. |
| `lat`, `lon` | Epicentre in decimal degrees. Needed for `radius_km`; without them no distance is computed. |
| `mag` | The magnitude. |
| `mag_type` | The magnitude scale of each row, for homogenisation to Mw (section 3.4). |
| `dt_days` | Days elapsed since the mainshock. The file is used as it stands. |

A row whose time, named position or magnitude cannot be read is dropped, and the report
says how many rows were read and how many were usable. A role you misspell is refused by
name, rather than ignored and left to change the numbers.

### 3.2 The clock

The package converts no time zone. Give the mainshock origin time on the same clock as
the catalogue. USGS exports are in UTC. Check the clock of any other source before you
copy its origin time from somewhere else. A column whose rows carry different UTC offsets
is refused, because there is no single clock to drop.

Dates such as `06.02.2023` are read month-first unless you set `dayfirst = true`. That
flag applies only to rows where the day and month could be confused: a stamp that states
its year first (`2023.02.06`, or an ISO stamp) is read as written, and the run warns how
many rows the flag did not touch.

### 3.3 The mainshock

A mapping with four keys: `t`, the origin time as text; `lat` and `lon` in degrees; and
`mw`, the moment magnitude. `mw` is used only for the radiated energy and the Bath
expectation. `t`, `lat` and `lon` define the window and the distance limit.

Events at or before `t` are excluded, so every analysed time is positive. Give `t` to the
precision of the catalogue's own timestamp for the mainshock. ComCat stamps carry
fractions of a second, and a mainshock time rounded down to the whole second leaves the
mainshock's own record inside its aftershock sample as the first event. One M 7.1 event is
a small share of a sample of thousands, but it sits far above the rest. In the 1999 Hector
Mine catalogue it moves b by about 0.001 at the completeness magnitude and by 0.1 at
M 4, and in the 2019 Ridgecrest catalogue it moves the goodness-of-fit completeness
estimate by one bin. The examples therefore use the timestamp as the catalogue prints it
(`1999-10-16 09:46:44.46`).

### 3.4 Magnitudes

Magnitudes are taken to be moment magnitudes unless you name a `mag_type` column. With one,
each row is homogenised to Mw by the Scordilis (2006) global relations: labels beginning
`mw` pass through, `ms` and `mb` are converted, and every other label (`ml`, `md`, a blank)
is used as reported, because no global relation is published for those scales. Matching
ignores case and surrounding spaces. The report counts what was converted:

```
magnitudes           column 'mag', scale column 'magType': 0 Ms and 3 mb converted to Mw by the Scordilis relations, 13 already Mw, 13510 used as published
```

On the elapsed-days shape nothing is converted, and `mag_type` is refused there rather than
accepted and ignored.

### 3.5 Window and distance limit

The window is `window_days` long (180 by default): an event is analysed if its elapsed
time is greater than zero and no more than the window. `radius_km` keeps only events
within that great-circle distance of the mainshock epicentre. It is off unless you set it,
because a spatial cut changes the result, and the report states which events it removed or
that none were cut.

## 4. The settings file

A TOML file with four sections, `[catalog]` (with a `[catalog.columns]` table inside it),
`[mainshock]`, `[analysis]` and `[constants]`. An unknown section, key or constant is refused by
name; a typo does not run a different analysis in silence. A relative `path` is read
relative to the folder that holds the settings file, so a run does not depend on where you
launch it from.

```toml
[catalog]
path = "hectormine.csv"
radius_km = 100

[catalog.columns]
datetime = "time"
lat = "latitude"
lon = "longitude"
mag = "mag"
mag_type = "magType"

[mainshock]
t = "1999-10-16 09:46:44.46"
lat = 34.6033
lon = -116.265
mw = 7.1

[analysis]
window_days = 180
n_boot = 200
n_fit_simulations = 600
seed = 0
```

### The [catalog] section

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `path` | text | required | The CSV file. Use forward slashes on Windows, or single quotes; a backslash inside double quotes is an escape character in TOML. |
| `dayfirst` | true or false | false | Read ambiguous dates day-first (section 3.2). |
| `radius_km` | number | none | Distance limit from the epicentre. Refused on the elapsed-days shape, which has no positions. |

### The [catalog.columns] section

Maps each role in section 3.1 to the name of a column in your file. `mag` is always
needed, plus either `datetime` or `date` with `time`, or `dt_days`.

### The [mainshock] section

`t`, `lat`, `lon` and `mw`, as in section 3.3. The section is required.

### The [analysis] section

All values are numbers; `true` and `false` are refused, because Python reads `true` as 1.

| Key | Default | Meaning |
| --- | --- | --- |
| `threshold` | the estimated Mc | Magnitude threshold of the b-value and the decay fit. Omit it to use the maximum-curvature estimate. |
| `window_days` | 180 | Aftershock window in days. |
| `dm` | 0.1 | Magnitude bin width. |
| `mc_correction` | 0.2 | Added to the maximum-curvature magnitude. |
| `n_boot` | 200 | Bootstrap resamples behind the uncertainty on p and c. Zero skips them. |
| `seed` | 0 | Seed of the bootstrap and of the decay-fit test, so a run repeats exactly. |
| `min_events_for_mc` | 50 | A window with this many events or fewer gets no completeness estimate. |
| `min_events_for_fit` | 100 | A sample at or above the threshold with this many events or fewer gets no b-value and no decay fit. |
| `n_fit_simulations` | 600 | Replicates behind the p-value of the decay-fit test. Zero falls back to the uncalibrated p-value, which carries no verdict. |

### The [constants] section

Optional. Any name from the appendix, set to a number. It changes a published constant for
that run only; the previous values are put back before the tool returns. The report lists
every constant that differs from its published value.

### Command-line options

```
tremor-lab run settings.toml [--fmd-out table.csv]
tremor-lab --version
```

`--fmd-out` also writes the frequency-magnitude table (section 5.2). Its path is relative
to the folder you run the command from, not to the settings file.

## 5. Outputs

### 5.1 The report

The command line prints one block of text. This is `examples/hectormine.toml` with the
default settings:

```text
Tremor Lab 1.1.2
catalogue            hectormine.csv
rows                 13526 read, 13526 usable
distance limit       100 km from the epicentre, 1 of the 13525 events inside the window removed; farthest kept 100 km
events in window     13524
completeness Mc      1.7  (maximum curvature)
threshold used       1.7  (8465 events in its completeness class, at or above 1.65)
Mc, other methods    goodness of fit 1.5, b stability 2.1
                     the methods span M 1.5 to 2.1; b depends on which is used
b-value              0.840 +/- 0.008  on 8465 events
Omori p              1.087 +/- 0.032
Omori c              4.830 +/- 0.572
Omori k              3106.8
decay fit test       KS 0.0136, p 0.002 +/- 0.002  (NOT consistent with a single Omori decay)
                     p from a parametric bootstrap, 600 replicates
mainshock energy     2.818e+15 J
Bath expectation     M 5.95
magnitudes           column 'mag', scale column 'magType': 0 Ms and 3 mb converted to Mw by the Scordilis relations, 13 already Mw, 13510 used as published
settings             window 180 days; magnitude bins 0.1 wide; Mc correction +0.2
                     threshold M 1.7, the completeness magnitude estimated from the catalogue
                     dates read month-first where the order is ambiguous
                     200 bootstrap resamples; 600 fit-test replicates; seed 0
                     constants at their published defaults
```

| Line | What it says |
| --- | --- |
| `rows` | Rows in the file, and how many could be read. Unreadable rows are counted, not hidden. |
| `distance limit` | Whether a spatial cut was made, how many events it removed, and the farthest event kept. |
| `events in window` | Events with elapsed time in (0, window]. |
| `completeness Mc` | The maximum-curvature estimate: the mode of the incremental distribution plus the correction. |
| `threshold used` | The threshold the b-value and the decay fit were taken at, how many events it keeps, and the lower edge of its bin, which is where the sample starts. |
| `Mc, other methods` | The goodness-of-fit and b-stability estimates, and the span of the three. A wide span means the b-value depends on a choice you have to report. |
| `b-value` | The Aki estimate with the half-bin offset, its Shi and Bolt standard error, and the sample size. |
| `Omori p`, `c`, `k` | The modified Omori-Utsu fit by maximum likelihood. The `+/-` values are bootstrap standard errors; k is derived from p and c and is the least well constrained of the three. |
| `decay fit test` | The Kolmogorov-Smirnov statistic of the time-rescaled residuals, its simulated p-value and the standard error of that p-value, with a verdict (section 6.4). |
| `mainshock energy`, `Bath expectation` | Radiated energy of the mainshock in joules, and the magnitude of the largest aftershock Bath's law expects. |
| `magnitudes` | What happened to the magnitudes (section 3.4). |
| `settings` | Every choice that moved a number, so the block can be pasted into a methods section. |

### 5.2 The frequency-magnitude table

With `--fmd-out table.csv`:

```text
magnitude,incremental,cumulative
1.0,151,13525
1.1,190,13374
1.2,360,13184
```

`magnitude` is the bin centre, `incremental` the events in that bin, and `cumulative` the
events at or above it, which is the quantity a Gutenberg-Richter plot draws against
magnitude. If the window holds no events there is nothing to tabulate, and the run says so
and exits non-zero rather than leave an older file at that path to be plotted as though
it belonged to this sequence.

### 5.3 The Python result

`analyze_case` returns a dictionary:

| Key | Content |
| --- | --- |
| `n_events` | Events in the window. |
| `mc` | Maximum-curvature completeness magnitude, or `None` for a sparse window. |
| `mc_methods` | Dictionary of the three estimates: `maximum_curvature`, `goodness_of_fit`, `b_stability`. A method that finds nothing holds `None`. |
| `threshold`, `n_above`, `sample_floor` | The threshold used, the events it keeps, and the lower edge of its bin. |
| `b_value` | `BValue(b, sigma, n)`, or `None`. |
| `omori` | `Omori(p, c, k, n)`, or `None`. |
| `omori_bootstrap` | `OmoriBootstrap(p_std, c_std, n_boot)`, or `None` when there is no bootstrap. |
| `omori_fit_test` | `FitTest(statistic, p_value, n, method, p_value_se)`, or `None`. |
| `omori_warning` | A sentence when the fit should not be believed (section 6.5), else `None`. |
| `b_stability` | The `BStability` curve of b against threshold, or `None` where none can be built. |
| `b_stability_note` | Why b-stability declined to give an estimate, else `None`. |
| `bootstrap_note` | Present only when the bootstrap could not run, with the reason. |
| `mainshock_energy_j`, `bath_mag` | Energy in joules and the Bath expectation. |
| `fmd` | `FMD(edges, inc, cum)`, or `None` for an empty window. |
| `n_unusable` | Rows dropped for a missing magnitude or time. |
| `note` | Why `b_value` and `omori` are `None`, when they are. |

### 5.4 Exit status

The command line exits 0 when it succeeds and 2 when it refuses. A refusal prints one line
to standard error beginning `tremor-lab:` and says what is wrong and, where it can, which
file or setting is responsible. Section 8 lists the messages.

## 6. What to expect

### 6.1 Sparse windows

A window with `min_events_for_mc` events or fewer gets no completeness estimate. A sample
at or above the threshold with `min_events_for_fit` events or fewer gets no b-value and no
decay fit. In both cases the numbers are `None` and `note` says why. The package does not force a
fit on a sequence too thin to support one. On the command line the report line reads
`b and decay not estimated:` followed by the reason, for example:

```text
b and decay          not estimated: only 13 events at or above M 3.5; a stable fit needs more than 100
```

### 6.2 Three completeness estimates

Maximum curvature finds where the incremental distribution peaks. The goodness-of-fit
method asks whether the data above a candidate look like a Gutenberg-Richter law, and
b-stability asks from where the slope stops moving. They answer different questions and
often disagree. In the Hector Mine run above they give 1.7, 1.5 and 2.1, and b is 0.840 at
the first and 0.987 at the last. The report prints the span for that reason, and a result
that quotes a b-value should quote the threshold beside it.

### 6.3 A method that returns "not estimated"

b-stability can find no stable region over the magnitudes available, or be unable to run
at a bin width wider than its averaging window allows. The report then shows `not
estimated` for that method and carries on with the others. It is a finding about that
catalogue or a setting to change, and it is not an error.

### 6.4 The decay-fit verdict

The Omori curve is fitted to the very times it is then tested against, so the textbook
Kolmogorov p-value is far too generous. The package simulates the null by parametric
bootstrap and reports that p-value with its own standard error. The verdict is one of
four sentences:

| Verdict | Meaning |
| --- | --- |
| `NOT consistent with a single Omori decay` | p plus twice its standard error is below 0.05. |
| `consistent with one Omori decay` | p minus twice its standard error is above 0.05. |
| `borderline; raise the replicate count to decide` | 0.05 lies within two standard errors of p. More replicates settle it (section 7.6). |
| `no verdict: this p-value is uncalibrated` | `n_fit_simulations` was zero. The asymptotic p-value is printed but should not be read. |

A small p-value means the sequence is not one decay. It does not say why; a large
aftershock that starts a sequence of its own is the common cause, and so is a catalogue
whose completeness changes during the window.

### 6.5 The CAUTION line

A fit that converged is still marked when it should not be believed: when the offset c has
collapsed onto the boundary of the model, so no interior maximum exists and c is not
identified, or when p falls outside 0.5 to 2.0, the range compiled from real sequences.
The fit is returned either way, with the sentence under `omori_warning`.

### 6.6 Repeatability

The bootstrap and the decay-fit test are seeded. The same file, settings and version give
the same report. The observation interval of the decay fit ends at the last event above
the threshold and not at the nominal window length; the two differ slightly, and the
choice moves p in the third decimal.

## 7. How-to recipes

### 7.1 Analyse a catalogue of your own from a settings file

Copy the example closest to your file, `examples/kahramanmaras.toml` for elapsed days or
`examples/hectormine.toml` for a USGS export, change the paths, columns and mainshock, and
run it. A raw KOERI export has separate date and time columns:

```toml
[catalog]
path = "koeri_export.csv"

[catalog.columns]
date = "Tarih"
time = "Saat"
lat = "Enlem"
lon = "Boylam"
mag = "Mag"
mag_type = "Tip"

[mainshock]
t = "2023-02-06 01:17:32"      # on the clock the export uses
lat = 37.1578
lon = 36.8092
mw = 7.8

[analysis]
threshold = 3.5
seed = 0
```

`tests/data/kahramanmaras_180d_koeri.csv` is a file of this shape; it names its magnitude
column `xM (Biggest Mag)`. Open your own file and read the real header before you write
the `[catalog.columns]` section, because exports change.

### 7.2 Download a USGS catalogue

The USGS service caps one query at 20,000 events, so a busy sequence is fetched in time
slices and joined. The header of `examples/ridgecrest.toml` holds a working recipe:

```
BASE="https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv"
BASE="$BASE&latitude=35.7695&longitude=-117.5993&maxradiuskm=100"
BASE="$BASE&minmagnitude=1.0&orderby=time-asc"
curl -o p1.csv "$BASE&starttime=2019-07-06T03:19:53&endtime=2019-07-16T00:00:00"
```

Repeat with later slices, then keep one header line and append the rest. Start the first
slice at or before the mainshock, so its own record is in the file, and take the origin
time `t` from that record, fractions of a second included (section 3.3). ComCat revises magnitudes and
adds events, so a download made later will differ slightly from the files in this
repository.

### 7.3 Run the analysis from Python

`read_catalog` reads a file and reduces it to the aftershock window; `analyze_case` runs
the estimators on that table. This is the same analysis as `examples/hectormine.toml`:

```python
from tremor_lab import analyze_case, read_catalog

mainshock = {"t": "1999-10-16 09:46:44.46", "lat": 34.6033, "lon": -116.265, "mw": 7.1}
columns = {"datetime": "time", "lat": "latitude", "lon": "longitude", "mag": "mag"}

catalog = read_catalog(
    "examples/hectormine.csv",
    columns=columns,
    mainshock=mainshock,
    window_days=180,
    mag_type_col="magType",
    radius_km=100,
)
result = analyze_case(catalog, mainshock, seed=0)
print(result["n_events"], result["mc"], result["threshold"])
print(result["b_value"])
print(result["omori"])
```

```text
13524 1.7 1.7
BValue(b=0.8397612929136761, sigma=0.008201957828422087, n=8465)
Omori(p=1.0867914401432452, c=4.83045054452574, k=3106.7737649687997, n=8465)
```

`read_catalog` returns a table with the columns `t`, `lat`, `lon`, `mw`, `dt_days` and
`dist_km`, and records in its `attrs` how many rows were read, parsed and removed by the
distance limit. `analyze_case` takes that table, the mainshock mapping and the same
settings as the `[analysis]` section, as keyword arguments. Section 5.3 lists what it
returns.

### 7.4 Read b against the threshold

The standard error is the scatter at one threshold and says nothing about how b moves as
the threshold moves. The stability curve shows that:

```python
from tremor_lab import b_stability

curve = b_stability(catalog["mw"].to_numpy(), dm=0.1)
for m, b, sigma, n in zip(curve.thresholds, curve.b, curve.sigma, curve.n):
    print(f"M >= {m:.1f}   b {b:.3f} +/- {sigma:.3f}   n {n}")
```

```text
M >= 1.5   b 0.765 +/- 0.006   n 11162
M >= 1.6   b 0.805 +/- 0.007   n 9806
M >= 1.7   b 0.840 +/- 0.008   n 8465
M >= 1.8   b 0.867 +/- 0.009   n 7178
M >= 1.9   b 0.881 +/- 0.010   n 5964
M >= 2.0   b 0.945 +/- 0.013   n 5192
M >= 2.1   b 0.987 +/- 0.015   n 4344
```

(the curve continues to M 4.0). Above a correctly estimated completeness magnitude b should
not depend on where the threshold is put. A curve that keeps climbing above the estimated
Mc means completeness was placed too low, or that the sample is not one Gutenberg-Richter
population. Report the threshold with every b-value.

### 7.5 Fit the decay above a higher threshold, or after the first day

The command line fits the decay from the mainshock. Python lets you choose the sample and
the start of the observation interval. `at_or_above` selects events at or above the lower
edge of a bin on the magnitude grid, which a plain `>=` on decimal magnitudes can get wrong
by rounding error:

```python
from tremor_lab import fit_omori
from tremor_lab.grid import at_or_above

mags = catalog["mw"].to_numpy()
times = catalog["dt_days"].to_numpy()

in_class = times[at_or_above(mags, 2.0 - 0.1 / 2)]   # the M 2.0 completeness class
fit = fit_omori(in_class)                            # observation starts at the mainshock
later = fit_omori(in_class, t_start=1.0)             # observation starts at day 1
print(fit)
print(later)
```

```text
Omori(p=1.2042484550790684, c=2.8178484985411942, k=2284.9618737148635, n=5192)
Omori(p=1.2216479877518964, c=3.0959373604082208, k=2471.332554625212, n=4643)
```

The offset c is the parameter that moves when early incompleteness is removed. Refit
at several thresholds before you interpret it.

### 7.6 Settle a borderline decay test

A simulated p-value is itself an estimate, with a standard error that shrinks as the
replicate count grows. When the verdict says `borderline`, raise the count. In a settings
file:

```toml
[analysis]
n_fit_simulations = 5000
```

In Python, either through `analyze_case(..., n_fit_simulations=5000)` or on a sample you
have chosen:

```python
from tremor_lab import omori_fit_test

test = omori_fit_test(in_class, fit, n_simulations=5000, seed=0)
print(test.statistic, test.p_value, test.p_value_se, test.method)
```

The cost is linear in the count: about seven milliseconds per replicate on 5,000 events.
Use the same seed when you want to repeat a result, and report the count.

### 7.7 Change a constant

For one call, pass the keyword argument. For a Python session, reassign the name in
`tremor_lab.constants`; every later call reads it, including calls made inside
`analyze_case`. From a settings file, use `[constants]`.

```python
from tremor_lab import bath_mag, constants, mc_maxcurvature

print(mc_maxcurvature(mags))                    # 1.7, using `mags` from section 7.5
print(mc_maxcurvature(mags, correction=0.0))    # 1.5, without the +0.2 correction
print(f"{bath_mag(7.1):.2f}")                   # 5.95, with the published deficit of 1.15
constants.DELTA_MB = 1.2                        # the whole session from here on
print(f"{bath_mag(7.1):.2f}")                   # 5.90
constants.DELTA_MB = 1.15                       # put it back
```

Three constants, `OMORI_C_FLOOR`, `OMORI_P_MIN` and `OMORI_P_MAX`, judge a fit rather than
enter it. They have no keyword argument and are changed by reassignment or from a settings
file.

### 7.8 Report the choices behind a number

The last block of the command-line report lists the window, bin width, completeness
correction, threshold and where it came from, the resample and replicate counts, the seed
and any changed constant. Paste it into a methods section or a table note. When you work
from Python, state the same choices yourself: the report block is the checklist.

### 7.9 Convert magnitudes on their own

```python
from tremor_lab import to_mw

print(to_mw([5.0, 6.5, 4.0, 3.2], ["mb", "Ms", "ml", "mww"]))
```

```text
[5.28  6.515 4.    3.2  ]
```

The `mb` value is converted, the `Ms` value is converted, and the `ml` and `mww` values
pass through unchanged.

### 7.10 Compare the two b-value estimators

The package reports the Aki estimator with Utsu's half-bin offset. The exact maximum
likelihood estimator for binned magnitudes (Tinti and Mulargia, 1987) is available beside
it, so the choice can be tested instead of assumed:

```python
from tremor_lab import b_value_aki, b_value_tinti

print(b_value_aki(mags, 1.7))
print(b_value_tinti(mags, 1.7))
```

On this catalogue the two differ by about 0.3 per cent at a bin width of 0.1. The README
explains the difference and how it shrinks with the bin width.

### 7.11 Reproduce the numbers and figures in the paper

```
python examples/regenerate_paper_figures.py
python examples/make_paper_figures.py
```

The first reads `tests/data/kahramanmaras_180d.csv`, `examples/hectormine.csv` and
`examples/ridgecrest.csv`, recomputes every number in the application section and the
supplementary tables with the package's own functions, and writes
`docs/tremor_lab_regeneration_data.json`. It takes about three minutes, nearly all of it
the calibrated decay-fit tests; `--skip-fit-test` leaves those out. The second draws
Figures 1 to 3 into `figures/` from that file and needs matplotlib. Fetching the
catalogues again returns different files (section 7.2), so regenerate from the committed
ones.

### 7.12 Check against an independent implementation

```
pip install seismostats
python examples/compare_with_seismostats.py
```

Install `seismostats` in a separate environment; it is not a dependency. The script states
which documented difference of convention accounts for each disagreement.

## 8. Messages and what they mean

Every refusal is one line on standard error, prefixed `tremor-lab:`, and exit status 2.

| Message begins | Cause and remedy |
| --- | --- |
| `<file>: missing [catalog] section` | A required section is absent, or its name is misspelt. `[catalog]` and `[mainshock]` are required. |
| `<file>: unknown section(s) [...]` | A section name the tool does not know. |
| `analyze_case() got an unexpected keyword argument 'x'` | An `[analysis]` key the tool does not know. Compare the spelling with section 4. |
| `<file>: [analysis] x must be a number` | A numeric setting was given as text or as `true` or `false`. |
| `<file>: unknown constant 'X' in [constants]` | Only the names in the appendix can be set. The message lists them. |
| `<file>: [constants] X must be a number` | A constant was given a non-number. |
| `these column names are not in the file: {...}` | A column in `[catalog.columns]` does not exist in the CSV. The message lists the file's real columns. |
| `unknown keys in [catalog]: [...]` | A misspelt key in the `[catalog]` section. |
| `unknown column roles [...]` | A role in `[catalog.columns]` the reader does not know. |
| `no magnitude column was given` or `no time column was given` | A required role is missing. |
| `radius_km cannot be applied to a catalogue configured with dt_days` | The elapsed-days shape has no positions. Remove `radius_km`, or use a raw export. |
| `mag_type cannot be applied to a catalogue configured with dt_days` | No magnitude is converted on that shape. |
| `the time column (...) mixes UTC offsets` | Put the file on one offset before reading it. |
| `[Errno 2] No such file or directory: ...` | The `path` is wrong. It is read relative to the settings file's folder. |
| `no frequency-magnitude table was written to ...` | The window held no events, so `--fmd-out` had nothing to write. |
| `the frequency-magnitude table could not be written to ...` | The output folder does not exist, or the file is open in another program. |

One message is a warning and does not stop the run: `dayfirst=True was asked for, but N of
M rows ... state the year first`. It says the flag changed fewer rows than you might have
expected, and how many.

## Appendix: the published constants

All 25 live in `tremor_lab.constants`. Each is a default, not a fixed rule of the package.
The values below are checked against the module by `tests/test_user_guide.py`.

| Name | Default | Meaning |
| --- | --- | --- |
| `ENERGY_A` | 1.5 | Slope of log10 E = a M + b, with E in joules. |
| `ENERGY_B` | 4.8 | Intercept of the same relation. |
| `DELTA_MB` | 1.15 | Bath magnitude deficit between a mainshock and its largest aftershock. |
| `MS_MW_BRANCH` | 6.1 | Ms at which the two Scordilis Ms-to-Mw branches meet. |
| `MS_MW_LOW_SLOPE` | 0.67 | Slope of the Ms-to-Mw relation below the branch. |
| `MS_MW_LOW_INTERCEPT` | 2.07 | Intercept below the branch. |
| `MS_MW_HIGH_SLOPE` | 0.99 | Slope above the branch. |
| `MS_MW_HIGH_INTERCEPT` | 0.08 | Intercept above the branch. |
| `MB_MW_SLOPE` | 0.85 | Slope of the Scordilis mb-to-Mw relation. |
| `MB_MW_INTERCEPT` | 1.03 | Intercept of the same relation. |
| `EARTH_RADIUS_KM` | 6371.0 | Sphere radius for great-circle distances. |
| `DM` | 0.1 | Magnitude bin width. |
| `MC_CORRECTION` | 0.2 | Added to the maximum-curvature completeness magnitude. |
| `SHI_BOLT_K` | 2.30 | Coefficient of the b-value standard error, the published rounding of ln 10. |
| `WINDOW_DAYS` | 180.0 | Default aftershock window in days. |
| `OMORI_C0` | 0.5 | Starting value of c in the decay search, in days. |
| `OMORI_P0` | 1.1 | Starting value of p in the decay search. |
| `GRID_TOLERANCE` | 1e-9 | How close a magnitude must come to a bin edge to count as on it. |
| `N_BOOT` | 200 | Bootstrap resamples for the uncertainty on p and c. |
| `N_FIT_SIMULATIONS` | 600 | Replicates that calibrate the decay-fit test. |
| `OMORI_C_FLOOR` | 0.001 | A fitted c below this is flagged as not identified. |
| `OMORI_P_MIN` | 0.5 | A fitted p below this is flagged. |
| `OMORI_P_MAX` | 2.0 | A fitted p above this is flagged. |
| `MIN_EVENTS_FOR_MC` | 50 | Events a window needs before Mc is estimated. |
| `MIN_EVENTS_FOR_FIT` | 100 | Events a sample needs before b and the decay are fitted. |
