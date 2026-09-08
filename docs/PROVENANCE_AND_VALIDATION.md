# Provenance, validation, and what may not be claimed

This is the canonical copy. It is the handover written for whoever drafts the
methods article, and it doubles as the project's design record: what each
estimator implements, how independent the implementations really are, what has
been checked against third-party code and what has not, the findings worth
publishing, and the list of things that must not be claimed.

It lives here rather than beside the thesis because this repository is what gets
pushed, archived and handed on; a design record outside it is a design record
that will be lost.

---

For whoever writes the software/methods article. Everything below is fact, checked
against the code and the test suite on 7 September 2026. Section 12 lists the things
that must **not** be claimed; read it before drafting.

Author: Aik Kazarian. Single-author methods paper. Companion to the PhD thesis
"Triggers: The Effect and Interaction of Earthquakes on the Example of Strong Events in
the Trans-Caucasus and Anatolian Region".

---

## 1. What the tool is, in one paragraph

Tremor Lab is a small open-source package that computes the standard estimators of
catalogue seismology for a single aftershock sequence: the magnitude of completeness,
the Gutenberg-Richter b-value, the modified Omori-Utsu decay, radiated energy and the
Bath energy screen, and magnitude homogenisation to the moment scale. It exists in two
forms that share no code: a Python package (the authority for published values) and a
single self-contained HTML page that runs the same estimators in a browser with nothing
installed. Both reproduce a set of locked reference values computed for the 2023
Kahramanmaras sequence, and the page demonstrates this to its reader on load.

## 2. The contribution, stated honestly

The estimators are standard and decades old. **The contribution is not new statistics.**
It is:

1. **Reduced friction.** The same analysis runs from Python, from one command line, or
   from a web page that needs no installation, no programming environment and no server.
2. **Reproducibility made concrete.** Every reported value is reproducible from the
   catalogue with one command, and the settings of a run are a file that can be archived
   beside the paper.
3. **Correctness demonstrated rather than asserted.** Three test layers, a mutation
   check that proves the tests can fail, a self-test that runs in front of the reader,
   and goodness-of-fit tests that are themselves calibrated rather than assumed
   (section 10) - which is how the tool came to reject its own headline decay fit.
4. **A cross-implementation check that is genuinely independent** — see section 5 and
   the caveat in section 7, and now a third-party check against `seismostats`
   (section 10a), which reproduces the completeness magnitude exactly and the
   b-value to 0.32 per cent. Read section 9b before leaning on it: the two
   in-house implementations agreed for months on a floating-point defect they
   both had, so agreement between them is evidence about the code, not proof
   of it.

Framing that works: this is infrastructure for teaching and for settings without a
computational stack, and a worked example of what "validated" should mean for a small
scientific tool. Framing that does not work: novelty of method.

## 3. What it computes, with equations and references

**Magnitude of completeness, Mc** — maximum-curvature method: the mode of the
incremental frequency-magnitude distribution, plus an additive correction, +0.2 by
default. *Wiemer and Wyss (2000).*

**Gutenberg-Richter b-value** — Aki maximum-likelihood estimator on events at or above
Mc, with the Shi and Bolt standard error:

    b     = 1 / (ln10 (mean(M) - (Mc - dM/2)))
    sigma = 2.30 b^2 sqrt( sum (M - mean(M))^2 / (n (n - 1)) )

The `- dM/2` is the lower edge of the completeness bin: a magnitude reported as Mc
stands for the interval Mc +/- dM/2. *Aki (1965); Shi and Bolt (1982).*

**Modified Omori-Utsu decay** — n(t) = k / (c + t)^p, by maximum likelihood on unbinned
post-mainshock times observed on (0, T]. *Ogata (1983); Utsu, Ogata and Matsu'ura
(1995).* Uncertainty on p and c from a nonparametric bootstrap over the occurrence
times: the elapsed times are resampled with replacement, holding both the event count and
the observation interval fixed at their observed values, and the figure reported is the
standard deviation of p and of c across the refits (200 resamples by default, seeded and
reproducible). *Efron (1979).* Two things must be said about it in the paper. It is a
standard deviation across resamples, **not a calibrated confidence interval**, and nothing
in the test suite checks its coverage. And it measures scatter conditional on the observed
event count, under the model's own premise that the times are independent draws from
k / (c + t)^p; where a sequence contains secondary bursts that premise fails and the
spread is optimistic. The reference catalogue is exactly such a sequence — it contains the
M 7.6 Elbistan event and its own aftershocks — so the quoted spread on p should be
described as indicative.

**Radiated energy** — log10 E = 1.5 M + 4.8, E in joules. *Gutenberg and Richter (1956);
Kanamori (1977).*

**Bath energy screen** — the ratio of a secondary event's energy to the energy of the
largest aftershock expected under Bath's law:
10^(1.5 (M_sec - (M_main - dMB))), with dMB = 1.15 by default. *Bath (1965).*

**Homogenisation to Mw** — Scordilis global relations. Ms in two branches
(0.67 Ms + 2.07 at or below Ms 6.1; 0.99 Ms + 0.08 above), mb linearly
(0.85 mb + 1.03). Labels beginning "mw" pass through; ml, md and unlabelled magnitudes
are taken as published, because no global relation exists for them. *Scordilis (2006).*

**Optimiser** — SciPy in the Python package, a hand-written Nelder-Mead simplex in the
browser page and in the spreadsheet version. *Nelder and Mead (1965).*

## 4. Three implementation details worth a paragraph each in the paper

**k is profiled out, not searched for.** In the Python package, only c and p are
searched. At any (c, p) the likelihood is maximised in closed form by k = n / I(c, p),
where I is the integrated rate over the observation interval. This removes one dimension
and any dependence on a starting value for k. A test asserts the profiled optimum is at
least as likely as a full three-parameter search and independent of the starting point.

**The observation interval ends at the last event, not at the nominal window.** On the
reference catalogue those differ — 179.596 days against 180 — and the choice moves p in
the third decimal. This is a real reporting decision that papers usually leave implicit.

**The fit can start after the origin.** `fit_omori` takes a `t_start`, so the hours in
which the network was still recovering can be excluded rather than absorbed into c.
Starting at zero is the default and leaves the reference values untouched. On the
reference catalogue c falls from 0.497 to 0.356 as the start moves out to one day while
p holds near 1.13 — the evidence that, of the two, c is the less identified.

## 5. Architecture, and how independent it really is

**Read this section carefully. Overstating the independence is the easiest way to make a
claim the code will not support.**

Four artefacts, written at different times, in **two numeric lineages** — one
Python/SciPy, one hand-written simplex:

| | Language | Optimiser | Lineage | Role |
|---|---|---|---|---|
| `SeismoSheet.gs` | Google Apps Script | hand-written Nelder-Mead | simplex | spreadsheet version, prior work |
| `web/tremor-lab.html` | JavaScript | hand-written Nelder-Mead, multi-start | simplex | the browser page |
| `thesis_analysis.py` | Python + SciPy | SciPy Nelder-Mead, 3 parameters | SciPy | the reference used for the thesis |
| `tremor_lab` (this) | Python + SciPy | SciPy Nelder-Mead, 2 parameters (k profiled) | SciPy | the package |

Code is carried over *within* a lineage rather than re-derived. `tremor_lab` inherits the
estimator bodies of `thesis_analysis.py` — `haversine_km`, `energy_joules`, the Bath
ratio and the Aki estimator are the same formulas with the constants promoted to
overridable arguments — and both fit Omori with
`scipy.optimize.minimize(method="Nelder-Mead")`. The browser page's numeric core is
transplanted from `SeismoSheet.gs`, not rewritten.

**There is therefore exactly one independent cross-check: between the two lineages.**
The claim that survives is: *two independently written codebases, using different
optimisers, agree to three decimals on a real sequence.* Do not write "three independent
implementations" or "four".

The page differs from the spreadsheet deliberately in five ways, all now stated in the
page's own footer. Three change the numerics: the bootstrap is seeded and therefore
reproducible (the spreadsheet used unseeded `Math.random`); times are handled in UTC, so
elapsed days do not depend on the viewer's time zone or on daylight saving; and the
simplex is restarted from several starting points, because one start was found to settle
short of the optimum on a real sequence. Two bring the page into line with the Python
package rather than with the spreadsheet: the degenerate-fit caution of section 8, and
the density gates (the spreadsheet declines below 40 events for Mc and 20 for b; the page
and the package require more than 50 and more than 100).

## 6. The locked reference values

2023 Kahramanmaras sequence, KOERI catalogue, 180-day window, threshold M >= 3.5:

| quantity | value |
|---|---|
| events in window | 3469 |
| Mc (maximum curvature) | 3.4 |
| b-value (Aki) | 0.844 +/- 0.019 on 1529 events |
| Omori | p = 1.161, c = 0.497, k = 359 |

At full precision the fit gives p 1.160897, c 0.496822, k 358.868. The observation
interval is 179.596 days.

These were produced by `thesis_analysis.py` and reproduced independently by
`SeismoSheet.gs`, whose optimiser is hand-written rather than SciPy's. The package
reproduces them, and so does the browser page — verified in a real browser, not by
inspection. `tests/test_regression.py` fails if any of them moves.

## 7. The most publishable finding: agreement on one catalogue proves nothing

This deserves its own subsection in the paper, because it is a general lesson, not a
bug report.

An adversarial audit found that the browser page and the Python package were selecting
**different sets of events** above the threshold. The page took every event at or above
`threshold - dM/2`; the package, the spreadsheet and the thesis script all take events
at or above `threshold` itself, with the half-bin correction applied later, inside the
Aki estimator where it belongs.

The two agree **only when the magnitudes and the threshold both sit on the same 0.1
grid** — which is exactly the condition of the locked reference case. Off that grid the
divergence is large: at a threshold of 4.05 on the same catalogue, one implementation
reported 625 events and b = 1.167, the other 498 events and b = 0.929. A 26% difference
in b, invisible to the reference test.

The methodological point: a cross-implementation check performed at a single, tidy
operating point can certify agreement that does not exist elsewhere. Cross-validation
has to be run across the parameter space, not at the published values.

The bug is fixed; the two now agree on off-grid thresholds, non-default bin widths, and
homogenised (therefore off-grid) magnitudes from USGS exports.

## 8. The second publishable finding: a converged fit that means nothing

On the 1988 Spitak catalogue at a threshold of M 4.2, both implementations converge, and
disagree slightly on p. Investigation showed why: the earliest event in that catalogue
sits 2.3148e-6 days — exactly 0.200 seconds — after the mainshock origin time. That drives the
offset c onto **zero, the boundary of the model**, where the likelihood has no interior
maximum: the negative log-likelihood keeps improving as c shrinks, from 392.78 at
c = 0.01 to 389.680 at c = 1e-12. p lands at ~0.31, far below anything reported for real
sequences.

Both implementations were reporting that number with no caveat. Both now mark it. The
tool reports the fit **and** states that c has collapsed onto the boundary and is not
identified by that catalogue, and that p lies outside the range compiled for real
sequences. The number is never suppressed — the caution is added.

Useful for the paper: an Omori fit can converge, report three tidy parameters, and be
meaningless. The diagnostic is cheap and almost never implemented.

## 9. The third publishable finding: completeness is uncertain by 0.7 magnitude units, and b follows it

The tool now estimates the completeness magnitude three ways. On the reference
catalogue they do not agree:

| method | Mc | b at that threshold |
|---|---|---|
| maximum curvature (Wiemer and Wyss 2000) | 3.4 | 0.851 +/- 0.018 |
| goodness of fit, 90 per cent (Wiemer and Wyss 2000) | 3.0 | 0.743 +/- 0.011 |
| b-value stability (Cao and Gao 2002; Woessner and Wiemer 2005) | 4.1 | 1.041 +/- 0.046 |

The published analysis uses a threshold of 3.5, which is none of these, and gives
b = 0.844 +/- 0.019.

The b-against-threshold curve shows why they disagree. b does not settle:

| threshold | 3.0 | 3.2 | 3.4 | 3.6 | 3.8 | 4.0 | 4.2 | 4.4 |
|---|---|---|---|---|---|---|---|---|
| b | 0.743 | 0.829 | 0.851 | 0.878 | 0.985 | 1.028 | 1.048 | 1.124 |

Above a correctly estimated completeness magnitude b should be flat. This one
climbs by 0.381 across the range, from 0.743 to 1.124, which is twenty times the
+/- 0.019 quoted at the published threshold. **The Shi and Bolt standard error describes
sampling scatter at one chosen threshold and nothing else.** Reporting
b = 0.844 +/- 0.019 without the curve overstates the precision by more than an
order of magnitude.

Two things must therefore appear in the paper.

1. **The stability curve as a figure**, with the three completeness estimates
   marked. The tool draws it; take the figure from the tool.
2. **A stated reason for the threshold of 3.5**, since the tool's own maximum
   curvature estimate is 3.4 and the manuscript is currently silent on the
   difference.

An honest formulation: *"b = 0.844 +/- 0.019 at a threshold of M 3.5, where the
error is the Shi and Bolt sampling uncertainty. Completeness estimates for this
catalogue range from M 3.0 to M 4.1 depending on method, over which range b
varies from 0.74 to 1.04; the choice of completeness magnitude dominates the
uncertainty in b."*

What this is **not**: evidence that the code is wrong. The estimator is correct
and the reference values stand. It is a property of the catalogue, and the tool
now measures it instead of hiding it.

## 9a. A fourth finding: the reference sequence is not a single Omori decay

With the fit test calibrated (see section 10), the 180-day Kahramanmaras window
at threshold M 3.5 gives KS 0.0210 at **p = 0.036 +/- 0.003, which rejects a
single modified Omori-Utsu decay at the 5 per cent level.**

That figure is from 5,000 replicates, and quote it with its Monte Carlo error,
because the p-value is an estimate and this one sits close to the threshold.
Three independent seeds at 5,000 gave 0.0372, 0.0358 and 0.0354. The shipped
default is 600 replicates, where the standard error is about 0.008 and the tool
reports the same case as *borderline* rather than choosing - which is the correct
caution at that precision, not a different answer. For the paper, run it at 5,000:
on the command line put `n_fit_simulations = 5000` under `[analysis]`, and in the
browser raise the "fit test replicates" box.

This is the physically expected answer and should be presented as a result, not
an embarrassment. The window contains the M 7.6 Elbistan event nine hours after
the mainshock, with an aftershock sequence of its own. One decay curve is not
the right model for two overlapping sequences, and the residual test says so.

Three consequences for the paper.

1. **p = 1.161 and c = 0.497 describe a mixture**, not a single sequence. They
   remain the correct maximum-likelihood estimates for the model as specified,
   and the locked reference values stand; but the model is rejected by its own
   residual test, and the paper must say so rather than leave a reader to find it.
2. **The bootstrap spread on p and c is optimistic**, because resampling assumes
   the times are independent draws from one decay, which the same test rejects.
3. **This is the natural argument for the declustering the tool does not do**,
   and for the follow-on work: fitting the two sequences separately, or an ETAS
   model, is the obvious next step and the paper can say so from evidence rather
   than from convention.

Note for honesty: an earlier version of this handover reported p = 0.506 and
called the fit adequate. That figure came from the uncalibrated test and is
wrong. Do not use it. A later version reported "about 0.03" from a 200-replicate
run whose own standard error was 0.011; the number was in the right place but
quoted with more confidence than 200 replicates support. Use 0.036 +/- 0.003.

## 9b. A fifth finding, and a methodological one: reported magnitudes are
decimal, and floating-point arithmetic is not

This one is about software rather than about the Earth, and it is the kind of
thing a methods paper is the right place to say.

Every estimator here compares a magnitude against a bin edge. Is this event at or
above the completeness threshold; which bin does it fall in. Magnitudes are
reported on a decimal grid - 3.4, 3.5, 4.1 - and none of those numbers exists
exactly in binary floating point, nor does an edge computed from them. Written
the obvious way, `m >= mc - dm / 2`, those comparisons are settled by rounding
error in the sixteenth decimal place.

At the default bin width of 0.1 the errors fell harmlessly on the reference
catalogue, which is why three earlier audits, a third-party cross-check and 210
tests all passed over it. At 0.2 - a documented setting and an obvious
sensitivity check for a referee - they did not:

| at dM 0.2 | before | after |
|---|---|---|
| events at or above Mc 4.2 | 394 | 498 |
| b at Mc 4.2 | 0.936 | 1.182 |
| b-stability curve, thresholds 4.0 to 4.6 | 1.143, 0.936, 0.995, 1.309 | 1.143, 1.182, 1.192, 1.309 |
| bins, for one event at every 0.1 step | 1, 3, 1, 2, 3, 2, 2, 1, 2 | 2 throughout |
| goodness-of-fit R at candidate 3.4 | 84.5 | 91.7 |

`4.2 - 0.2 / 2` evaluates to 4.1000000000000005, so every event reported at
exactly M 4.1 failed its own completeness threshold and an entire magnitude class
was dropped from the sample. b was wrong by 0.25, six times the standard error
quoted beside it. The b-stability curve - the diagnostic whose whole purpose is to
reveal a mis-set completeness magnitude - acquired a dip that was an artefact of
arithmetic, and the goodness-of-fit statistic was scored seven points low against
a threshold of 90.

Three things to say about it in the paper.

1. **The published values are unaffected.** They are on the dM 0.1 path, and every
   locked number - 3,469 events, Mc 3.4, b 0.844 +/- 0.019 on 1,529, p 1.161,
   c 0.497, k 358.9 - is unchanged, before and after. This is a defect in the
   sensitivity analysis a referee would ask for, not in the headline result.
2. **Agreement between implementations did not catch it, because both had it.**
   The Python package and the JavaScript page made the same mistake for the same
   reason, so they agreed with each other while both being wrong. That is worth
   stating plainly: cross-implementation agreement is evidence, but it is evidence
   of shared assumptions as much as of correctness, and it is blind to anything
   inherited from the arithmetic both sit on.
3. **The fix is one tolerance, applied everywhere a magnitude meets an edge.**
   `constants.GRID_TOLERANCE`, 1e-9, is a published default like any other. Any
   two magnitudes a real catalogue distinguishes differ by at least 0.001, six
   orders of magnitude more, so it cannot merge two genuinely different values.

Fixing it exposed a second divergence of the same kind. `span / dm` is exactly 2.5
at dM 0.2, and Python's `round` sends a half down while JavaScript's `Math.round`
sends it up, so the two implementations averaged the b-stability curve over
different windows and reported different completeness magnitudes - 4.0 against 4.2
- from curves that were identical row for row. Both now round halves up, as the
binning does. The two implementations agree at every bin width tested: 3.2 at
dM 0.05, 4.1 at 0.1, 4.2 at 0.2, 4.0 at 0.25, 4.5 at 0.5.

A sentence that works: *"Magnitude comparisons are performed against the decimal
grid on which magnitudes are reported rather than in binary floating point; at bin
widths where a completeness threshold coincides with a reported magnitude, the
naive comparison excludes that magnitude class entirely and biases b by several
times its standard error."*

## 10. How correctness is demonstrated

**Layer 1, analytic.** Expected values are exact arithmetic on a published relation: the
energy exponent at Mw 7.8, the Bath ratio at a one-magnitude anomaly, both Scordilis
conversions, one degree of latitude in km.

**Layer 2, synthetic recovery.** Magnitudes drawn from a Gutenberg-Richter law with known
b, complete from the lower edge of the Mc bin and reported on the 0.1 grid — which is
what the half-bin correction assumes; a naive construction sits 10% low, and that
mistake is documented in the test. Occurrence times drawn from an Omori process with
known p, c and k by inverse-CDF sampling. Both recovered across several seeds.

**Layer 3, regression.** The locked values above, on the bundled catalogue.

**Model adequacy.** The decay is checked by the Ogata transformed-time
Kolmogorov-Smirnov residual test, with one subtlety that must be stated in the
paper. The parameters are estimated from the very times being tested, so the
fitted curve hugs the data and the textbook Kolmogorov p-value does not apply.
Measured on 300 sequences drawn from the model, the uncalibrated form rejected at
the 5 per cent level in 0 per cent of cases, with a mean p-value of 0.87: it
could hardly ever fail. The tool therefore simulates the null distribution by
parametric bootstrap, refitting each replicate; that rejects at 7.5 per cent
against a 5 per cent target, with a mean p-value of 0.473.

The calibration itself reproduces across the two implementations, which is worth
stating because it is the part a reader is most likely to doubt. On the reference
catalogue, 600 replicates in Python and 600 in the browser's JavaScript - different
random number generators, different optimisers - give a mean simulated KS statistic
of 0.01445 against 0.01440, a 97.5th percentile of 0.02170 against 0.02181, and the
same p-value of 0.0416. What is reproduced is not only the point estimate but the
simulated null distribution behind its p-value.

Because the p-value is an estimate, the tool reports its binomial standard error
beside it and refuses a verdict when the two are within two standard errors of
0.05. It also refuses a verdict outright on the uncalibrated asymptotic p-value,
which on this catalogue reads 0.508 where the calibrated one reads 0.036: the
number is shown, since hiding it would be worse, but nothing may be concluded
from it.

The result on the reference catalogue changes accordingly, and this is section
9a below. The test has a negative control in the suite: a rate that rises with
time cannot be Omori and is rejected, and so is a second sequence injected into
the window. A uniform sequence is deliberately not used as the control, because
the model covers it at p = 0, where the rate is constant.

**Mutation testing.** The suite was checked by deliberately breaking the estimators ten
ways and confirming a test fails each time: deleting the cos(latitude) factor from the
haversine; dropping the half-bin offset; changing the Shi and Bolt coefficient from 2.30
to 2.0; setting the Mc correction to zero; fixing the observation interval at 180 days;
swapping the two Scordilis branches; changing the Earth radius by 1%; swapping p and c in
the returned result; and two plumbing mutations in the command-line tool and the
bootstrap seed. **All ten are caught.** Three of them passed silently before the audit,
against a suite that was then 136 tests and green.

This is the part most worth writing up: a green test suite is not evidence until you have
shown it can go red.

## 10a. The independent cross-check, and the estimator choice it exposed

`examples/compare_with_seismostats.py` runs the reference catalogue through
`seismostats`, an independent package from the Swiss Seismological Service at ETH
Zurich, written by different people from the same literature. This is the
third-party comparison a referee will ask for, and the result is good:

| quantity | Tremor Lab | seismostats | |
|---|---|---|---|
| maximum-curvature Mc | 3.4 | 3.4 | identical |
| events above M 3.5 | 1529 | 1529 | identical |
| b-value | 0.844132 | 0.846804 | 0.32 per cent apart |
| completeness by b-stability | 4.1 | 4.0 | one bin apart |

**The b-value gap is not an error in either package.** They implement two
different published estimators:

- Tremor Lab reports **Aki (1965) with Utsu's half-bin offset**,
  b = 1 / (ln10 (mean(M) - (Mc - dM/2))). This is what the thesis and the
  spreadsheet use, and it is what the locked reference value 0.844 is.
- `seismostats` reports **Tinti and Mulargia (1987)**,
  b = ln(1 + dM / mean(M - Mc)) / (dM ln10), the exact maximum-likelihood
  solution for magnitudes reported on a grid.

The half-bin form is the first-order approximation of the exact solution, and the
two converge as the bin narrows: 0.32 per cent apart at dM 0.1, 0.004 per cent at
dM 0.01. The package now implements both; `b_value_tinti` reproduces
`seismostats` to six decimals, which is the strongest single piece of correctness
evidence the project has.

Two honest caveats to carry into the paper.

1. **The comparison covers completeness and the b-value only.** `seismostats`
   does not implement the modified Omori-Utsu fit, the energy relation or the
   Bath screen, so none of those has been checked against anything third-party.
   Say what was compared; do not let it stand for the whole tool.
2. **Standard errors were not compared directly.** Both use Shi and Bolt (1982),
   but `seismostats` takes ln(10) as the coefficient where Tremor Lab takes the
   published rounding of 2.30. That accounts for the remainder of the sigma
   difference and is an adjustable constant, not a defect.

A sentence that works: *"The completeness magnitude and event count reproduce
exactly against an independent implementation; the b-value agrees to 0.32 per
cent, the residual being the documented difference between the Aki estimator with
Utsu's half-bin correction and the exact binned maximum-likelihood solution of
Tinti and Mulargia, which converge as the magnitude bin narrows."*

## 11. Facts and figures

- Package: `tremor-lab` 1.0.0, MIT licence, Python 3.11+, tagged `v1.0.0`.
- Runtime dependencies: NumPy (>=2.0,<3), SciPy (>=1.13,<2), pandas (>=2.2,<3). Nothing else.
- Source: about 2,600 lines across 10 modules. **284 tests**.
- Browser page: one file, 2,660 lines, 196 KB, pure ASCII, no external requests except
  Google Fonts. Runs offline from a double-click. Eight settings are editable on it:
  the window, the threshold, the bin width, the Mc correction, the Bath deficit, the
  bootstrap count, the fit-test replicate count and the seed. The rest of the published
  constants are fixed in the page and adjustable in the package, which is the authority
  for published values in any case.
- Public API: 30 names, including `b_stability`, `mc_b_stability`,
  `mc_goodness_of_fit`, `omori_fit_test` and `b_value_tinti`. 25 published
  constants, of which three - the two Omori plausibility bounds and the offset
  floor - judge a fit rather than entering one and so have no keyword argument. Twenty of them are overridable three
  ways — per call by keyword, per session by reassignment, or from a settings file; the
  three fit-quality bounds (`OMORI_C_FLOOR`, `OMORI_P_MIN`, `OMORI_P_MAX`) have no
  keyword argument and are overridable the latter two ways only.
- **Speed: do not quote a figure from this handover.** The tool analyses the
  3,469-event reference catalogue and draws all three charts in well under a second in a
  desktop browser, and the self-test costs more than one analysis rather than less,
  because it repeats the reference fit and adds a second Omori fit to 4,000 synthetic
  times. Neither implementation carries timing instrumentation, so any number in the
  paper must be measured on the machine the paper names, with the method stated.
- The package declines to estimate at 50 events or fewer in the window (no Mc) and at
  100 events or fewer above threshold (no b-value, no decay fit) — the gates are strictly
  "more than" — reporting the reason instead of forcing a fit.

## 12. What must NOT be claimed

Read this before drafting. Each of these is a real trap.

1. **Do not present the Van, Spitak or Racha numbers as results.** During testing those
   catalogues were run with mainshock times, epicentres and windows *chosen for the test*,
   not taken from the thesis. They demonstrate that the software reads real KOERI and
   USGS exports and that the two implementations agree — nothing more. Any published
   value for those sequences must be recomputed with the thesis's own windowing rule.
2. **Do not claim the tool validates Bath's law**, or any other law. It computes a
   comparison; a single sequence neither confirms nor refutes.
3. **Do not describe the browser page as the authority.** It runs the JavaScript
   implementation. The Python package is the authority for published values. The two
   agree on the reference catalogue and on the real exports tested; say that, not more.
4. **Do not quote b = 0.844 +/- 0.019 as though that error bar were the
   uncertainty in b.** It is the sampling scatter at one threshold. The
   completeness magnitude is uncertain from M 3.0 to M 4.1 by method, and b
   varies from 0.74 to 1.04 across that range. See section 9.
5. **Do not quote k without a caveat, and do not quote a 43,000-event figure** — no such
   case exists in the repository. Of the three decay parameters k is much the weakest: it
   is derived from the fitted c and p rather than searched for. The only k-recovery test
   in the repository (`tests/test_omori.py::test_the_productivity_constant_is_recovered`)
   draws 8,678 synthetic events from p = 1.15, c = 0.5, k = 2000 and asserts only that k
   returns within 10%. At that test's own seed k comes back 1.3% low; across seeds 0-4 the
   k error reaches 6.4% while p stays within 1.4%. p is the robust number.
6. **Do not claim a DOI.** None exists yet, and the software availability statement
   cannot be written without one. The repository is public at
   github.com/HaikKaz/tremor-lab and its CI has run and passed, so both of those
   may be stated. See section 12.
7. **Do not describe the probabilistic anomaly method** (spatial rarity, the simulation).
   It is deliberately outside this tool and is a separate layer to be added once its
   specification is frozen.
8. **Do not say the estimators are novel.** See section 2.
9. **Do state the third-party comparison, and state its limits.** The estimators
   have now been checked against `seismostats` (Swiss Seismological Service, ETH
   Zurich): see section 10a. What has been compared is the completeness magnitude,
   the b-value and the b-stability criterion. The Omori-Utsu decay, the energy
   relations and the Bath screen have **not** been compared with anything
   third-party, because `seismostats` does not implement them. Do not let the
   b-value comparison imply the decay fit was independently checked.
10. **Do not let the self-test stand for validation of a reader's own analysis.** It
   checks the estimators against known values on a bundled catalogue. It says nothing
   about whether the reader chose a sensible window, threshold or mainshock.
11. **Do not present "284 tests" as coverage.** It is a count, not a measure. What can
    honestly be said is stronger and more specific: ten deliberate breakages of the
    estimators were each caught by at least one test (section 10).
12. **Do not describe the bundled fixture as raw data.** `kahramanmaras_180d.csv` is a
    derived file: it carries elapsed days and magnitudes only, already windowed to 180
    days, with no timestamps. Its provenance from the original KOERI export should be
    stated in the paper, and the export itself archived with the release.
13. **Do not claim the Scordilis relations are range-checked.** They are applied without
   bounds, matching the spreadsheet implementation, and the two Ms branches are mildly
   discontinuous across the uncalibrated 6.1-6.2 gap (6.157 against 6.119 at Ms 6.1).
   This is a deliberate fidelity choice and is documented.

## 13. Status and what is outstanding

Done: package, tests, command-line tool, browser page, README, MIT licence,
`CITATION.cff`, local git repository tagged v1.0.0, published as a private web page.

Outstanding before the paper can cite the software:

- **A public repository.** Done: github.com/HaikKaz/tremor-lab, 38 commits.
- **A Zenodo DOI.** None yet. Required for the software availability statement.
- **ORCID and affiliation** are filled in: ORCID 0009-0007-1842-1590, National
  Academy of Sciences of Armenia.

**Name.** The author's name is **Aik Kazarian**, as on the passport and on the
ORCID record, and that is what `CITATION.cff` and the Zenodo metadata carry. Earlier
drafts of the thesis and of this handover used "Haik"; use "Aik" in the paper, the
author list and the availability statement, so the byline, the ORCID record and the
software archive all match. A single author whose name is spelled two ways across a
paper and its cited software is exactly the sort of thing that splits a citation
record in two.
- **Continuous integration** runs and passes. `.github/workflows/tests.yml` lints,
  runs the suite on Python 3.11, 3.12 and 3.13, reproduces the reference values both
  directly and through the full catalogue pipeline, and runs the `seismostats`
  cross-check in a separate job. Its first execution, on commit 5c123d3, passed all
  four jobs in 63 seconds. This is worth a sentence in the paper: it means the
  locked values reproduce on a clean checkout with freshly resolved dependencies,
  on three Python versions, and not only on the author's machine.

Until the DOI exists, the availability statement cannot be written. Draft around it, or
leave a marked placeholder.

## 14. Limitations to state in the paper

The estimators inherit their known behaviour. Maximum-curvature Mc is sensitive to
binning and to short-term aftershock incompleteness in the hours after a mainshock. The
Aki b-value assumes completeness above the chosen threshold. The Omori fit needs a dense
catalogue: where the catalogue above threshold is too thin the tool reports the absence
rather than producing a fit, and two of the four case-study catalogues (Spitak, Racha)
fall in that category — which is itself a finding, already reported in the thesis.
Magnitude homogenisation is applied only where a scale column exists; where it does not,
magnitudes are used as published and the tool says so.

## 15. References to verify before submission

Aki (1965); Bath (1965); Efron (1979); Gutenberg and Richter (1944, 1956);
Kanamori (1977); Nelder and Mead (1965); Ogata (1983); Scordilis (2006);
Shi and Bolt (1982); Utsu, Ogata and Matsu'ura (1995); Wiemer and Wyss (2000).

Every primary source must be checked against the original before citation; the list
above records what the code implements, not a verified bibliography.

---

## Files the writer may want

- `tremor_lab/README.md` — the fullest existing prose description; much of the software
  section can be adapted from it.
- `tremor_lab/web/tremor-lab.html` — the browser tool; open it and press the example
  button to see the whole thing work.
- `tremor_lab/examples/reproduce_reference.py` — one command, prints the eight reference
  values beside the computed ones.
- `tremor_lab/tests/test_regression.py` — the test that locks the reference values.
- `SeismoSheet_Description.md` — the earlier spreadsheet write-up; the framing in its
  sections 1 and 6 is reusable.
