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

**Drafting the paper? Read section 0 first.** It states the settled position in one
place. The sections after it are chronological and a few of them correct each other,
so reading straight through will hand you superseded conclusions.

Author: Aik Kazarian. Single-author methods paper. Companion to the PhD thesis
"Triggers: The Effect and Interaction of Earthquakes on the Example of Strong Events in
the Trans-Caucasus and Anatolian Region".

---

## 0. Drafting brief: read this section first, and possibly only this

**If you are drafting the manuscript, this section is your brief.** Sections 1 to 15
below are the evidence base and the project's design record. They were written across
several working sessions in chronological order, so they read as an investigation rather
than as a position: section 9a is corrected by 9d, and part of 9c is corrected by 10c.
**Where two sections disagree, the later one wins**, and this brief already states the
settled position. Go into the detail sections for numbers, quotations and conditions, not
for the argument.

Author: **Aik Kazarian** (this exact spelling, matching the passport, the ORCID record
and the software archive; earlier thesis drafts say "Haik"). ORCID 0009-0007-1842-1590,
National Academy of Sciences of Armenia. Single author.

### 0.1 The claim the paper defends

A paper should defend one claim, or at most two. This work has produced eight candidate
findings and the temptation is to report them all. Do not. The recommendation, which is
mine and which you may overturn if you disagree, is this.

**Primary claim.** *A reported b-value is not a reproducible statement about a sequence
unless the magnitude threshold is reported with it, because b varies with that threshold
by roughly an order of magnitude more than its conventional uncertainty; and the standard
maximum-curvature estimator selects a threshold below the range in which b is stable.*

That claim is quantitative, it is falsifiable, it is supported by three catalogues
including a control designed to break it, and it has been independently confirmed by
published work using different data and different software. It also has an immediate
practical consequence, which is what makes it worth a methods paper: **report the
stability curve, not a single b.**

**Secondary claim, illustrative.** *In an earthquake doublet, catalogue completeness
degrades twice, once per mainshock, and the second degradation is more severe than the
first.* Measured on Kahramanmaras: Mc rises to 3.8 after the M 7.8, is recovering by nine
hours, then collapses to 4.4 when the M 7.6 Elbistan event arrives at 0.380 days. This is
a specific observation about doublets that a single mainshock cannot produce, and it
motivates the primary claim by showing how far completeness can move inside one window.

**Tremor Lab is the instrument, not the subject.** The software is what makes the
analysis reproducible and is cited and archived, but a paper whose contribution is "we
wrote a tool" is weaker than one whose contribution is a measured result the tool
establishes. Lead with the finding. Describe the tool in Data and Methods, and let the
availability statement carry the rest.

### 0.2 Where to send it

**Seismological Research Letters (SRL)** is the natural home: it publishes methods and
software contributions, it published Gulia et al. (2020) which this work cites and
partly replicates, and the length suits. **BSSA** is the alternative if the framing
leans further toward the measurement than the method. Both are SSA journals using
author-date citation and requiring a **Data and Resources** section in place of a generic
data-availability statement. Section 0.10 drafts that section.

### 0.3 The evidence, mapped to IMRaD

**Introduction.** The gap is concrete and should be stated as such: completeness
magnitude is routinely estimated by maximum curvature and then used as a fixed threshold
for a b-value, and the b-value is then reported with a Shi and Bolt (1982) standard
error that describes sampling scatter at that threshold only. Nothing in that pipeline
propagates the uncertainty in the threshold itself. The literature on time-varying
completeness is well established and should be cited rather than rediscovered: Kagan
(2004), Helmstetter et al. (2006), Gulia et al. (2020), Tan (2025). What is not
established, and what this paper supplies, is a quantitative statement of what that
omission costs.

**Data and Methods.** Three catalogues, chosen so that one is a control (see 0.4). State
for each: source, access date, magnitude type, homogenisation, spatial and temporal
selection, and the completeness treatment. Magnitude type matters and reviewers check it:
the Kahramanmaras catalogue is homogenised to Mw via Scordilis (2006), the two USGS
catalogues carry mixed types homogenised the same way, and the published comparisons in
section 10c use ML in two cases, which is why several of them are not directly
comparable. Name every estimator with its source (0.9) and name the software with its
version and DOI.

**Results.** Report the stability curves and the completeness bands without interpreting
them. The three tables in 0.4 are the core of this section. Report KS statistics with
sample sizes, not verdicts alone.

**Discussion.** Three things. First, the control: the effect appears in a single
mainshock with no unusual features, so it is a property of the estimator and not of
doublets. Second, the independent confirmations in 10c, which are strong and should be
given room. Third, and this is the part that requires courage, the correction in 0.6.

**Conclusions.** Short. The recommendation to report stability curves, and the honest
limits of what three catalogues establish.

### 0.4 The numbers, with their conditions

All b-values are Aki (1965) maximum likelihood with the Utsu (1965) half-bin offset,
uncertainties by Shi and Bolt (1982), bin width 0.1. All Mc by maximum curvature (Wiemer
and Wyss, 2000) with the +0.2 correction unless stated. b-stability is Cao and Gao (2002)
as formulated by Woessner and Wiemer (2005).

**Catalogue A, Kahramanmaras 2023.** KOERI, homogenised to Mw, 3,469 events M 3.0 to 7.6,
180 days from the M 7.8 of 6 February 2023. Note the catalogue is truncated at M 3.0,
which matters (0.7).

| threshold | 3.2 | 3.4 | 3.5 | 3.8 | 4.0 | 4.1 | 4.3 | 4.4 | 4.5 |
|---|---|---|---|---|---|---|---|---|---|
| b | 0.829 | 0.851 | **0.844** | 0.985 | 1.028 | **1.041** | 1.048 | 1.124 | 1.137 |
| sigma | 0.014 | 0.018 | **0.019** | 0.030 | 0.040 | **0.046** | 0.058 | 0.073 | 0.087 |
| n | 2697 | 1874 | **1529** | 954 | 625 | 498 | 309 | 258 | 201 |

Mc by three methods: goodness-of-fit **3.0**, maximum curvature **3.4**, b-stability
**4.1**. The thesis analyses at M 3.5, where b is 0.844 and still climbing.

**Catalogue B, Hector Mine 1999. This is the control.** USGS ComCat, 13,525 events
M >= 1.0 within 100 km, 180 days from the M 7.1 of 16 October 1999. A single mainshock
with no M 6 or larger event in the following fortnight. Landers 1992 was considered and
rejected as a control because its M 6.3 Big Bear event three hours later makes it a
doublet as well.

| threshold | 1.5 | 1.7 | 1.9 | 2.1 | 2.2 | 2.4 | 2.6 | 2.8 | 3.0 |
|---|---|---|---|---|---|---|---|---|---|
| b | 0.764 | **0.839** | 0.879 | **0.985** | 1.003 | 0.990 | 0.968 | 0.980 | 0.989 |
| n | 11163 | 8466 | 5965 | 4345 | 3526 | 2217 | 1384 | 900 | 585 |

Maximum curvature returns **1.7**; b-stability returns **2.1**. **The plateau is the
single most important result in the paper**: b sits at 0.98 plus or minus 0.02 across
nine consecutive thresholds spanning 0.9 magnitude units on thousands of events, which
proves the climb below it is a real bias and not an artefact of shrinking samples. The
estimator underestimates the threshold by 0.4 magnitude units, and b at its answer is
0.146 low against a quoted standard error of 0.008, a factor of eighteen.

**Catalogue C, Ridgecrest 2019.** USGS ComCat, 28,963 events M >= 1.0 within 100 km, 180
days from the M 7.1 of 6 July 2019. A doublet, the M 7.1 preceded by an M 6.4 about
thirty-four hours earlier.

| threshold | 1.1 | 1.3 | 1.8 | 2.2 | 2.6 | 3.0 | 3.5 |
|---|---|---|---|---|---|---|---|
| b | 0.695 | **0.733** | 0.794 | 0.820 | 0.865 | 0.924 | 1.341 |

Maximum curvature returns 1.3. **b-stability returns no answer at all**: b never
stabilises anywhere in the available range, and the implementation returns None rather
than inventing a threshold. That refusal is worth a sentence, because it is the correct
behaviour when the method's premise fails.

**Completeness in time, Kahramanmaras.** The secondary claim.

| band | 0 to 0.20 d | 0.20 to 0.38 d | **0.38 to 0.60 d** | 0.60 to 1.0 d | 1 to 2 d | 2 to 5 d | 30 to 180 d |
|---|---|---|---|---|---|---|---|
| Mc | 3.8 | 3.7 | **4.4** | 3.6 | 3.4 | 3.2 | 3.3 |
| n | 155 | 106 | 120 | 214 | 347 | 497 | 968 |

The M 7.6 Elbistan event falls at 0.380 days, inside the band with the worst
completeness. The published M 3.5 threshold is below Mc for roughly the first day, during
which 425 events enter the sample.

**Time and threshold are separable at Ridgecrest and not at Kahramanmaras.** At
Ridgecrest, dropping the first day raises b by about 0.15 at every threshold, and raising
the threshold raises it by about 0.15 at every start time, roughly additively, which is
why it never plateaus. At Kahramanmaras the two cuts remove overlapping event sets and
only the threshold effect survives. Full grids are in sections 9e and 10c.

**Decay fits.** Refitting above the early completeness magnitude drops the Omori offset c
monotonically in both catalogues, confirming that c absorbs the missing early population:
Kahramanmaras 0.497 to 0.195 days across thresholds 3.5 to 4.4, Hector Mine 4.814 to
0.710 across 1.7 to 2.6. **It does not rescue the fit.** The Kolmogorov-Smirnov departure
grows and the rejection mostly stands. Report this as a half-confirmed prediction, which
is what it is. Note also that p drifts upward with threshold, by 0.03 and 0.19
respectively, so p needs its threshold stated for the same reason b does.

### 0.5 The independent confirmations, which are strong

These are verified first-hand, with the source located and the text quoted. Full detail,
conditions and verbatim quotations are in section 10c.

**The threshold finding is independently confirmed for this same region.** Hainzl,
Kumazawa and Ogata (2024), on AFAD background seismicity from 2000 to 2022, using the
same Aki estimator, report that "the b estimate becomes stable within its uncertainties
for M c > 3.5 and scatters around 1.07". Section 0.4 finds stability at Mc 4.1 and a
value of 1.04 from the aftershock catalogue. Different data, twenty-three years, same two
conclusions. **Cite this rather than presenting the finding as new.** What is new here is
the control.

**The size of the bias is independently measured.** The same paper fits a standard ETAS
model and an incompleteness-corrected ETASI model to one catalogue and gets b = 0.56
against b = 0.87. A bias of 0.31 on identical data, from a completely different method.

**A like-for-like agreement at the same threshold.** Ali and Abdelrahman (2024) reach
Mc 4.4 for this sequence on an independent IRIS catalogue of 471 events using ZMAP. At
that threshold their c = 0.204 plus or minus 0.058 against our 0.195, with b and p within
about one standard error. This is a stronger external check than the `seismostats`
comparison in section 10a, because that one shared our input data and this one does not.
Do not tabulate k: it scales with the number of events above threshold and their
catalogue is not ours.

**Time-varying completeness is confirmed but the numbers are not comparable.** Tan (2025)
reports Mc falling from about 3.0 on the first day to about 1.5 by mid-March on an AFAD
catalogue of 50,085 events. Same shape as ours, offset by nearly a magnitude unit,
because their catalogue is complete far below M 3.0 and ours is truncated at it. Cite for
the phenomenon and the recovery timescale, not for the values.

### 0.6 The correction you must not quietly drop

Sections 9c and 9e report b estimated band by band at each band's own maximum-curvature
Mc, giving b rising with time: 0.737 to 1.018 for Kahramanmaras, 0.881 to 1.004 for
Ridgecrest. **The published time-resolved estimates run the other way.** Hainzl, Kumazawa
and Ogata (2024), using the Ogata-Katsura (1993) estimator which fits time-varying b
jointly with a time-varying detection function, report b of about 1.2 for early
aftershocks decaying to about 0.85 later, a 50 per cent coseismic increase. That is also
the premise of the Gulia and Wiemer Foreshock Traffic Light System.

Our trend has the wrong sign, and the reason is the paper's own thesis: the early bands
have the worst completeness, so they carry the largest downward bias. The Kahramanmaras
early figure rests on 72 events with an Mc estimated from those same 72 events.

**Report this.** It is better evidence for the primary claim than agreement would have
been, because the tool reproduces the artefact and the diagnosis predicts it. A reviewer
who finds this independently and sees it unacknowledged will reject the paper. A reviewer
who sees it reported and explained will trust everything else.

### 0.7 Limitations to state, not bury

- **The Kahramanmaras catalogue is truncated at M 3.0**, and the incremental counts in
  the four lowest bins are 393, 379, 427 and 396, a plateau rather than a peak. Maximum
  curvature returns 3.4 from four bins sitting immediately above a truncation floor. The
  estimate is poorly constrained and the paper should say so. If the original KOERI export
  extends below M 3.0, rerunning on it would settle the question rather than bound it.
- **Three catalogues is three catalogues.** Two doublets and one single mainshock, all
  M 7.1 to 7.8, two of them Californian. The claim generalises to the extent that a
  reviewer believes three cases generalise. Say so, and resist the urge to write "in
  general".
- **b stabilising near 1.0 in two catalogues does not mean b is 1.0.** It means the
  estimate stops depending on the threshold there. Whether the catalogue is complete at
  that threshold is a separate question this data cannot answer.
- Everything in section 14, which lists the estimators' inherited behaviour.

### 0.8 What must not be claimed

Section 12 is a list of sixteen specific traps and should be read in full before drafting.
The four that matter most here:

1. **Never report a b-value or an Omori p without its threshold** (item 14).
2. **Do not present the band-by-band b values as a measurement of how b evolved**
   (item 15, and 0.6 above).
3. **Do not claim the decay-fit rejection is caused by the doublet structure** (item 7).
   The control rejects just as hard. If the paper wants a doublet claim, cite Tan (2025)
   and Rodriguez-Perez and Zuniga (2025), who fit the two ruptures separately and find
   they differ, rather than resting it on the residual test.
4. **Do not quote the Shi and Bolt error as the uncertainty in b** (item 4). That is the
   entire point of the paper.

### 0.9 Citations

**Verified first-hand and ready to cite**, with full details and quotations in sections
10b and 10c: Ali and Abdelrahman (2024) doi:10.3390/fractalfract8050252; Hainzl, Kumazawa
and Ogata (2024) doi:10.1093/gji/ggae006; Tan (2025) doi:10.1007/s11600-024-01419-y;
Rodriguez-Perez and Zuniga (2025) doi:10.1007/s11600-024-01428-x; Convertito, Tramelli
and Godano (2024) doi:10.1038/s41598-023-50837-3; Gulia, Wiemer and Vannucci (2020)
doi:10.1785/0220190307; Hainzl (2022) doi:10.1785/0120210146; Mizrahi, Nandan and Wiemer
(2021) doi:10.1029/2021JB022379.

**Method references** for Data and Methods: Aki (1965); Utsu (1965); Shi and Bolt (1982);
Wiemer and Wyss (2000); Cao and Gao (2002); Woessner and Wiemer (2005); Tinti and
Mulargia (1987); Ogata (1988); Ogata and Katsura (1993); Scordilis (2006); Bath (1965).

**Checked by search agent, not read first-hand.** Treat as indicative and verify before
submission: Herrmann and Marzocchi (2021) SRL 92(2); Huang, Tang and Feng (2022); Kagan
(2004) doi:10.1785/012003098; Helmstetter, Kagan and Jackson (2006)
doi:10.1785/0120050067; Lippiello et al. (2019) doi:10.3390/geosciences9080355. The
Helmstetter relation mc(t) = M - 4.5 - 0.75 log10(t) is widely attributed but its
constants were not verified against the paper; read Figure 6 before using it.

**Forty-one further candidate values were never verified**, because the verification run
was cut short. They are listed nowhere in this document precisely so they cannot be cited
by accident. Do not add citations from memory. If a claim needs support that is not in
the verified list, mark it "[citation needed]" and say what kind of source would carry it.

### 0.10 Data and Resources, draft text

SSA journals require this section. It must be accurate, and reviewers check it.

> The Kahramanmaras catalogue was obtained from the Kandilli Observatory and Earthquake
> Research Institute (KOERI) regional catalogue. The Ridgecrest and Hector Mine catalogues
> were obtained from the U.S. Geological Survey Advanced National Seismic System
> Comprehensive Catalog (ComCat) at https://earthquake.usgs.gov/fdsnws/event/1/ (last
> accessed September 2026), selected within 100 km of each mainshock epicentre at
> M >= 1.0 over 180 days; the exact query strings are distributed with the software in
> `examples/ridgecrest.toml` and `examples/hectormine.toml`. Magnitudes were homogenised
> to Mw following Scordilis (2006). All analyses were performed with Tremor Lab v1.1.1,
> openly available under the MIT licence at https://github.com/HaikKaz/tremor-lab and
> archived at https://doi.org/10.5281/zenodo.22653579; the version used here is
> https://doi.org/10.5281/zenodo.22661044. The cross-check reported in the Discussion used
> `seismostats` (Swiss Seismological Service, ETH Zurich).

Two cautions. Check the KOERI sentence against what can actually be documented; section 13
records that the bundled file is derived rather than raw, and the original export is not
archived. And do not write that the software "reproduces every value reported here" unless
that is still true of the final manuscript.

### 0.11 Figures the paper needs

1. **The three stability curves on one panel**, b against threshold, with the
   maximum-curvature and b-stability thresholds marked on each. This figure carries the
   primary claim on its own, and the Hector Mine plateau is the thing a reader should see
   first.
2. **Completeness against time for Kahramanmaras**, with the two mainshocks marked. This
   carries the secondary claim, and the spike at 0.38 days is the whole point.
3. Optionally, **c against threshold** for both catalogues, which is a clean monotone
   result and supports the incompleteness interpretation of the Omori offset.

Follow the journal's figure conventions. Do not put a figure in that only restates a
table.

### 0.12 Style

The author dislikes text that reads as machine-generated, and reviewers in this field are
suspicious of it. **No em dashes.** Avoid "moreover", "furthermore", "it is important to
note", "underscore", "shed light on", and the "not just X, but Y" construction. Do not end
every paragraph with a sentence restating its opening. Vary sentence and paragraph length.
Prefer a concrete subject and the active voice outside Methods, where the passive is
conventional and fine.

State results in the past tense and established facts in the present. Report effect sizes
and actual values rather than bare significance: "KS 0.0136 on 8,466 events, p = 0.003" is
a result, "significantly different" is not.

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
Mc - dM/2, the lower edge of the completeness bin and the bound the estimator
itself uses in its denominator (section 7 shows why selecting at Mc instead is
biased off the grid), with the Shi and Bolt standard error:

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
`threshold - dM/2`; the package, the spreadsheet and the thesis script took events at or
above `threshold` itself, while still putting `threshold - dM/2` in the estimator's
denominator.

A later audit showed the second of those is the one that is wrong, so **both now select
at `threshold - dM/2`**, and this section has been corrected accordingly. The reason is
worth stating, because it is not obvious. The Aki estimator with Utsu's correction
models magnitudes as an exponential whose lower bound is `Mc - dM/2`; selecting the
sample at `Mc` while telling the estimator the bound is `Mc - dM/2` is a mismatch
between the data and the model. On magnitudes reported on the dM grid the mismatch is
invisible, because no reported value lies between the two bounds. Off the grid it is
not. Simulated with 2,000,000 events drawn from a b = 1.0 law, complete from M 3.0, at
a threshold of 3.5:

| rule | on a 0.1 grid | off the grid |
|---|---|---|
| select at `Mc`, bound `Mc - dM/2` | 0.9952 | **0.8966** |
| select at `Mc - dM/2`, bound `Mc - dM/2` | 0.9952 | **0.9996** |
| select at `Mc`, bound `Mc` (unbinned form) | 1.1240 | 0.9998 |

The two rules are identical on gridded data and the self-consistent one is unbiased on
both, so it is the one to describe in the paper. The third line is the reminder that the
half-bin term is not optional on binned data: dropping it puts b 12 per cent high.

The two agreed **only when the magnitudes and the threshold both sat on the same 0.1
grid** — which is exactly the condition of the locked reference case. Off that grid the
divergence was large: at a threshold of 4.05 on the same catalogue, one implementation
reported 625 events and b = 1.167, the other 498 events and b = 0.929. A 26% difference
in b, invisible to the reference test.

The methodological point: a cross-implementation check performed at a single, tidy
operating point can certify agreement that does not exist elsewhere. Cross-validation
has to be run across the parameter space, not at the published values.

The bug is fixed; the two now agree on off-grid thresholds, non-default bin widths, and
homogenised (therefore off-grid) magnitudes from USGS exports. The command-line report
states the bound it sampled at — "1529 events in its completeness class, at or above
3.45" — so the rule never has to be inferred from a threshold alone.

A second lesson sits underneath the first. This section described the resolution
backwards for a while: the code was corrected and the record was not, so the canonical
handover endorsed the rule the code had been changed away from. A design record is only
as good as the last time somebody checked it against the code, and nothing in the test
suite can check prose.

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

**Read section 9d before using this.** The obvious reading - that the rejection is
caused by the M 7.6 Elbistan event nine hours after the mainshock, giving two
overlapping sequences where one curve is fitted - was tested against a control and
is not supported. A single mainshock with no comparable second event is rejected
just as hard. The rejection is real; the doublet explanation for it is not
established, and the paper must not assert it.

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

## 9b. A fifth finding, methodological: magnitudes are decimal, arithmetic is not

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

## 9c. A second catalogue: Ridgecrest 2019, and what it changes

Everything above was measured on one sequence. On 8 September 2026 the tool was run
on a second, independent one: the 2019 Ridgecrest sequence in California, from the
USGS ComCat catalogue - 28,964 events at M >= 1.0 within 100 km of the M 7.1
mainshock over 180 days. Ridgecrest is also a doublet, an M 6.4 followed 33.8 hours
later by the M 7.1, so it tests the same structure as Kahramanmaras with a different
network, a different agency's export format and a completeness magnitude two
magnitude units lower.

Three results, in order of what they are worth to the paper.

### The single-decay rejection generalises

The transformed-time residual test gives **KS 0.0240 at p = 0.002 +/- 0.002 on 600
replicates: a single modified Omori-Utsu decay is rejected**, more decisively than on
Kahramanmaras (p = 0.036). Two independent doublets, two rejections. Section 9a
argued from one catalogue that one decay curve is the wrong model for two overlapping
sequences; this is the second case, and it turns an observation about one sequence
into a claim that can be made generally.

### A single completeness magnitude is the wrong model for a dense catalogue

This is the more important finding, and it is a limitation of the method rather than
of the software. Maximum curvature over the whole 180 days gives Mc 1.3. Estimated
band by band, the completeness magnitude is nothing like constant:

| time since mainshock | events | Mc | b at that Mc |
|---|---|---|---|
| 0 to 29 minutes | 122 | 3.7 | 0.881 +/- 0.104 |
| 29 minutes to 2.4 hours | 384 | 3.2 | 1.098 |
| 2.4 hours to 1 day | 3,101 | 2.2 | 0.894 |
| 1 to 7 days | 8,123 | 1.6 | 0.909 +/- 0.014 |
| 7 to 30 days | 8,171 | 1.3 | 0.901 +/- 0.011 |
| 30 to 180 days | 9,062 | 1.3 | 1.004 +/- 0.013 |

The network cannot detect small events while the ground is still shaking, so the
first hours are missing exactly the events the decay fit needs. Mc falls from 3.7 to
1.3 - a factor of about 250 in the rate of detectable events.

The consequence for b is large and easy to state. Analysed the ordinary way, at one
Mc for the whole window, **b = 0.733 +/- 0.005**. Estimated in each later band at
that band's own completeness, b is close to **0.94**. The single-Mc figure is low by
0.205, which is **43 times its own quoted standard error**. The Shi and Bolt error
does not know about it, because it describes sampling scatter at one threshold and
nothing else - the same point section 9 makes about the choice of Mc, here in a much
sharper form.

**Read sections 9e and 10c before using the per-band b values in the table above.** The
claim that a single Mc biases b low survives - it is confirmed independently, and the
size of the bias is published. What does not survive is the *trend across the bands*.
Estimating b in each band at that band's own maximum-curvature Mc gives b rising with
time, and the published time-resolved estimates for both sequences run the other way,
because the early bands carry the largest downward bias from exactly the effect this
section is describing. Quote the whole-window figure and the diagnosis; do not present
the band-by-band b column as a measurement of how b evolved.

It also explains the offset. The fit returns **c = 1.415 +/- 0.075 days**, where
values below half a day are usual. A missing early population flattens the start of
the decay, and c absorbs it.

**For the paper:** Kahramanmaras never showed this because the published analysis
used a threshold of M 3.5, far above completeness at every time. A catalogue recorded
to M 1.3 exposes it immediately. Either state the limitation, or restrict the decay
fit to a threshold above the worst-case early completeness - which for Ridgecrest
would be about M 3.7, and would cost most of the data.

### Cross-implementation agreement, on data neither had seen

The package and the browser page were run on the same 28,964 events, both at their
defaults. Not the catalogue either was built against:

| | package | browser page |
|---|---|---|
| events in window | 28,964 | 28,964 |
| Mc | 1.3 | 1.3 |
| events in the completeness class | 20,645 | 20,645 |
| b | 0.733108 +/- 0.004719 | 0.733108 +/- 0.004719 |
| Omori p | 1.172065 | 1.172067 |
| Omori c | 1.414299 | 1.414303 |
| Omori k | 6,659.3 | 6,659.3 |

b and its standard error agree to six decimals; p and c differ in the sixth
significant figure, which is the two optimisers, SciPy's Nelder-Mead against a
hand-written simplex. This is a stronger statement than the reference-catalogue
agreement in section 5, because the reference values were what both implementations
were built and tuned against, and these were not.

One difference is a capability gap rather than a disagreement: the package applied a
100 km distance limit and the browser page has no control for one, so the page
analysed one more event. Say so if the page's numbers are quoted.

### What is not claimed

The completeness and b-value results have since been compared with the published
literature: see section 10b, which confirms the time-varying completeness against
Gulia et al. (2020) using the same method, and the b-value consequence against two
further groups. The Omori parameters remain unchecked against anyone, because the
published Ridgecrest decay work is ETAS-based and its p and c are not the same
quantities as a single modified Omori-Utsu fit.

## 9d. The control that corrects section 9a, and what the test actually measures

Sections 9a and 9c report that a single modified Omori-Utsu decay is rejected for
two doublets, Kahramanmaras and Ridgecrest, and read that as evidence that one
curve cannot describe two overlapping sequences. That reading needed a control:
short-term aftershock incompleteness and secondary aftershocks would produce the
same rejection in *any* sequence, doublet or not.

The control is the 1999 Hector Mine earthquake, M 7.1 in California - a clean
single mainshock with no M 6 or larger event anywhere in the following fortnight.
Landers 1992 was considered first and rejected as a control, because the M 6.3 Big
Bear event three hours later makes it a doublet too.

USGS ComCat, 13,525 events at M >= 1.0 within 100 km over 180 days, Mc 1.7:

| sequence | structure | n | KS | 95th pct of the null | p | verdict |
|---|---|---|---|---|---|---|
| Kahramanmaras 2023 | doublet | 1,529 | 0.0210 | 0.0194 | 0.033 | rejected |
| Ridgecrest 2019 | doublet | 20,644 | 0.0240 | 0.0054 | 0.003 | rejected |
| **Hector Mine 1999** | **single** | 8,466 | **0.0136** | 0.0084 | **0.003** | **rejected** |

**The single mainshock is rejected as decisively as the doublets.** The doublet
explanation in section 9a is therefore not supported, and the paper must not make it.

The table also shows what the test is really responding to. The critical value falls
from 0.0194 to 0.0054 as the catalogue grows from 1,529 events to 20,644, roughly as
one over the square root of n. Hector Mine has the **smallest** departure from the
model of the three - KS 0.0136 against Kahramanmaras's 0.0210 - and is rejected far
more decisively, because its sample is five times larger. Kahramanmaras's departure
would be overwhelmingly rejected at Ridgecrest's sample size; Hector Mine's would
pass at Kahramanmaras's.

That is the ordinary large-sample behaviour of a goodness-of-fit test, and it has a
straightforward reading here. The modified Omori-Utsu law describes one generation of
aftershocks. Real sequences contain aftershocks of aftershocks, which is the entire
reason ETAS exists, so the law is known to be an approximation. With a few thousand
events that approximation is detectable, and the test detects it. Rejection at n =
20,000 is not evidence that a sequence is unusual; it is evidence that the catalogue
is large.

### What to write instead

1. **Report the effect size, not only the verdict.** KS 0.0136 to 0.0240 across three
   sequences means the fitted curve tracks the data closely in all of them. A table
   of KS, n and p says something; "rejected" on its own does not.
2. **Do not attribute the rejection to the doublet structure.** The control rules
   that out as a sufficient explanation. If the paper wants that claim it needs a
   different test - fitting the two sequences separately and showing the residuals
   improve, which this tool does not do.
3. **The honest general statement** is that the modified Omori-Utsu law is a good but
   imperfect description of all three sequences, and that a residual test on a modern
   catalogue of thousands of events will detect the imperfection. That is a useful
   methodological point about the test, and it is defensible.
4. **The large c values point at incompleteness.** Ridgecrest returns c = 1.415 days
   and Hector Mine c = 4.814 days, where below half a day is usual. Section 9c shows
   completeness moving by 2.4 magnitude units over the first month at Ridgecrest. A
   missing early population flattens the start of the decay and c absorbs it. That is
   a more likely explanation for the misfit than sequence structure, and it is
   testable: refit above the early-time completeness magnitude and see whether c falls.

This control cost twenty minutes and changed a headline finding from something a
referee could have dismantled into something defensible. It is worth saying in the
paper that it was done.

## 9e. Kahramanmaras: completeness moves, and the thesis's b-value is measured below the stable range

Section 9c found on Ridgecrest that completeness moves sharply after a mainshock. The
same measurement on the thesis's own catalogue confirms it, and then leads somewhere
more consequential: **the b-value the thesis reports is estimated at a threshold where
b has not stabilised.** A control settles which of the two possible explanations holds.

### Completeness collapses twice, once per mainshock

Bundled KOERI catalogue, 3,469 events, M 3.0 to 7.6, maximum curvature in each band:

| time after the M 7.8 | events | Mc |
|---|---|---|
| 0 to 0.20 days | 155 | 3.8 |
| 0.20 to 0.38 days | 106 | 3.7 |
| **0.38 to 0.60 days** (spans the M 7.6 at 0.380 d) | 120 | **4.4** |
| 0.60 to 1.0 days | 214 | 3.6 |
| 1 to 2 days | 347 | 3.4 |
| 2 to 5 days | 497 | 3.2 |
| 30 to 180 days | 968 | 3.3 |

**The doublet produces two degradations, and the second is worse than the first.**
Completeness falls to M 3.8 after the M 7.8, is recovering by nine hours, then collapses
again to M 4.4 when the M 7.6 Elbistan event arrives - worse, because the network is
already saturated. This is a clean, specific observation about this sequence, and it is
something a doublet does that a single mainshock cannot. It is worth reporting on its
own.

The published analysis uses a fixed **M 3.5** threshold for all 180 days, so for roughly
the first day it is analysing below completeness, and 425 events enter the sample from a
period the catalogue could not support.

### But the first day is not the explanation

The obvious move is to drop the incomplete head of the sequence and re-estimate. Doing
that does shift b, from 0.844 to 0.947, which is 5.5 standard errors. It is tempting to
stop there and report the higher number.

**That would be wrong.** If the first day were the problem, then after removing it b
would no longer depend on where the magnitude threshold is put. It still does:

| start | thr 3.5 | thr 3.8 | thr 4.0 | thr 4.4 |
|---|---|---|---|---|
| 0 (as published) | 0.844 +/- 0.019 (1529) | 0.985 (954) | 1.028 (625) | 1.124 (258) |
| 1 day | 0.947 +/- 0.024 (1104) | 1.140 (652) | 1.178 (393) | 1.211 (139) |
| 2 days | 0.944 (937) | 1.156 (560) | 1.213 (339) | 1.275 (117) |
| 7 days | 0.923 (641) | 1.120 (385) | 1.166 (236) | 1.249 (88) |

The climb with threshold is present at every start time, and is if anything steeper
after the time cut. Cutting time and raising the threshold remove overlapping sets of
events, so both shift b for the same underlying reason. Presenting 0.947 as the
corrected value would just be reporting a second arbitrary point on the same slope.

### b does not stabilise until about M 4.1

The b-stability method of Cao and Gao (2002) and Woessner and Wiemer (2005) exists for
exactly this, and the package implements it. Over the whole window:

| Mc | 3.2 | 3.4 | **3.5** | 3.8 | 4.0 | **4.1** | 4.3 | 4.5 |
|---|---|---|---|---|---|---|---|---|
| b | 0.829 | 0.851 | **0.844** | 0.985 | 1.028 | **1.041** | 1.048 | 1.137 |
| +/- | 0.014 | 0.018 | **0.019** | 0.030 | 0.040 | **0.046** | 0.058 | 0.087 |

b rises monotonically by about 0.2 and only flattens near the top. `mc_b_stability`
returns **Mc 4.1**. Hainzl, Kumazawa and Ogata (2024) report exactly this for the same
region from independent data - b stable only above a cut-off of 3.5, scattering around
1.07 - see section 10c. The published M 3.5 sits in the middle of the steep part - the
threshold is 0.6 magnitude units below where b stops moving.

The three estimators disagree accordingly, and now it is clear which is which: maximum
curvature **3.4**, goodness-of-fit **3.0**, b-stability **4.1**. Section 9 already
reported that spread; this says the top of it is the one to believe.

### The control: this is the estimator, not the sequence

Two readings are available. Either maximum curvature is underestimating completeness, or
b genuinely varies with magnitude in this sequence. Hector Mine settles it - a single
mainshock, a different network, a different decade, 13,525 events:

| Mc | 1.5 | **1.7 (maxcurv)** | 1.9 | **2.1 (b-stab)** | 2.2 | 2.4 | 2.6 | 2.8 | 3.0 |
|---|---|---|---|---|---|---|---|---|---|
| b | 0.764 | **0.839** | 0.879 | **0.985** | 1.003 | 0.990 | 0.968 | 0.980 | 0.989 |
| n | 11,163 | 8,466 | 5,965 | 4,345 | 3,526 | 2,217 | 1,384 | 900 | 585 |

**The same climb, and then an unambiguous plateau** - b sits at 0.98 +/- 0.02 across nine
consecutive thresholds spanning 0.9 magnitude units, on thousands of events. That plateau
is what a correctly-thresholded catalogue looks like, and it proves the shape is not an
artefact of running out of data.

Maximum curvature returns 1.7 there, where b is 0.839. The stable value is 0.985. **The
estimator underestimates completeness by about 0.4 magnitude units and biases b low by
about 0.15, in a sequence with no doublet and no unusual feature at all.** Kahramanmaras
shows the same thing, at the same size, for the same reason. This is known behaviour -
maximum curvature is documented as biased low, and the +0.2 correction this package
applies is an attempt to patch it that evidently does not go far enough here.

Note also where both sequences land: 0.98 and about 1.04, either side of the canonical
Gutenberg-Richter value of 1.0, reached independently from two unrelated catalogues.

### What to write

1. **Report the two-stage completeness collapse.** It is measured, specific to this
   doublet, and nothing else in this work is as clean.
2. **Do not report a single b without saying what threshold it is at.** b = 0.844 is the
   value at M 3.5. It is real, correctly computed and reproducible - and it is measured
   0.6 units below where b stops changing.
3. **Give the stability curve, not one number.** The honest statement is that b rises
   from 0.83 to about 1.04 as the threshold rises to M 4.1 and is flat above it, and that
   the Shi and Bolt error - 0.019 at M 3.5 - describes sampling scatter at a fixed
   threshold and is silent about a 0.2 shift driven by threshold choice. This is the same
   point section 9 makes, in its sharpest available form.
4. **Cite the control.** A referee will ask whether the effect is peculiar to a doublet.
   It is not, and Hector Mine shows so on 13,525 events.
5. **Do not claim b is "really" 1.04.** What is shown is that the estimate depends on the
   threshold and stabilises near 1.0. Whether the true b is that value depends on whether
   the catalogue is complete at M 4.1, which this data cannot establish.

### Three catalogues, three different failure modes

Running the same stability curve on all three sequences separates two effects that had
been conflated, and shows the method failing in three distinguishable ways.

| | events | maxcurv Mc | b-stability Mc | b at maxcurv | b where stable |
|---|---|---|---|---|---|
| Hector Mine (single) | 13,525 | 1.7 | **2.1** | 0.839 | 0.985, flat over 0.9 units |
| Kahramanmaras (doublet) | 3,469 | 3.4 | **4.1** | 0.851 | 1.041, flattening at the data edge |
| Ridgecrest (doublet) | 28,963 | 1.3 | **None** | 0.733 | never stabilises |

**Ridgecrest never stabilises at all.** b climbs continuously from 0.695 at Mc 1.1 to
1.341 at Mc 3.5 without a plateau anywhere, and `mc_b_stability` correctly returns None
rather than inventing a threshold. That is the right behaviour and worth saying: the
estimator refuses when its own premise fails.

### This confirms section 9c rather than replacing it

Section 9c attributed Ridgecrest's low single-Mc b of 0.733 to averaging over a window
in which completeness moved by 2.4 magnitude units. Section 9e's threshold effect is a
different mechanism, so the two had to be separated. Varying both:

| start | thr 1.3 | thr 1.8 | thr 2.2 | thr 2.6 |
|---|---|---|---|---|
| 0 | 0.733 | 0.794 | 0.820 | 0.865 |
| 1 day | 0.886 | 0.988 | 1.016 | 1.066 |
| 7 days | 0.950 | 1.003 | 1.050 | 1.129 |

Dropping the first day raises b by about 0.15 **at every threshold**, and raising the
threshold raises it by about 0.15 **at every start time**. The two effects are roughly
additive and roughly equal here, so **both sections are right and neither subsumes the
other**: Ridgecrest suffers from time-averaging and from an under-set threshold at once,
which is why it never reaches a plateau.

Kahramanmaras is the opposite case - there the time cut and the threshold cut remove
overlapping sets of events, so they are not separable, and the threshold effect is the
one that survives. The practical rule for the paper is that **b needs both its threshold
and its time window stated**, and that the stability curve should be shown for each.

### The decay fit: the prediction half-holds

Section 9d predicted that the large fitted c comes from a missing early population, and
that refitting above the early-time completeness would bring c down. Tested on both:

| | thr | p | c | KS | p-value |
|---|---|---|---|---|---|
| Kahramanmaras | 3.5 | 1.161 | **0.497** | 0.0210 | 0.042 |
| | 3.8 | 1.147 | 0.364 | 0.0336 | 0.002 |
| | 4.0 | 1.143 | 0.255 | 0.0497 | 0.002 |
| | 4.4 | 1.195 | **0.195** | 0.0742 | 0.002 |
| Hector Mine | 1.7 | 1.086 | **4.814** | 0.0136 | 0.002 |
| | 2.0 | 1.204 | 2.808 | 0.0109 | 0.080 |
| | 2.2 | 1.276 | 2.093 | 0.0181 | 0.002 |
| | 2.6 | 1.276 | **0.710** | 0.0235 | 0.032 |

**c falls monotonically in both, by a factor of 2.5 and of 6.8.** The prediction about c
is confirmed, and confirmed on a control as well as on the doublet - c is absorbing
missing early events, exactly as section 9d argued.

**The fit is not rescued.** The KS departure grows as the threshold rises and the
rejection mostly stands. Say the incompleteness explanation accounts for the offset
parameter and not for the misfit, and leave it there. Note too that p drifts upward with
threshold in both - by 0.03 and 0.19 - so p is threshold-dependent as well, and any
reported p needs its threshold stated for the same reason b does.

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
which on this catalogue reads 0.506 in the package and 0.508 in the browser (they
use different optimisers, so the fitted curve differs in the fourth decimal) where
the calibrated one reads 0.036: the
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

## 10b. The published literature on Ridgecrest, and what it confirms

Section 9c said the Ridgecrest results had not been compared with the published
literature and must not be presented as agreeing with anyone. They have now. Every
citation below was fetched and checked; the one that carries the most weight was
read directly from the paper by the author of this document rather than taken on
trust, and is marked.

### The time-varying completeness result is confirmed, by the same method

**Gulia, L., Wiemer, S. and Vannucci, G. (2020),** *Pseudo-prospective evaluation of
the Foreshock Traffic Light System in Ridgecrest and implications for aftershock
hazard assessment*, Seismological Research Letters 91(5), 2828-2842,
doi:10.1785/0220190307. **Verified first-hand:** the PDF was downloaded from
zenodo.org/records/5076000 and the text extracted. Their line 234 states they use
"the Maximum Curvature method (Wiemer and Wyss, 2000)" - the estimator this package
implements. Lines 187-190 read, verbatim:

> "after the Mw6.4 Mc increased from the background value (Mc=1.2) to about 1.8,
> before dropping back to a near-to-background value within 12 hours. After the
> Mw7.1 event, it increased to between 3.3 and 3.5, then recovered within three days
> to near-to-background values."

Against section 9c, measured independently with this tool on USGS ComCat:

| time after the M 7.1 | this tool | Gulia et al. (2020) |
|---|---|---|
| first half hour | Mc 3.7 | rises to 3.3-3.5 |
| 29 min to 2.4 h | Mc 3.2 | |
| 2.4 h to 1 day | Mc 2.2 | |
| 1 to 7 days | Mc 1.6 | recovered near background within 3 days |
| 7 days on | Mc 1.3 | background Mc 1.2 |

The same phenomenon, the same method, the same magnitudes. Our peak is 0.2 to 0.4
higher than theirs, which is what a narrower early window should give: they use
300-event moving windows, we use a fixed 29-minute band. **This is independent
published confirmation of the finding in section 9c, and it should be cited rather
than the finding presented as new.** What is new here is not that completeness
varies - that is established - but the size of the consequence for b, below.

### The b-value consequence is corroborated

Two independent groups report the same time dependence in b that section 9c
attributes to the missing early population.

- **Hainzl, S. (2022),** *ETAS-Approach Accounting for Short-Term Incompleteness of
  Earthquake Catalogs*, BSSA 112(1), doi:10.1785/0120210146. Reports the apparent b
  below 0.5 immediately after the mainshock, rising to about 1, converging after
  about a day for Ridgecrest - and attributes the depression to missing small events
  rather than to a real change in b.
- **Huang, K., Tang, L. and Feng, W. (2022),** on the SCEDC catalogue with an
  Aki/Utsu maximum-likelihood estimator with the half-bin correction - the same
  estimator as this package - report mean b recovering from 0.59 to 0.86.

Our single-Mc figure of 0.733 sits between the depressed early values and the
recovered ones, which is exactly what an average over a period of changing apparent
completeness should give, and our band-by-band values of 0.90 to 1.00 sit with the
recovered ones. **Do not present 0.733 as a measurement of b for this sequence.**
It is an artefact of applying one completeness magnitude to a window in which
completeness moved by 2.4 magnitude units, and the literature says so independently.

### Whole-catalogue completeness

Two studies report a maximum-curvature Mc for the Ridgecrest region close to ours.
These were verified by the search agents rather than first-hand, so treat them as
indicative until checked:

| source | catalogue and window | Mc |
|---|---|---|
| this tool | USGS ComCat, 180 days from the M 7.1, 100 km | 1.3 |
| Herrmann and Marzocchi (2021), SRL 92(2) | SCSN, Apr-Dec 2019, 100 km | MAXC+0.2 = 1.10 |
| Huang, Tang and Feng (2022) | SCEDC, Jul-Dec 2019 | MAXC = 1.10 |

One bin apart, on a different catalogue over a different window - ours begins at the
mainshock and so contains proportionally more of the incomplete early period, which
raises Mc. Comparable, not identical, and the difference is explained.

### What the literature says about the decay-fit rejection

Section 9d concluded from the Hector Mine control that rejection of a single
Omori-Utsu decay is not evidence of doublet structure. The literature supports the
more specific reading offered there. Work on Ridgecrest routinely uses ETAS rather
than a single Omori law, and the incompleteness-corrected variants exist precisely
because the uncorrected fit is biased: Hainzl's ETASI and the PETAI approach of
**Mizrahi, L., Nandan, S. and Wiemer, S. (2021),** *Embracing Data Incompleteness for
Better Earthquake Forecasting*, JGR Solid Earth, doi:10.1029/2021JB022379, which
reports that fifteen minutes after the M 7.1 an event below M 3.0 had almost no
chance of being detected and M 3.5 roughly even odds. That is the missing population
section 9d blames for the large fitted c.

### References to cite for short-term aftershock incompleteness

The phenomenon has a name and a literature; the paper should cite it rather than
describe it from scratch.

- **Kagan, Y. Y. (2004),** *Short-Term Properties of Earthquake Catalogs and Models of
  Earthquake Source*, BSSA 94(4), 1207-1228, doi:10.1785/012003098. The standard
  early reference. **No numeric value has been read from it** - the full text was not
  reachable - so cite it for the phenomenon, not for a number.
- **Helmstetter, A., Kagan, Y. Y. and Jackson, D. D. (2006),** BSSA 96(1), 90-106,
  doi:10.1785/0120050067. Usually cited for the empirical relation
  mc(t) = M - 4.5 - 0.75 log10(t) for southern California. **The constants were not
  verified first-hand** - the full text is paywalled and the attribution comes from
  later papers citing it. If that formula enters the manuscript, read Figure 6 first.
- **Lippiello, E., Cirillo, A., Godano, C., Papadimitriou, E. and Karakostas, V.
  (2019),** *Post Seismic Catalog Incompleteness*, Geosciences 9(8), 355,
  doi:10.3390/geosciences9080355. Open access, and the cleanest available statement
  of the standard parameterisation.

### What is still not compared

The Omori parameters. Published Ridgecrest decay work is almost entirely ETAS, whose
p and c are not the same quantities as a single modified Omori-Utsu fit and must not
be set beside ours in a table. A like-for-like comparison would need a published
single-Omori fit on a comparable window and threshold, and none was found. Say the
comparison covers completeness and the b-value, as with `seismostats` in section 10a,
and that the decay parameters remain unchecked against anyone.

## 10c. The published Kahramanmaras literature: three confirmations and one correction

Section 10b did this for Ridgecrest. This does it for the thesis's own sequence. Fifty-one
candidate published values were found; **ten were verified first-hand** - the source
located, the paper read, the number and its conditions quoted verbatim. The remaining
forty-one were still being checked when the run hit a usage limit, so **they are
unverified, not rejected, and must not be cited from this document.**

Three of the ten bear directly on section 9e. One contradicts a result in section 9c.

### 1. A like-for-like match at the same threshold, from an independent catalogue

**Ali, S. M. and Abdelrahman, K. (2024),** *Analysis of the Fractal Dimension, b-value,
Slip Ratio, and Decay Rate of Aftershock Seismicity Following the 6 February 2023
(Mw 7.8 and 7.5) Turkiye Earthquakes*, Fractal and Fractional 8(5), 252,
doi:10.3390/fractalfract8050252. IRIS catalogue, 471 events M 3.3 to 7.8, 6 Feb 2023 to
10 Jan 2024, 35-39 N / 34-41 E. Maximum-likelihood fits in ZMAP 6.0. Their maximum
curvature returns **Mc 4.4** for the whole sequence.

That is the same threshold section 9e reaches by a different route, so the two are
directly comparable:

| at Mc 4.4 | Ali and Abdelrahman (2024) | this tool |
|---|---|---|
| b | 1.21 +/- 0.03 | 1.124 +/- 0.073 |
| p | 1.1 +/- 0.04 | 1.195 |
| **c** | **0.204 +/- 0.058** | **0.195** |
| k | 76.75 +/- 8.84 | 49.9 |

**c agrees to 0.009 days**, and b and p agree within about one standard error. Different
catalogue, different agency, different software, same threshold, same estimators. This is
the strongest external check in this document - stronger than the `seismostats`
cross-check in section 10a, which used the same input data.

Three caveats. **k is not comparable** and must not be tabulated: it scales with the
number of events above threshold, and their 471-event catalogue is not ours. The paper's
b uncertainty is internally inconsistent - the abstract says +/- 0.1, the body +/- 0.03;
prefer the body figure or note the conflict. And the paper misreads its own k as a time
in days ("the decrease in aftershock activity began between approximately 68 to 86 days
after"); k is a productivity constant, so quote their numbers but not that sentence.

### 2. The b-stability threshold is independently confirmed, in this same region

This is the citation section 9e most needs.

**Hainzl, S., Kumazawa, T. and Ogata, Y. (2024),** *Aftershock forecasts based on
incomplete earthquake catalogues: ETASI model application to the 2023 SE Turkiye
earthquake sequence*, Geophysical Journal International 236(3), 1609-1620,
doi:10.1093/gji/ggae006. On AFAD background seismicity for the same region from 2000 to
2022, using the Aki (1965) maximum-likelihood estimator - the same estimator this package
implements - they plot b against cut-off magnitude and report, verbatim:

> "The b estimate becomes stable within its uncertainties for M c > 3.5 and scatters
> around 1.07."

Section 9e, measured independently on the aftershock catalogue, finds b stabilising at
**Mc 4.1** at a value of **1.04**. An independent group, on different data spanning
twenty-three years, reaches the same two conclusions: b in this region only stabilises
above a high cut-off, and it stabilises near 1.05. **Cite this rather than presenting
section 9e as a new observation.** What is new in 9e is the Hector Mine control showing
the effect is a property of the estimator rather than of this region.

### 3. The size of the incompleteness bias is published, on identical data

The same paper fits a standard ETAS model and an incompleteness-corrected ETASI model to
one catalogue (9,438 events, m >= 2, Mc 1.95, estimated network blind time 162 s):

| | b |
|---|---|
| ETAS, uncorrected | 0.56 |
| ETASI, corrected for short-term incompleteness | **0.87** |

**A downward bias of 0.31 on identical data, attributable entirely to unhandled
incompleteness.** Section 9e measures 0.15 to 0.20 from threshold choice alone on two
catalogues. Same direction, same order of magnitude, arrived at by a completely different
method. This is the strongest published support for the paper's central methodological
claim.

### 4. Time-varying completeness is confirmed, but the absolute numbers are not comparable

**Tan, O. (2025),** *Long-term Aftershock Properties of the Catastrophic 6 February 2023
Kahramanmaras (Turkiye) Earthquake Sequence*, Acta Geophysica 73, 1023-1040,
doi:10.1007/s11600-024-01419-y. AFAD national catalogue, ML, 50,085 events, 6 Feb to 31
Oct 2023:

> "decreases gradually from ~ 3.0 on the first day to ~ 2.0 two weeks later ... The Mc
> converges to the constant value of ~ 1.5 after mid-March and is 1.6 for the nine-month
> catalog data."

The same shape as section 9e - high on the first day, recovering over weeks to a
constant - from an independent group on an independent catalogue. **The absolute values
are not comparable and must not be tabulated side by side.** Their first-day Mc is 3.0
against our 3.8, and their long-window Mc is 1.6 against our 3.4, because the AFAD
catalogue is complete far below M 3.0 while the bundled KOERI extract is truncated at
M 3.0. For our catalogue the truncation floor, not the network, sets the long-window
figure. Cite Tan for the phenomenon and the recovery timescale, not for the numbers.

### 5. The correction: our band-by-band b trend has the wrong sign

This is the finding that costs something, and it must go in the paper.

The same Hainzl, Kumazawa and Ogata (2024) paper applies the Ogata-Katsura (1993)
estimator, which fits a time-varying b jointly with a time-varying detection function
instead of assuming a fixed threshold. It reports **b about 1.2 for the early aftershocks,
decaying to about 0.85 for the later ones** - a 50 per cent coseismic increase in b, which
is also the premise of the Gulia and Wiemer Foreshock Traffic Light System cited in
section 10b.

Sections 9c and 9e report the opposite. Estimated band by band at each band's own maximum
curvature Mc, this tool gives:

| | early | late |
|---|---|---|
| Kahramanmaras | 0.737 +/- 0.080 (29 min - 2.4 h) | 1.018 +/- 0.040 (30-180 d) |
| Ridgecrest | 0.881 +/- 0.104 (first 29 min) | 1.004 +/- 0.013 (30-180 d) |

**b rising with time, where the literature has it falling.** The explanation is the one
section 9e establishes, applied to our own numbers: an Aki estimate at a threshold below
true completeness is biased *low*, and completeness is *worst* in the early bands - so the
early bands carry the largest downward bias. The Kahramanmaras early figure rests on 72
events with an Mc itself estimated from those 72 events. Hainzl et al. avoid this by
modelling detection probability explicitly rather than thresholding.

**Do not present the band-by-band b trend as a result.** State that estimating b in each
band at that band's maximum-curvature Mc produces a trend of the wrong sign, that the
published time-resolved estimates run the other way, and that this is a demonstration of
the very bias the paper is about rather than a competing measurement. It is a better
illustration of the point than a clean agreement would have been - the tool reproduces the
artefact, and the diagnosis explains it.

### 6. The two ruptures do decay differently, which section 9d should acknowledge

Section 9d showed that rejection of a single Omori-Utsu decay is not by itself evidence of
doublet structure, since the single-mainshock Hector Mine control is rejected as hard. That
stands. But two groups fit the two ruptures separately and find they differ:

| | EAFZ (Pazarcik, Mw 7.8) | Cardak (Elbistan, Mw 7.6) |
|---|---|---|
| Tan (2025), Mc 2.0 | p = 1, c = 5 d, k = 2600 | p = 0.7, c = 5 d, k = 930 |
| Rodriguez-Perez and Zuniga (2025) | p = 1.25 +/- 0.08 | p = 1.14 +/- 0.09 |

Tan additionally reports that the single-Omori model fits the Cardak sequence poorly after
mid-June. Our combined p = 1.161 sits between both published pairs, which is what a fit to
the union of two sequences should give.

So the honest position is narrower than either extreme: **the doublet structure is real and
published fits resolve it, but the KS rejection is not the evidence for it** - a control
without a doublet rejects too. If the paper wants the doublet claim, cite Tan and
Rodriguez-Perez for it rather than the residual test.

**Rodriguez-Perez, Q. and Zuniga, F. R. (2025),** *Statistical and source characterization
of the 2023 Kahramanmaras Turkiye earthquake sequence*, Acta Geophysica 73, 1241-1260,
doi:10.1007/s11600-024-01428-x. Note for comparison purposes that this paper states no
numeric Mc anywhere, which limits what can be done with its b-values.

### What is still not compared

Forty-one of the fifty-one candidate values were never verified. Two further papers were
read and are usable but bear less directly: **Convertito, V., Tramelli, A. and Godano, C.
(2024),** Scientific Reports 14, 1596, doi:10.1038/s41598-023-50837-3 (on-fault b of 0.7
to 0.8 at the Mw 7.8 nucleation segment; most grid cells mc 1.5 to 2.2), and Tan's
whole-catalogue b of about 0.8 at Mc about 1.6 - which is *not* comparable to our 0.844 at
M 3.5 despite the numerical coincidence, because section 9e shows b at those two thresholds
should not be equal.

## 11. Facts and figures

- Package: `tremor-lab` 1.1.1, MIT licence, Python 3.11+, released as `v1.1.1` and
  archived at doi.org/10.5281/zenodo.22661044.
- Runtime dependencies: NumPy (>=2.0,<3), SciPy (>=1.13,<2), pandas (>=2.2,<3). Nothing else.
- Source: about 2,700 lines across 10 modules. **294 tests**, green in continuous
  integration on Python 3.11, 3.12 and 3.13.
- Browser page: one file, about 2,800 lines, 205 KB, pure ASCII. The only things it
  fetches are the Google Fonts stylesheet and the font files it points at; with no network it
  falls back to system fonts and every number is still computed. Runs from a
  double-click, with nothing installed. Eight settings are editable on it:
  the window, the threshold, the bin width, the Mc correction, the Bath deficit, the
  bootstrap count, the fit-test replicate count and the seed. The rest of the published
  constants are fixed in the page and adjustable in the package, which is the authority
  for published values in any case.
- Public API: 30 names, including `b_stability`, `mc_b_stability`,
  `mc_goodness_of_fit`, `omori_fit_test` and `b_value_tinti`. 25 published
  constants, of which three - the two Omori plausibility bounds and the offset
  floor - judge a fit rather than entering one and so have no keyword argument. The
  other 22 are overridable three ways: per call by keyword, per session by
  reassignment, or from a settings file. The three bounds
  (`OMORI_C_FLOOR`, `OMORI_P_MIN`, `OMORI_P_MAX`) take the latter two only.
- **Speed: do not quote a figure from this handover.** The tool analyses the
  3,469-event reference catalogue and draws all four charts in well under a second in a
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
   uncertainty in b.** It is the sampling scatter at one threshold. The completeness
   magnitude is uncertain from M 3.0 to M 4.1 by method, and b climbs from 0.83 to
   1.04 across that range without stabilising until the top of it. Sections 9, 9e and
   10c; Hainzl, Kumazawa and Ogata (2024) report the same threshold dependence for
   this region independently.
5. **Do not quote k without a caveat, and do not quote a 43,000-event figure** — no such
   case exists in the repository. Of the three decay parameters k is much the weakest: it
   is derived from the fitted c and p rather than searched for. The only k-recovery test
   in the repository (`tests/test_omori.py::test_the_productivity_constant_is_recovered`)
   draws 8,678 synthetic events from p = 1.15, c = 0.5, k = 2000 and asserts only that k
   returns within 10%. At that test's own seed k comes back 1.3% low; across seeds 0-4 the
   k error reaches 6.4% while p stays within 1.4%. p is the robust number.
6. **The DOI exists; cite the right one.** Cite the concept DOI,
   10.5281/zenodo.22653579, which always resolves to the latest version. The DOI of
   the current release is 10.5281/zenodo.22661044 (v1.1.1); use a version DOI only
   where the exact release matters, and use 1.1.1 rather than 1.1.0 - the estimators
   are the same, but 1.1.0 mis-reads a European CSV export and takes the b-value
   sample at a bound its own estimator does not assume. Do not cite the
   GitHub URL in place of the DOI: a repository can be renamed or deleted.
7. **Do not claim the decay-fit rejection is caused by the doublet structure.** It
   was tested: a single mainshock is rejected just as hard (section 9d). Report the
   rejection with its KS statistic and sample size, and say that the modified
   Omori-Utsu law is an approximation which a large catalogue detects.
8. **Do not describe the probabilistic anomaly method** (spatial rarity, the simulation).
   It is deliberately outside this tool and is a separate layer to be added once its
   specification is frozen.
9. **Do not say the estimators are novel.** See section 2.
10. **Do state the third-party comparison, and state its limits.** The estimators
   have now been checked against `seismostats` (Swiss Seismological Service, ETH
   Zurich): see section 10a. What has been compared is the completeness magnitude,
   the b-value and the b-stability criterion. The Omori-Utsu decay, the energy
   relations and the Bath screen have **not** been compared with anything
   third-party, because `seismostats` does not implement them. Do not let the
   b-value comparison imply the decay fit was independently checked.
11. **Do not let the self-test stand for validation of a reader's own analysis.** It
   checks the estimators against known values on a bundled catalogue. It says nothing
   about whether the reader chose a sensible window, threshold or mainshock.
12. **Do not present "294 tests" as coverage.** It is a count, not a measure. What can
    honestly be said is stronger and more specific: ten deliberate breakages of the
    estimators were each caught by at least one test (section 10).
13. **Do not describe the bundled fixture as raw data.** `kahramanmaras_180d.csv` is a
    derived file: it carries elapsed days and magnitudes only, already windowed to 180
    days, with no timestamps. Its provenance from the original KOERI export should be
    stated in the paper, and the export itself archived with the release.
14. **Never report a b-value, or an Omori p, without the threshold it was measured
   at.** Both are threshold-dependent in every catalogue tested: b by about 0.2 and p
   by up to 0.19 between Mc and Mc + 1. A bare "b = 0.844" is not a reproducible
   statement about the sequence, it is a statement about a choice. Section 9e.
15. **Do not present the band-by-band b values as a measurement of how b evolved.**
   Estimated at each band's own maximum-curvature Mc they trend the wrong way against
   the published time-resolved estimates, because the early bands carry the largest
   downward bias. Report them as a demonstration of the bias, with the published
   result cited, or not at all. Sections 9c and 10c.
16. **Do not claim the Scordilis relations are range-checked.** They are applied without
   bounds, matching the spreadsheet implementation, and the two Ms branches are mildly
   discontinuous across the uncalibrated 6.1-6.2 gap (6.157 against 6.119 at Ms 6.1).
   This is a deliberate fidelity choice and is documented.

## 13. Status and what is outstanding

Everything is done. The software is public, released, archived and citable.

- **Repository:** github.com/HaikKaz/tremor-lab, current release `v1.1.1`.
- **DOI:** concept 10.5281/zenodo.22653579, which always resolves to the latest
  version and is the one to cite. Version 1.1.1 is 10.5281/zenodo.22661044; version
  1.1.0 is 10.5281/zenodo.22653580 and should not be cited in preference to it. The
  1.1.0 archive also carries a Software Heritage identifier,
  `swh:1:dir:d09344f80b18a2ba5ad8408d0140ee9450182652`.
- **ORCID:** 0009-0007-1842-1590, affiliation National Academy of Sciences of
  Armenia, both in `CITATION.cff` and in the Zenodo metadata.
- **Continuous integration** runs and passes: `.github/workflows/tests.yml` lints,
  runs the suite on Python 3.11, 3.12 and 3.13, reproduces the reference values
  both directly and through the full catalogue pipeline, and runs the `seismostats`
  cross-check in a separate job. This is worth a sentence in the paper: it means the
  locked values reproduce on a clean checkout with freshly resolved dependencies, on
  three Python versions, and not only on the author's machine.

**One thing for the author, not the writer.** The name is **Aik Kazarian**, as on
the passport and the ORCID record, and that is what `CITATION.cff` and the Zenodo
metadata carry. Earlier thesis drafts use "Haik". Use "Aik" in the byline, the
author list and the availability statement, so the paper, the ORCID record and the
software archive agree. A single author spelled two ways across a paper and its
cited software is how one citation record becomes two.

The availability statement can now be written. One that is accurate:

> *Tremor Lab v1.1.1 is openly available under the MIT licence at
> https://github.com/HaikKaz/tremor-lab and archived at
> https://doi.org/10.5281/zenodo.22653579. The version used in this work is
> v1.1.1 (https://doi.org/10.5281/zenodo.22661044). It reproduces every value
> reported here from the bundled catalogue with a single command, and its test
> suite runs on Python 3.11 to 3.13 in continuous integration.*

Two cautions on that wording. "Reproduces every value reported here" is true of the
values in section 4 and false of anything computed outside this tool, so check it
against the final manuscript before using it. And it says the suite *runs*, not that
the software is correct; section 12 lists what may not be claimed.

## 14. Limitations to state in the paper

The estimators inherit their known behaviour. **Maximum-curvature Mc is biased low**,
and sections 9e and 10c measure how much that costs: on Hector Mine it returns 1.7 where
b only stabilises at 2.1, and the b estimated at its answer is 0.15 below the stable
value. The +0.2 correction this package applies does not close the gap. Mc is also
sensitive to binning and to short-term aftershock incompleteness in the hours after a
mainshock - section 9c measures that on Ridgecrest, where completeness moves from M 3.7
in the first half hour to M 1.3 after a month.

The Aki b-value assumes completeness above the chosen threshold, and **b is
threshold-dependent in every catalogue tested** - by about 0.2 between Mc and Mc + 1,
which is an order of magnitude larger than the Shi and Bolt error at either end. The
same is true of Omori p, by up to 0.19. Neither should be reported without its
threshold. The Omori fit needs a dense
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
