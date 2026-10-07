"""Gutenberg-Richter b-value and its uncertainty."""

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import chi2 as chi2_distribution
from scipy.stats import norm

from tremor_lab import constants
from tremor_lab.completeness import fmd
from tremor_lab.grid import at_or_above


class BValue(NamedTuple):
    """b-value, its standard error, and the number of events it was estimated from."""

    b: float
    sigma: float
    n: int


def b_value_aki(
    mags: ArrayLike,
    mc: float,
    dm: float | None = None,
    shi_bolt_k: float | None = None,
) -> BValue:
    """
    b-value by the Aki (1965) maximum-likelihood estimator.

    b = 1 / (ln10 (mean(M) - (Mc - dm/2)))

    with the Shi and Bolt (1982) standard error

    sigma = k b^2 sqrt(sum (M - mean(M))^2 / (n (n - 1))).

    Events below Mc - dm/2 are discarded. The half-bin offset is the lower edge of the
    completeness bin: a magnitude reported as Mc stands for the interval Mc +/- dm/2,
    so the smallest complete magnitude is Mc - dm/2 rather than Mc.

    That offset is NOT in Aki (1965), whose estimator is b = 1/(ln10 (mean(M) - Mc)).
    It is the correction for binned magnitudes attributed to Utsu (1965), as in the
    documentation of seismostats, and it is not cosmetic: on the reference catalogue
    at Mc 3.5 it moves b from 0.935 to 0.844, about five standard errors. A referee
    checking this estimator against Aki alone will not find the term, so both papers
    are cited. Some papers cite Utsu (1966) for the size of the bias that binning
    causes; the attribution has not been checked against the 1965 paper itself.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    mc : float
        Completeness magnitude.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM` (0.1). Pass 0.0 for unbinned
        magnitudes, where no half-bin offset applies.
    shi_bolt_k : float, optional
        Coefficient of the standard error. Defaults to `constants.SHI_BOLT_K` (2.30),
        the published rounding of ln(10).

    Returns
    -------
    BValue
        Named tuple of b, sigma and n.

    Raises
    ------
    ValueError
        If fewer than two events lie at or above the completeness threshold, where the
        standard error is undefined.

    References
    ----------
    Aki, K. (1965), for the maximum-likelihood estimator. Utsu, T. (1965), for the
    half-bin correction applied to binned magnitudes. Shi, Y. and Bolt, B. A.
    (1982), for the standard error.
    """
    dm = constants.DM if dm is None else dm
    shi_bolt_k = constants.SHI_BOLT_K if shi_bolt_k is None else shi_bolt_k
    threshold = mc - dm / 2
    m = np.asarray(mags, float)
    m = m[at_or_above(m, threshold)]
    n = m.size
    if n < 2:
        raise ValueError(
            f"b-value needs at least two events at or above {threshold}; got {n}"
        )
    mean_m = m.mean()
    # b_value_tinti refuses this and b_value_aki did not: with every magnitude on
    # the threshold the denominator is zero and b came back as an infinity, or,
    # if a caller passes a threshold above the data, negative - a b-value that
    # says event numbers rise with magnitude.
    if mean_m <= threshold:
        raise ValueError(
            f"the mean magnitude {mean_m:.4f} does not exceed the completeness "
            f"threshold {threshold:.4f}, so b is not defined here; the sample is "
            f"either a single magnitude bin or was taken above the data"
        )
    b = 1.0 / (np.log(10) * (mean_m - threshold))
    sigma = shi_bolt_k * b**2 * np.sqrt(((m - mean_m) ** 2).sum() / (n * (n - 1)))
    return BValue(float(b), float(sigma), n)


class BStability(NamedTuple):
    """b and its standard error as a function of the threshold applied."""

    thresholds: NDArray[np.float64]
    b: NDArray[np.float64]
    sigma: NDArray[np.float64]
    n: NDArray[np.int64]


def b_stability(
    mags: ArrayLike,
    thresholds: ArrayLike | None = None,
    dm: float | None = None,
    min_events: int = 50,
) -> BStability:
    """
    b against the threshold it was estimated at.

    Above a correctly estimated completeness magnitude the b-value should not
    depend on where the threshold is put; a curve that keeps climbing is the
    standard sign that completeness has been placed too low, or that the sample
    is not a single Gutenberg-Richter population. The Shi and Bolt error is the
    scatter at one threshold and says nothing about this, so the curve carries
    information the headline figure cannot.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    thresholds : array_like, optional
        Thresholds to evaluate. Defaults to a grid of 2.5 magnitude units
        starting at the mode of the incremental distribution.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM`.
    min_events : int, optional
        Thresholds retaining fewer events than this are dropped, since b is not
        meaningful there. 50 by default.

    Returns
    -------
    BStability
        Arrays of threshold, b, sigma and n, one entry per usable threshold.

    References
    ----------
    Cao, A. and Gao, S. S. (2002). Woessner, J. and Wiemer, S. (2005).
    """
    dm = constants.DM if dm is None else dm
    m = np.asarray(mags, float)
    if thresholds is None:
        # Anchored at the mode of the incremental distribution, not at the
        # smallest magnitude: completeness cannot lie below the mode, which is
        # the premise of the maximum-curvature method, and one anomalously small
        # event moves the minimum while leaving the mode where it was.
        edges, inc, _ = fmd(m, dm)
        start = float(edges[int(np.argmax(inc))])
        # Halves up, as everywhere else: 2.5/0.2 is exactly 12.5, and Python's
        # round sends that down where JavaScript's Math.round sends it up, so the
        # package and the browser page searched different numbers of thresholds.
        steps = int(np.floor(2.5 / dm + 0.5 + constants.GRID_TOLERANCE))
        thresholds = np.round(start + dm * np.arange(steps + 1), 10)
    thresholds = np.atleast_1d(np.asarray(thresholds, float))

    kept, bs, sigmas, ns = [], [], [], []
    for threshold in thresholds:
        if int(at_or_above(m, threshold - dm / 2).sum()) < min_events:
            continue
        try:
            estimate = b_value_aki(m, threshold, dm=dm)
        except ValueError:
            # No b exists at this threshold - every magnitude in the sample sits
            # on the completeness bin's lower edge, so the mean does not exceed
            # it. That is a point to leave out of the curve, not a reason to
            # abandon the curve; the browser page has always skipped it.
            continue
        kept.append(threshold)
        bs.append(estimate.b)
        sigmas.append(estimate.sigma)
        ns.append(estimate.n)
    return BStability(
        np.array(kept), np.array(bs), np.array(sigmas), np.array(ns, dtype=np.int64)
    )


def mc_b_stability(
    mags: ArrayLike,
    dm: float | None = None,
    span: float = 0.5,
    min_events: int = 50,
) -> float | None:
    """
    Completeness magnitude by the b-value stability method.

    The lowest threshold at which b has stopped changing: b there differs from
    the mean of b over the next `span` magnitude units by no more than its own
    standard error. Where maximum curvature answers "where is the peak of the
    incremental distribution", this answers "from where onward does the slope
    stop moving", and the two disagreeing is itself worth reporting.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM`.
    span : float, optional
        Width of the averaging window above each candidate, 0.5 by default.
    min_events : int, optional
        Minimum events above a threshold for it to be a candidate.

    Returns
    -------
    float or None
        The completeness magnitude, or None if b never stabilises over the
        range available.

    References
    ----------
    Cao, A. and Gao, S. S. (2002). Woessner, J. and Wiemer, S. (2005).
    """
    dm = constants.DM if dm is None else dm
    # How many bins the averaging window spans, rounded the same way magnitudes
    # are binned: halves up, on the grid. Python's round() is banker's rounding
    # and sends 2.5 down to 2 where JavaScript's Math.round sends it up to 3, so
    # at dm 0.2 the package and the browser page averaged over different windows
    # and reported different completeness magnitudes, 4.0 against 4.2, from
    # identical b-stability curves. Where span is not a whole number of bins the
    # window is the next bin up: 0.5 at dm 0.2 averages over 0.6.
    steps = int(np.floor(span / dm + 0.5 + constants.GRID_TOLERANCE))
    if steps < 1:
        raise ValueError(
            f"a bin width of {dm} rounds the averaging window of {span} down "
            f"to no bins at all, so there is nothing to average over and "
            f"stability cannot be tested; use a narrower dm or a wider span"
        )
    curve = b_stability(mags, dm=dm, min_events=min_events)
    for i in range(len(curve.thresholds) - steps):
        window = curve.b[i : i + steps + 1]
        if abs(window.mean() - curve.b[i]) <= curve.sigma[i]:
            return float(curve.thresholds[i])
    return None


class BShift(NamedTuple):
    """The change in b between two nested thresholds, with its own sampling error."""

    mc_low: float
    mc_high: float
    b_low: float
    b_high: float
    shift: float
    se_constant_b: float
    z_constant_b: float
    p_constant_b: float
    se_boot: float | None
    ci_boot: tuple[float, float] | None
    n_low: int
    n_high: int
    n_boot: int


def b_shift(
    mags: ArrayLike,
    mc_low: float,
    mc_high: float,
    dm: float | None = None,
    n_boot: int = 2000,
    confidence: float = 0.95,
    seed: int = 0,
    shi_bolt_k: float | None = None,
) -> BShift:
    """
    How far b moves between two thresholds, judged against the error of that move.

    The standard error quoted beside a b-value describes that one estimate. The
    shift b(Mc_high) - b(Mc_low) has an error of its own, and neither threshold's
    Shi and Bolt error is it: the sample above the higher threshold is a subset of
    the sample above the lower one, so the two estimates are correlated. Dividing
    a shift by the error of one end gives a ratio that has no probabilistic meaning
    and that can be wrong by a factor of several in either direction.

    Two errors are returned, because they answer different questions.

    ``se_constant_b`` is the error of the shift if b were the same at both
    thresholds, the hypothesis of a stable b. The estimator at the lower threshold
    uses more events and is efficient under that hypothesis, so the two estimates
    differ by noise that is uncorrelated with the better one (Hausman 1978) and

        var(b_high - b_low) = var(b_high) - var(b_low)
                            = sigma_low^2 (n_low / n_high - 1)

    where sigma_low is the Shi and Bolt error at the lower threshold. The ratio of
    the shift to this error, ``z_constant_b``, is a test of that hypothesis.

    ``se_boot`` and ``ci_boot`` come from resampling the events above the lower
    threshold with replacement and recomputing both estimates on every resample.
    They hold whether or not b is constant, which is what a confidence interval for
    the size of a shift that is not zero needs, and they make no use of the Shi and
    Bolt formula.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    mc_low, mc_high : float
        The two thresholds, ``mc_low < mc_high``.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM`.
    n_boot : int, optional
        Bootstrap resamples, 2000 by default; an interval needs far more than the
        200 that suffice for a standard error. Pass 0 to skip the bootstrap, in which
        case ``se_boot`` and ``ci_boot`` are None.
    confidence : float, optional
        Level of the percentile interval, 0.95 by default.
    seed : int, optional
        Seed of the resampling generator, so a reported interval is reproducible.
    shi_bolt_k : float, optional
        Coefficient of the standard error. Defaults to `constants.SHI_BOLT_K`.

    Returns
    -------
    BShift
        Thresholds, b at each, the shift, ``se_constant_b``, ``z_constant_b`` and its
        two-sided normal p-value ``p_constant_b``, the bootstrap error and interval,
        the event count at each threshold, and the number of resamples used.

    Raises
    ------
    ValueError
        If the thresholds are not in increasing order, if either leaves fewer than two
        events, if both select the same events, if ``n_boot`` is 1 or negative, or if
        ``confidence`` is not strictly between 0 and 1.

    References
    ----------
    Hausman, J. A. (1978). Efron, B. (1979). Shi, Y. and Bolt, B. A. (1982).
    """
    dm = constants.DM if dm is None else dm
    if not mc_low < mc_high:
        raise ValueError(
            f"mc_low must be below mc_high; got {mc_low} and {mc_high}. The shift is "
            f"taken from the lower threshold to the higher one"
        )
    if n_boot < 0:
        raise ValueError(f"n_boot cannot be negative; got {n_boot}")
    if n_boot == 1:
        raise ValueError(
            "n_boot of 1 cannot give a spread. Use 0 to skip the bootstrap, or a "
            "number of resamples large enough for an interval"
        )
    if not 0.0 < confidence < 1.0:
        raise ValueError(
            f"confidence must lie strictly between 0 and 1; got {confidence}"
        )

    low = b_value_aki(mags, mc_low, dm=dm, shi_bolt_k=shi_bolt_k)
    high = b_value_aki(mags, mc_high, dm=dm, shi_bolt_k=shi_bolt_k)
    if low.n == high.n:
        raise ValueError(
            f"the thresholds {mc_low} and {mc_high} select the same {low.n} events, "
            f"so there is no shift to measure; widen the gap between them"
        )
    shift = high.b - low.b
    se = low.sigma * float(np.sqrt(low.n / high.n - 1.0))
    z = shift / se
    p_value = float(2.0 * norm.sf(abs(z)))

    se_boot = None
    interval = None
    used = 0
    if n_boot:
        floor_low = mc_low - dm / 2
        floor_high = mc_high - dm / 2
        sample = np.asarray(mags, float)
        sample = sample[at_or_above(sample, floor_low)]
        rng = np.random.default_rng(seed)
        shifts = np.full(n_boot, np.nan)
        for i in range(n_boot):
            resample = sample[rng.integers(0, sample.size, sample.size)]
            above = resample[at_or_above(resample, floor_high)]
            if above.size < 2:
                continue
            gap_low = resample.mean() - floor_low
            gap_high = above.mean() - floor_high
            if gap_low <= 0 or gap_high <= 0:
                continue
            shifts[i] = (1.0 / gap_high - 1.0 / gap_low) / np.log(10)
        shifts = shifts[np.isfinite(shifts)]
        used = int(shifts.size)
        if used < 2:
            raise ValueError(
                "fewer than two bootstrap resamples gave a usable shift; the sample "
                "above the higher threshold is too small to resample"
            )
        tail = 100.0 * (1.0 - confidence) / 2.0
        lo, hi = np.percentile(shifts, [tail, 100.0 - tail])
        se_boot = float(np.std(shifts, ddof=1))
        interval = (float(lo), float(hi))

    return BShift(
        float(mc_low),
        float(mc_high),
        low.b,
        high.b,
        float(shift),
        float(se),
        float(z),
        p_value,
        se_boot,
        interval,
        low.n,
        high.n,
        used,
    )


# Chi-squared tail probability below which a candidate threshold is rejected without
# simulating its null; see `b_plateau`.
_OBVIOUS_REJECTION = 1e-10


class BPlateau(NamedTuple):
    """Where b stops depending on the threshold, by a test rather than by eye."""

    onset: float | None
    b: float | None
    sigma: float | None
    n: int | None
    statistic: float | None
    dof: int | None
    p_value: float | None
    alpha: float
    method: str
    thresholds: NDArray[np.float64]
    statistic_profile: NDArray[np.float64]
    dof_profile: NDArray[np.int64]
    p_profile: NDArray[np.float64]


def b_plateau(
    mags: ArrayLike,
    thresholds: ArrayLike | None = None,
    dm: float | None = None,
    min_events: int = 50,
    alpha: float = 0.05,
    min_window: int = 3,
    n_simulations: int = 500,
    seed: int = 0,
    shi_bolt_k: float | None = None,
) -> BPlateau:
    """
    The lowest threshold above which b shows no detectable dependence on the threshold.

    `mc_b_stability` calls b stable when it sits within one of its own standard
    errors of its mean over the next half magnitude unit. That rule has no level
    attached to it. This one tests, at a stated level, for each candidate threshold
    the hypothesis that b is the same at every higher threshold in the usable range,
    and the onset is the lowest candidate at which the hypothesis is not rejected.
    The level applies to each candidate; the rule that takes the lowest candidate
    that passes has no single error rate of its own.

    The statistic is built from the Hausman (1978) comparison of nested
    maximum-likelihood estimates. Under a constant b the estimate at the candidate
    threshold, which uses the most events, is efficient, and the change in b from one
    threshold to the next higher one is uncorrelated with the change before it. Each
    step therefore contributes an independent term, and over the steps from the
    candidate upwards

        T = sum over steps of (b_j - b_(j-1))^2 / (sigma_i^2 n_i (1/n_j - 1/n_(j-1)))

    where i is the candidate threshold, n counts the events at each threshold and
    sigma_i is the Shi and Bolt error at the candidate. Over a single step T is the
    square of `b_shift`'s ``z_constant_b``.

    T is referred to a chi-squared distribution with one degree of freedom per step
    only when ``n_simulations`` is 0, and that reference is approximate. The steps near
    the top of the range rest on a few dozen events, their differences are skewed, and
    the sum of their squares has a heavier upper tail than chi-squared: on 2,000
    catalogues of 8,500 events simulated with a constant b the approximate test
    rejected 10.1 per cent at the 5 per cent level (Wilson 95 per cent interval 8.8
    to 11.4; ``examples/synthetic_incompleteness.py`` measures it). This is the
    same situation as the decay-fit test in `omori_fit_test`, and it is handled
    the same way. By default the null distribution is simulated: events are drawn
    above the candidate from the exponential law the candidate's own b implies, the
    whole curve is recomputed on each draw, and the p-value is the share of draws
    whose T is at least as large as the observed one.

    What the onset is and is not. It is a decision rule, not an estimate with a
    confidence level: the lowest candidate that passes is chosen, so a threshold can
    pass by luck where the test has little power, which is why a candidate needs at
    least ``min_window`` steps above it. A threshold that passes has no detectable
    dependence of b on the threshold above it; that is weaker than a proof that
    completeness holds there, and a b that drifts slowly over a range with few events
    will pass. Read the profile as well as the onset.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    thresholds : array_like, optional
        Thresholds to test, as for `b_stability`.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM`.
    min_events : int, optional
        Thresholds retaining fewer events than this are dropped, 50 by default.
    alpha : float, optional
        Level at which constancy is rejected, 0.05 by default.
    min_window : int, optional
        Fewest steps above a candidate for it to be tested, 3 by default, so a
        plateau needs at least four thresholds. A candidate near the top of the range
        has almost nothing above it to disagree with.
    n_simulations : int, optional
        Draws used to simulate the null distribution of T at each candidate, 500 by
        default. Zero selects the chi-squared reference, which is approximate and
        reported as such. The cost grows with the number of events and thresholds; a
        catalogue of 30,000 events takes about a quarter of a minute.
    seed : int, optional
        Seed of the simulation, so a reported p-value is reproducible.
    shi_bolt_k : float, optional
        Coefficient of the standard error. Defaults to `constants.SHI_BOLT_K`.

    Returns
    -------
    BPlateau
        The onset threshold, and b, its error and n there; T, the degrees of freedom
        and the p-value at the onset; ``alpha``; how the p-value was obtained; and the
        whole profile, one entry per tested threshold, with NaN or 0 where a candidate
        has too few steps above it. The onset and everything attached to it are None
        where no candidate passes. With a simulated null the p-value of a candidate is
        (1 + draws with T at least as large) / (1 + draws), so it cannot fall below
        1 / (1 + n_simulations).

    Raises
    ------
    ValueError
        If ``alpha`` is not strictly between 0 and 1, ``min_window`` is below 1, or
        ``n_simulations`` is negative.

    References
    ----------
    Hausman, J. A. (1978). Cao, A. and Gao, S. S. (2002). Woessner, J. and Wiemer,
    S. (2005). Ogata, Y. (1988) for the simulated null.
    """
    dm = constants.DM if dm is None else dm
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must lie strictly between 0 and 1; got {alpha}")
    if min_window < 1:
        raise ValueError(f"min_window must be at least 1; got {min_window}")
    if n_simulations < 0:
        raise ValueError(f"n_simulations cannot be negative; got {n_simulations}")
    coefficient = constants.SHI_BOLT_K if shi_bolt_k is None else shi_bolt_k
    curve = b_stability(mags, thresholds=thresholds, dm=dm, min_events=min_events)
    t, b, _, n = curve
    # b_stability used the published coefficient; the variance scale needs the one
    # asked for, and sigma is linear in it.
    sigma = curve.sigma * coefficient / constants.SHI_BOLT_K
    count = t.size
    statistic = np.full(count, np.nan)
    dof = np.zeros(count, dtype=np.int64)
    p = np.full(count, np.nan)
    if count >= 2:
        step_b = np.diff(b)
        # Variance of each step in units of sigma_i^2 n_i. Two thresholds that keep the
        # same events differ by nothing and carry no information; they are skipped.
        step_v = 1.0 / n[1:] - 1.0 / n[:-1]
        informative = step_v > 0.0
        terms = np.where(
            informative, step_b**2 / np.where(informative, step_v, 1.0), 0.0
        )
        tail_terms = np.cumsum(terms[::-1])[::-1]
        tail_steps = np.cumsum(informative[::-1].astype(np.int64))[::-1]
        for i in range(count - 1):
            steps = int(tail_steps[i])
            if steps < 1:
                continue
            statistic[i] = tail_terms[i] / (sigma[i] ** 2 * n[i])
            dof[i] = steps

    tested = [i for i in range(count) if dof[i] >= min_window]
    if n_simulations:
        m = np.asarray(mags, float)
        rng = np.random.default_rng(seed)
        for i in tested:
            # A statistic this far into the chi-squared tail is beyond anything the
            # heavier simulated tail reaches too, and a catalogue with a long climb
            # in b would otherwise spend its whole run simulating candidates that
            # fail by hundreds. Those get the smallest p-value a simulation returns.
            if chi2_distribution.sf(statistic[i], dof[i]) < _OBVIOUS_REJECTION:
                p[i] = 1.0 / (1.0 + n_simulations)
                continue
            p[i] = _simulated_plateau_p(
                m, t[i:], dm, n[i:], statistic[i], n_simulations, coefficient, rng
            )
        method = f"simulated null, {n_simulations} draws"
    else:
        for i in tested:
            p[i] = float(chi2_distribution.sf(statistic[i], dof[i]))
        method = "chi-squared, approximate"

    onset = next((i for i in tested if p[i] >= alpha), None)
    if onset is None:
        return BPlateau(
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            alpha,
            method,
            t,
            statistic,
            dof,
            p,
        )
    return BPlateau(
        float(t[onset]),
        float(b[onset]),
        float(sigma[onset]),
        int(n[onset]),
        float(statistic[onset]),
        int(dof[onset]),
        float(p[onset]),
        alpha,
        method,
        t,
        statistic,
        dof,
        p,
    )


def _simulated_plateau_p(
    mags: NDArray[np.float64],
    window: NDArray[np.float64],
    dm: float,
    observed_counts: NDArray[np.int64],
    observed: float,
    n_simulations: int,
    coefficient: float,
    rng: np.random.Generator,
) -> float:
    """Monte Carlo p-value of T for the window of thresholds starting at its first.

    Events are drawn above the first threshold's lower edge from the exponential law
    whose mean excess is the one observed there, as many as were observed, and the
    curve over the window is recomputed on each draw. The thresholds are those of the
    observed curve, so a threshold that kept at least ``min_events`` events in the data
    is evaluated on every draw whatever it keeps there. A draw leaving fewer than two
    events at a threshold cannot give a b and is not counted.
    """
    floors = window - dm / 2
    n0 = int(observed_counts[0])
    sample = mags[at_or_above(mags, floors[0])]
    mean_excess = float(sample.mean() - floors[0])
    exceeded = 0
    used = 0
    for _ in range(n_simulations):
        x = floors[0] + mean_excess * rng.exponential(size=n0)
        # Which interval between consecutive floors each event falls in; summing the
        # intervals from the top gives the count and total above every floor.
        slot = np.searchsorted(floors, x, side="right") - 1
        counts = np.bincount(slot, minlength=floors.size)
        totals = np.bincount(slot, weights=x, minlength=floors.size)
        n_above = np.cumsum(counts[::-1])[::-1]
        if n_above[-1] < 2:
            continue
        gap = np.cumsum(totals[::-1])[::-1] / n_above - floors
        if np.any(gap <= 0):
            continue
        b = 1.0 / (np.log(10) * gap)
        sigma_sq_n = (coefficient * b[0] ** 2) ** 2 * x.var(ddof=1)
        step_v = 1.0 / n_above[1:] - 1.0 / n_above[:-1]
        informative = step_v > 0.0
        terms = np.where(
            informative, np.diff(b) ** 2 / np.where(informative, step_v, 1.0), 0.0
        )
        used += 1
        if terms.sum() / sigma_sq_n >= observed:
            exceeded += 1
    if used == 0:
        return float("nan")
    return (1.0 + exceeded) / (1.0 + used)


def b_value_tinti(
    mags: ArrayLike,
    mc: float,
    dm: float | None = None,
    shi_bolt_k: float | None = None,
) -> BValue:
    """
    b-value by the exact maximum-likelihood estimator for binned magnitudes.

    b = ln(1 + dm / mean(M - Mc)) / (dm ln10)

    Where `b_value_aki` applies Utsu's half-bin offset to the Aki estimator,
    which is a first-order approximation, this is the exact solution for
    magnitudes reported on a grid of width dm. The two converge as dm shrinks:
    on the reference catalogue they differ by 0.32 per cent at dm 0.1 and by
    0.004 per cent at dm 0.01.

    This is the estimator the independent package `seismostats` uses, and
    `examples/compare_with_seismostats.py` checks that this implementation
    reproduces it. The package's own reported values use `b_value_aki`, which is
    the convention of the thesis and of the spreadsheet implementation; this
    function exists so the choice can be tested rather than assumed.

    Parameters
    ----------
    mags : array_like
        Event magnitudes.
    mc : float
        Completeness magnitude.
    dm : float, optional
        Magnitude bin width. Defaults to `constants.DM` (0.1).
    shi_bolt_k : float, optional
        Coefficient of the standard error. Defaults to `constants.SHI_BOLT_K`.

    Returns
    -------
    BValue
        Named tuple of b, sigma and n.

    Raises
    ------
    ValueError
        If fewer than two events lie at or above the completeness threshold.

    References
    ----------
    Tinti, S. and Mulargia, F. (1987). Shi, Y. and Bolt, B. A. (1982).
    """
    dm = constants.DM if dm is None else dm
    shi_bolt_k = constants.SHI_BOLT_K if shi_bolt_k is None else shi_bolt_k
    threshold = mc - dm / 2
    m = np.asarray(mags, float)
    m = m[at_or_above(m, threshold)]
    n = m.size
    if n < 2:
        raise ValueError(
            f"b-value needs at least two events at or above {threshold}; got {n}"
        )
    mean_above = (m - mc).mean()
    if mean_above <= 0:
        raise ValueError("mean magnitude does not exceed the completeness magnitude")
    b = np.log1p(dm / mean_above) / (dm * np.log(10))
    mean_m = m.mean()
    sigma = shi_bolt_k * b**2 * np.sqrt(((m - mean_m) ** 2).sum() / (n * (n - 1)))
    return BValue(float(b), float(sigma), n)
