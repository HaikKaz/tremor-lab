"""Published constants used by the estimators, all of them overridable.

Every value here is a default, not a fixed law of the package.

Most can be changed for a single call, by passing the matching keyword argument::

    energy_joules(7.8, a=1.44, b=5.24)
    mc_maxcurvature(mags, dm=0.05, correction=0.0)

Change it for a whole session, by reassigning the name here. Functions read these at
call time, so the new value applies to every later call, including calls made inside
`analyze_case`::

    from tremor_lab import constants
    constants.DELTA_MB = 1.2

The command-line tool exposes the same names in the ``[constants]`` section of its
settings file, so they can be changed without writing Python.

The three fit-quality bounds, ``OMORI_C_FLOOR``, ``OMORI_P_MIN`` and ``OMORI_P_MAX``,
are the exception: they have no keyword argument, because they judge a fit rather than
enter it. They are changed by reassignment here or from the settings file.
"""

# Radiated energy: log10(E) = ENERGY_A * M + ENERGY_B, with E in joules.
# Gutenberg and Richter (1956); Kanamori (1977).
ENERGY_A = 1.5
ENERGY_B = 4.8

# Bath magnitude deficit between a mainshock and its largest aftershock.
# Central value of the reported 1.1-1.2 range. Bath (1965).
DELTA_MB = 1.15

# Scordilis (2006) global Ms-to-Mw relation, in two branches meeting at MS_MW_BRANCH.
MS_MW_BRANCH = 6.1
MS_MW_LOW_SLOPE = 0.67
MS_MW_LOW_INTERCEPT = 2.07
MS_MW_HIGH_SLOPE = 0.99
MS_MW_HIGH_INTERCEPT = 0.08

# Scordilis (2006) global mb-to-Mw relation.
MB_MW_SLOPE = 0.85
MB_MW_INTERCEPT = 1.03

# Sphere radius used for great-circle distances.
EARTH_RADIUS_KM = 6371.0

# Magnitude bin width of the frequency-magnitude distribution.
DM = 0.1

# Additive correction applied to the maximum-curvature completeness magnitude.
# Wiemer and Wyss (2000).
MC_CORRECTION = 0.2

# Coefficient of the Shi and Bolt (1982) b-value standard error. It stands in for
# ln(10) = 2.3026 as published.
SHI_BOLT_K = 2.30

# Default aftershock window in days.
WINDOW_DAYS = 180.0

# Starting point of the Omori-Utsu search, as (c, p) in days and dimensionless.
OMORI_C0 = 0.5
OMORI_P0 = 1.1

# Bootstrap resamples used for the Omori uncertainty.
N_BOOT = 200

# Replicates used to calibrate the decay fit test. The parameters are estimated
# from the data being tested, so the textbook Kolmogorov p-value does not apply
# and the null distribution has to be simulated. 600 rather than 200 because the
# p-value is itself an estimate: at 200 its standard error is about a quarter of
# a p-value near 0.05, and two runs of the same test land on opposite sides of
# the threshold by luck. Costs about three seconds on a 1,500-event sequence.
N_FIT_SIMULATIONS = 600

# Bounds used only to judge whether a fitted decay is worth believing; they never
# constrain the fit itself. c below the floor means the offset has collapsed onto
# the boundary of the model, where the likelihood has no interior maximum and c is
# not identified. The p range brackets the values compiled from real sequences by
# Utsu, Ogata and Matsu'ura (1995); outside it the fit is reported with a caution.
OMORI_C_FLOOR = 1e-3
OMORI_P_MIN = 0.5
OMORI_P_MAX = 2.0

# Minimum event counts below which an estimate is reported as unavailable rather than
# fitted to a catalogue too sparse to support it.
MIN_EVENTS_FOR_MC = 50
MIN_EVENTS_FOR_FIT = 100

# A snapshot taken here, before any caller or settings file can reassign a
# name, so a report can say which values were changed from the published ones.
_PUBLISHED = {
    name: value
    for name, value in list(globals().items())
    if name.isupper() and not name.startswith("_")
}
