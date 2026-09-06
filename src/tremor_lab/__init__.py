"""Aftershock-sequence statistics: completeness, b-value, Omori-Utsu decay, energy.

Published constants live in `tremor_lab.constants` and can be overridden per call or
for a whole session; see that module for how.
"""

from tremor_lab import constants
from tremor_lab.analysis import analyze_case
from tremor_lab.bvalue import (
    BStability,
    BValue,
    b_stability,
    b_value_aki,
    mc_b_stability,
)
from tremor_lab.catalog import CATALOG_COLUMNS, haversine_km, read_catalog, window
from tremor_lab.completeness import (
    FMD,
    GoodnessOfFit,
    fmd,
    mc_goodness_of_fit,
    mc_maxcurvature,
)
from tremor_lab.magnitude import (
    bath_mag,
    bath_ratio,
    energy_joules,
    mb_to_mw,
    ms_to_mw,
    to_mw,
)
from tremor_lab.omori import (
    FitTest,
    Omori,
    OmoriBootstrap,
    bootstrap_omori,
    fit_omori,
    omori_fit_test,
    omori_nll,
)

__version__ = "1.0.0"

__all__ = [
    "CATALOG_COLUMNS",
    "FMD",
    "BStability",
    "BValue",
    "FitTest",
    "GoodnessOfFit",
    "Omori",
    "OmoriBootstrap",
    "analyze_case",
    "b_stability",
    "b_value_aki",
    "bath_mag",
    "bath_ratio",
    "bootstrap_omori",
    "constants",
    "energy_joules",
    "fit_omori",
    "fmd",
    "haversine_km",
    "mb_to_mw",
    "mc_b_stability",
    "mc_goodness_of_fit",
    "mc_maxcurvature",
    "ms_to_mw",
    "omori_fit_test",
    "omori_nll",
    "read_catalog",
    "to_mw",
    "window",
]
