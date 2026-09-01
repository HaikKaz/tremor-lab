"""Aftershock-sequence statistics: completeness, b-value, Omori-Utsu decay, energy.

Published constants live in `tremor_lab.constants` and can be overridden per call or
for a whole session; see that module for how.
"""

from tremor_lab import constants
from tremor_lab.analysis import analyze_case
from tremor_lab.bvalue import BValue, b_value_aki
from tremor_lab.catalog import CATALOG_COLUMNS, haversine_km, read_catalog, window
from tremor_lab.completeness import FMD, fmd, mc_maxcurvature
from tremor_lab.magnitude import (
    bath_mag,
    bath_ratio,
    energy_joules,
    mb_to_mw,
    ms_to_mw,
    to_mw,
)
from tremor_lab.omori import (
    Omori,
    OmoriBootstrap,
    bootstrap_omori,
    fit_omori,
    omori_nll,
)

__version__ = "1.0.0"

__all__ = [
    "CATALOG_COLUMNS",
    "FMD",
    "BValue",
    "Omori",
    "OmoriBootstrap",
    "analyze_case",
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
    "mc_maxcurvature",
    "ms_to_mw",
    "omori_nll",
    "read_catalog",
    "to_mw",
    "window",
]
