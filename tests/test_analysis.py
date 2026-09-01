"""The per-case driver: what it reports, and what it refuses to report."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tremor_lab import constants
from tremor_lab.analysis import analyze_case
from tremor_lab.catalog import CATALOG_COLUMNS, read_catalog

DATA = Path(__file__).parent / "data"

MAINSHOCK = {"t": "2023-02-06 01:17:32", "lat": 37.0, "lon": 37.0, "mw": 7.8}


@pytest.fixture(scope="module")
def kahramanmaras():
    return pd.read_csv(DATA / "kahramanmaras_180d.csv")


def sparse_catalog(n, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "dt_days": rng.uniform(0.01, 180.0, n),
            "mw": np.round(3.0 + rng.exponential(0.43, n), 1),
        }
    )


def test_the_reader_produces_exactly_what_the_analysis_consumes():
    # The contract between the two halves of the package: whatever read_catalog
    # returns must satisfy analyze_case without further preparation.
    columns = {
        "date": "Tarih",
        "time": "Saat",
        "lat": "Enlem",
        "lon": "Boylam",
        "mag": "Mag",
    }
    cat = read_catalog(DATA / "raw_catalog_koeri.csv", columns, MAINSHOCK)
    assert tuple(cat.columns) == CATALOG_COLUMNS
    result = analyze_case(cat, MAINSHOCK)
    assert result["n_events"] == 4


def test_a_catalogue_without_the_required_columns_is_rejected():
    with pytest.raises(KeyError, match="dt_days"):
        analyze_case(pd.DataFrame({"mw": [3.0, 4.0]}), MAINSHOCK)


def test_energy_and_bath_expectation_come_from_the_mainshock(kahramanmaras):
    result = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=3.5, n_boot=0)
    assert np.log10(result["mainshock_energy_j"]) == pytest.approx(16.5)
    assert result["bath_mag"] == pytest.approx(6.65)


def test_the_estimated_completeness_is_used_when_no_threshold_is_given(kahramanmaras):
    result = analyze_case(kahramanmaras, MAINSHOCK, n_boot=0)
    assert result["mc"] == 3.4
    assert result["threshold"] == 3.4


def test_a_fixed_threshold_overrides_the_estimate(kahramanmaras):
    result = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=4.0, n_boot=0)
    assert result["mc"] == 3.4
    assert result["threshold"] == 4.0
    assert result["n_above"] < 1529


def test_the_window_is_applied_to_the_supplied_catalogue(kahramanmaras):
    short = analyze_case(kahramanmaras, MAINSHOCK, window_days=30, n_boot=0)
    assert short["n_events"] < 3469


def test_a_sparse_catalogue_reports_no_completeness_magnitude():
    result = analyze_case(sparse_catalog(20), MAINSHOCK)
    assert result["mc"] is None
    assert result["b_value"] is None and result["omori"] is None
    assert "needs more than 50" in result["note"]


def test_a_catalogue_too_thin_above_the_threshold_reports_no_fit():
    result = analyze_case(sparse_catalog(80), MAINSHOCK, mc_threshold=3.0)
    assert result["mc"] is not None
    assert result["b_value"] is None and result["omori"] is None
    assert "a stable fit " in result["note"]


def test_the_density_rules_are_configurable():
    forced = analyze_case(
        sparse_catalog(80),
        MAINSHOCK,
        mc_threshold=3.0,
        min_events_for_fit=10,
        n_boot=0,
    )
    assert forced["b_value"] is not None and forced["omori"] is not None


def test_the_density_rules_follow_the_constants(monkeypatch):
    monkeypatch.setattr(constants, "MIN_EVENTS_FOR_MC", 500)
    assert analyze_case(sparse_catalog(80), MAINSHOCK)["mc"] is None


def test_the_bootstrap_can_be_switched_off(kahramanmaras):
    result = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=3.5, n_boot=0)
    assert result["omori_bootstrap"] is None


def test_the_bootstrap_is_reproducible_from_its_seed(kahramanmaras):
    first = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=3.5, n_boot=8, seed=3)
    again = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=3.5, n_boot=8, seed=3)
    assert first["omori_bootstrap"] == again["omori_bootstrap"]
    # Repeating one seed proves determinism but not that the seed is used at all;
    # a different seed must reach the resampler and give a different spread.
    other = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=3.5, n_boot=8, seed=9)
    assert other["omori_bootstrap"] != first["omori_bootstrap"]


def test_the_frequency_magnitude_table_covers_every_event(kahramanmaras):
    result = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=3.5, n_boot=0)
    assert result["fmd"].inc.sum() == result["n_events"]
    assert result["fmd"].cum[0] == result["n_events"]


def test_the_bin_width_reaches_the_completeness_estimate(kahramanmaras):
    coarse = analyze_case(kahramanmaras, MAINSHOCK, dm=0.5, n_boot=0)
    assert coarse["mc"] != 3.4


def test_a_sound_fit_carries_no_caution(kahramanmaras):
    result = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=3.5, n_boot=0)
    assert result["omori_warning"] is None


def test_a_fit_whose_offset_collapses_is_flagged(monkeypatch, kahramanmaras):
    # Raising the floor above the fitted c is the cheapest way to exercise the
    # boundary case deterministically; the real trigger is a catalogue whose
    # first events sit on the origin time.
    monkeypatch.setattr(constants, "OMORI_C_FLOOR", 1.0)
    result = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=3.5, n_boot=0)
    assert "collapsed" in result["omori_warning"]


def test_a_decay_exponent_outside_the_reported_range_is_flagged(monkeypatch, kahramanmaras):
    monkeypatch.setattr(constants, "OMORI_P_MIN", 1.5)
    result = analyze_case(kahramanmaras, MAINSHOCK, mc_threshold=3.5, n_boot=0)
    assert "outside 1.5" in result["omori_warning"]


def test_rows_with_a_missing_magnitude_are_dropped_and_counted(kahramanmaras):
    holed = kahramanmaras.copy()
    holed.loc[10, "mw"] = np.nan
    holed.loc[20, "dt_days"] = np.nan
    result = analyze_case(holed, MAINSHOCK, mc_threshold=3.5, n_boot=0)
    assert result["n_unusable"] == 2
    assert result["n_events"] == 3467
    assert result["mc"] == 3.4
