"""Regression against the locked reference values for the 2023 Kahramanmaras sequence.

These values were produced by the validated Python/SciPy reference implementation and
reproduced independently by the Google Apps Script implementation, whose optimiser is
a hand-written Nelder-Mead simplex rather than SciPy's. Two independent optimisers
agreeing to three decimals on a real sequence is what makes them reference values.

If a change to this package makes this file fail, the change is wrong, not the numbers.

Locked: 3469 events in a 180-day window; Mc 3.4; b 0.844 +/- 0.019 on 1529 events at or
above M 3.5; Omori p 1.161, c 0.497, k 359. At full precision the fit gives
p 1.160897, c 0.496822, k 358.868.
"""

from pathlib import Path

import pandas as pd
import pytest

from tremor_lab import analyze_case, b_value_aki, fit_omori, mc_maxcurvature

DATA = Path(__file__).parent / "data"

MAINSHOCK = {"t": "2023-02-06 01:17:32", "lat": 37.1578, "lon": 36.8092, "mw": 7.8}
THRESHOLD = 3.5


@pytest.fixture(scope="module")
def catalog():
    return pd.read_csv(DATA / "kahramanmaras_180d.csv")


@pytest.fixture(scope="module")
def result(catalog):
    return analyze_case(
        catalog, MAINSHOCK, mc_threshold=THRESHOLD, window_days=180, n_boot=0
    )


def test_the_fixture_is_the_catalogue_the_values_were_locked_against(catalog):
    assert len(catalog) == 3469
    assert list(catalog.columns) == ["dt_days", "mw", "lat", "lon"]


def test_events_in_window(result):
    assert result["n_events"] == 3469


def test_completeness_magnitude(catalog):
    assert mc_maxcurvature(catalog["mw"].to_numpy()) == 3.4


def test_completeness_magnitude_through_the_driver(result):
    assert result["mc"] == 3.4


def test_b_value(catalog):
    mags = catalog.loc[catalog["mw"] >= THRESHOLD, "mw"].to_numpy()
    b = b_value_aki(mags, THRESHOLD)
    assert b.n == 1529
    assert b.b == pytest.approx(0.844, abs=0.0005)
    assert b.sigma == pytest.approx(0.019, abs=0.0005)


def test_b_value_through_the_driver(result):
    assert result["n_above"] == 1529
    assert result["b_value"].n == 1529
    assert result["b_value"].b == pytest.approx(0.844, abs=0.0005)
    assert result["b_value"].sigma == pytest.approx(0.019, abs=0.0005)


def test_omori_decay(catalog):
    times = catalog.loc[catalog["mw"] >= THRESHOLD, "dt_days"].to_numpy()
    fit = fit_omori(times)
    assert fit.n == 1529
    assert fit.p == pytest.approx(1.161, abs=0.001)
    assert fit.c == pytest.approx(0.497, abs=0.001)
    assert fit.k == pytest.approx(359, abs=1.0)


def test_omori_decay_through_the_driver(result):
    assert result["omori"].n == 1529
    assert result["omori"].p == pytest.approx(1.161, abs=0.001)
    assert result["omori"].c == pytest.approx(0.497, abs=0.001)
    assert result["omori"].k == pytest.approx(359, abs=1.0)


def test_the_observation_interval_is_the_last_event_not_the_window_length(catalog):
    # The locked fit ran to the last event above threshold, at 179.596 days. Fitting
    # to a nominal 180 days would move p in the third decimal.
    times = catalog.loc[catalog["mw"] >= THRESHOLD, "dt_days"].to_numpy()
    assert times.max() == pytest.approx(179.596, abs=0.001)
    assert fit_omori(times, t_end=180.0).p != fit_omori(times).p
