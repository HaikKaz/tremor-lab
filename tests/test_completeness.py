"""Frequency-magnitude binning and the maximum-curvature completeness magnitude."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tremor_lab import constants
from tremor_lab.completeness import fmd, mc_goodness_of_fit, mc_maxcurvature

DATA = Path(__file__).parent / "data"


def test_incremental_counts_land_in_the_expected_bins():
    edges, inc, _ = fmd([3.0, 3.1, 3.1, 3.2])
    assert edges == pytest.approx([3.0, 3.1, 3.2])
    assert inc.tolist() == [1, 2, 1]


def test_cumulative_counts_events_at_or_above_each_bin():
    _, _, cum = fmd([3.0, 3.1, 3.1, 3.2])
    assert cum.tolist() == [4, 3, 1]


def test_every_event_is_counted_exactly_once():
    mags = np.round(np.random.default_rng(0).uniform(3.0, 7.0, 500), 1)
    _, inc, _ = fmd(mags)
    assert inc.sum() == len(mags)


def test_magnitudes_on_the_decimal_grid_land_on_their_own_bin():
    # The reason for binning by nearest centre: 3.2 has no exact binary form, so a
    # half-open interval can place it in the neighbouring bin.
    mags = np.round(np.arange(3.0, 7.05, 0.1), 1)
    edges, inc, _ = fmd(mags)
    assert inc.tolist() == [1] * len(mags)
    assert edges == pytest.approx(mags)


def test_bin_width_is_configurable():
    edges, inc, _ = fmd([3.0, 3.05, 3.1], dm=0.05)
    assert edges == pytest.approx([3.0, 3.05, 3.1])
    assert inc.tolist() == [1, 1, 1]


def test_default_bin_width_follows_the_constant(monkeypatch):
    monkeypatch.setattr(constants, "DM", 1.0)
    edges, inc, _ = fmd([3.0, 3.1, 4.0])
    assert edges == pytest.approx([3.0, 4.0])
    assert inc.tolist() == [2, 1]


def test_completeness_is_the_modal_magnitude_plus_the_correction():
    assert mc_maxcurvature([3.0, 3.1, 3.1, 3.2]) == pytest.approx(3.3)


def test_the_correction_can_be_switched_off():
    assert mc_maxcurvature([3.0, 3.1, 3.1, 3.2], correction=0.0) == pytest.approx(3.1)


def test_default_correction_follows_the_constant(monkeypatch):
    monkeypatch.setattr(constants, "MC_CORRECTION", 0.0)
    assert mc_maxcurvature([3.0, 3.1, 3.1, 3.2]) == pytest.approx(3.1)


def test_completeness_is_returned_on_the_bin_grid_without_floating_point_dust():
    mc = mc_maxcurvature([3.0, 3.2, 3.2, 3.5])
    assert mc == 3.4
    assert repr(mc) == "3.4"


def test_the_two_binning_conventions_agree_on_the_reference_catalogue():
    # thesis_analysis.py binned with np.histogram over half-open intervals; this
    # package assigns each magnitude to its nearest bin centre. The change is
    # numerical safety, not a change of method, so both must give the same answer.
    mags = pd.read_csv(DATA / "kahramanmaras_180d.csv")["mw"].to_numpy()
    bins = np.arange(np.floor(mags.min()), np.ceil(mags.max()) + 0.1, 0.1)
    counts, edges = np.histogram(mags, bins=bins)
    legacy_mc = round(float(edges[np.argmax(counts)] + 0.2), 2)
    assert mc_maxcurvature(mags) == legacy_mc == 3.4


def test_goodness_of_fit_on_the_reference_catalogue():
    mags = pd.read_csv(DATA / "kahramanmaras_180d.csv")["mw"].to_numpy()
    result = mc_goodness_of_fit(mags)
    assert result.mc == pytest.approx(3.0)
    assert result.r_value == pytest.approx(90.83, abs=0.05)
    # R is 100 (1 - sum|observed - synthetic| / sum observed), a residual
    # fraction, so it can never exceed 100. Dropping the absolute value would.
    assert np.all(result.r_values <= 100.0)
    assert result.candidates.size == result.r_values.size


def test_goodness_of_fit_returns_nothing_when_the_level_is_unreachable():
    mags = pd.read_csv(DATA / "kahramanmaras_180d.csv")["mw"].to_numpy()
    result = mc_goodness_of_fit(mags, confidence=99.9)
    assert result.mc is None
    assert result.r_value is None
    assert result.r_values.size > 0


def test_goodness_of_fit_prefers_the_lowest_candidate_that_reaches_the_level():
    mags = pd.read_csv(DATA / "kahramanmaras_180d.csv")["mw"].to_numpy()
    strict = mc_goodness_of_fit(mags, confidence=95.0)
    assert strict.mc is not None
    assert strict.mc > mc_goodness_of_fit(mags, confidence=90.0).mc
    assert strict.r_value >= 95.0


def test_one_anomalously_small_event_does_not_move_the_candidate_list():
    # A fixed number of candidates counted up from the catalogue minimum would be
    # pushed off the data entirely by a single stray event, and the method would
    # return nothing for a catalogue it had just handled.
    mags = pd.read_csv(DATA / "kahramanmaras_180d.csv")["mw"].to_numpy()
    clean = mc_goodness_of_fit(mags)
    with_outlier = mc_goodness_of_fit(np.append(mags, 0.0))
    assert with_outlier.mc == clean.mc
    assert with_outlier.r_value == pytest.approx(clean.r_value, abs=0.01)


def test_a_tie_between_two_modal_bins_is_broken_by_taking_the_lower():
    """Two equally populated bins is a real case, and it needs a stated rule.

    The lower one wins, which is the conservative choice: a completeness
    magnitude set too low is visible in the b-stability curve, while one set too
    high silently throws away events. The reference catalogue's mode is unique,
    so selecting the last maximal bin instead of the first left the suite green.
    """
    # Two events in the bin at 3.0 and two in the bin at 3.4, nothing between.
    mags = np.array([3.0, 3.0, 3.4, 3.4])
    assert mc_maxcurvature(mags, dm=0.1, correction=0.0) == 3.0
