"""Magnitudes are decimal; binary floating point is not. These pin the difference.

Reported magnitudes lie on a decimal grid, and neither 4.1 nor a bin edge computed
from it exists exactly in binary. Every comparison of a magnitude against an edge
was therefore being decided by rounding error in the sixteenth decimal place.

At the default bin width of 0.1 the errors happened to fall the harmless way on
the reference catalogue, so nothing here went red. At 0.2 they did not, and the
consequences were not subtle: an entire magnitude class dropped out of its own
completeness sample, b wrong by six times its own standard error, and a
maximum-curvature mode that was an artefact of arithmetic rather than a property
of the catalogue.

Everything here is at bin widths other than the default, because that is where the
defects lived and where nothing was looking.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tremor_lab import constants
from tremor_lab.bvalue import b_stability, b_value_aki, b_value_tinti, mc_b_stability
from tremor_lab.completeness import _edge_decimals, fmd, mc_goodness_of_fit
from tremor_lab.grid import at_or_above, bin_index

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def mags():
    return pd.read_csv(DATA / "kahramanmaras_180d.csv")["mw"].to_numpy()


# --------------------------------------------------------------- the helpers


def test_a_magnitude_sitting_exactly_on_a_threshold_is_at_or_above_it():
    # 4.2 - 0.2/2 is 4.1000000000000005 in binary, so the plain comparison says
    # an event reported at M 4.1 is below its own completeness threshold.
    threshold = 4.2 - 0.2 / 2
    assert not 4.1 >= threshold
    assert at_or_above([4.1], threshold)[0]


def test_the_tolerance_cannot_swallow_a_magnitude_a_catalogue_would_distinguish():
    # The finest resolution any agency reports is 0.01, and 0.001 is already
    # beyond it. The tolerance is six orders of magnitude smaller again.
    assert not at_or_above([4.099], 4.1)[0]
    assert not at_or_above([4.0999], 4.1)[0]


def test_a_magnitude_half_way_between_two_bin_centres_goes_up():
    # Centres run 3.0, 3.2, 3.4, 3.6. M 3.5 sits exactly between 3.4 and 3.6 and
    # M 3.3 exactly between 3.2 and 3.4; both go up. The spreadsheet rounds up, so
    # the package and the browser page do too. Before the fix the division decided
    # it: 3.5 came out exactly 2.5 and went up, while 3.3 came out
    # 1.4999999999999991 and went down.
    assert bin_index([3.5], lo=3.0, dm=0.2)[0] == 3
    assert bin_index([3.3], lo=3.0, dm=0.2)[0] == 2


def test_the_tolerance_is_a_published_constant_a_reader_can_change(mags):
    assert constants.GRID_TOLERANCE in constants._PUBLISHED.values()
    # Widened far enough, it changes an answer - which is what makes it a real
    # setting rather than decoration.
    assert at_or_above([4.0], 4.1, tolerance=0.2)[0]
    assert not at_or_above([4.0], 4.1)[0]


# ------------------------------------------------------------------ the b-value


@pytest.mark.parametrize(
    ("mc", "expected_n"),
    [(4.2, 498), (4.4, 309), (5.2, 40), (5.4, 28)],
)
def test_no_magnitude_class_is_dropped_from_its_own_completeness_sample(
    mags, mc, expected_n
):
    # Each of these thresholds lands mc - dm/2 on a reported magnitude. Before the
    # grid comparison, every event at that magnitude was excluded: at mc 4.2 the
    # sample was 394 events instead of 498, and b came out 0.936 instead of 1.182.
    assert b_value_aki(mags, mc, dm=0.2).n == expected_n


def test_the_two_b_estimators_select_the_same_events(mags):
    for mc in (3.4, 4.0, 4.2, 4.4, 5.2):
        assert b_value_aki(mags, mc, dm=0.2).n == b_value_tinti(mags, mc, dm=0.2).n


def test_the_b_stability_curve_has_no_step_that_only_arithmetic_explains(mags):
    # The dropped magnitude classes put a dip in the diagnostic whose whole
    # purpose is to reveal a mis-set Mc: 1.1431, 0.9355, 0.9951, 1.3087, where
    # the surrounding trend is monotone.
    curve = b_stability(mags, dm=0.2)
    pairs = zip(curve.thresholds, curve.b, strict=True)
    window = [b for t, b in pairs if 3.95 < t < 4.65 and np.isfinite(b)]
    assert window == sorted(window), "b jumps down and back up between adjacent bins"


def test_stability_cannot_be_claimed_when_there_is_nothing_to_average_over(mags):
    # When the window rounds to zero bins it collapses onto the single value it
    # is being compared against, the difference is exactly zero, and the lowest
    # candidate comes back as though it had passed a test that never ran. It now
    # refuses instead. At dm 1.0 the window still covers one bin, which is a real
    # if crude comparison, and the honest answer there is that nothing stabilised.
    for dm in (1.2, 1.5, 2.0):
        with pytest.raises(ValueError, match="nothing to average over"):
            mc_b_stability(mags, dm=dm)
    assert mc_b_stability(mags, dm=1.0) is None


# ------------------------------------------------------------------- the binning


def test_equal_counts_at_every_magnitude_give_equal_bins(mags):
    # One event at each 0.1 step. Every bin of width 0.2 holds two of them, except
    # the lowest, which is half a bin. Before the fix: 1, 3, 1, 2, 3, 2, 2, 1, 2.
    grid = np.round(np.arange(3.0, 4.7, 0.1), 1)
    assert list(fmd(grid, dm=0.2).inc) == [1] + [2] * 8


@pytest.mark.parametrize("dm", [0.05, 0.1, 0.2, 0.25, 0.5])
def test_the_cumulative_curve_agrees_with_counting_the_events_directly(mags, dm):
    # The goodness-of-fit test estimates b from the events at or above a bin's
    # lower edge and scores it against the cumulative count in that bin. If the
    # two disagree it is scoring one sample against another sample's curve.
    edges, _, cum = fmd(mags, dm=dm)
    for edge, count in zip(edges, cum, strict=True):
        assert int(at_or_above(mags, edge - dm / 2).sum()) == count


def test_the_goodness_of_fit_statistic_is_scored_against_its_own_sample(mags):
    # At dm 0.2 this read 84.5 where a self-consistent calculation gives 91.7, an
    # error of seven points in a statistic whose accept threshold is 90.
    result = mc_goodness_of_fit(mags, dm=0.2)
    at_34 = list(np.round(result.candidates, 1)).index(3.4)
    assert result.r_values[at_34] == pytest.approx(91.7, abs=0.5)


@pytest.mark.parametrize(
    ("dm", "decimals"), [(0.1, 2), (0.2, 2), (0.05, 3), (0.25, 3), (0.125, 4)]
)
def test_a_bin_label_carries_enough_digits_to_be_the_edge_it_names(dm, decimals):
    # From how the width is written, not from its logarithm, which allowed two
    # decimals for 0.125 and so labelled the bin standing for 3.125 as 3.13 -
    # and that label was returned as the completeness magnitude.
    assert _edge_decimals(dm) == decimals


def test_bin_labels_are_the_values_they_stand_for_at_an_awkward_width():
    grid = np.round(np.arange(3.0, 4.0, 0.125), 3)
    edges = fmd(grid, dm=0.125).edges
    assert np.allclose(edges, np.round(edges, 3))
    assert 3.125 in [float(e) for e in edges]


# ------------------------------------------------------- and the default is safe


def test_the_locked_reference_values_are_unchanged_by_any_of_this(mags):
    from tremor_lab import mc_maxcurvature

    assert mc_maxcurvature(mags) == 3.4
    b = b_value_aki(mags[mags >= 3.5], 3.5)
    assert b.n == 1529
    assert b.b == pytest.approx(0.844132, abs=5e-6)
    assert b.sigma == pytest.approx(0.018540, abs=5e-6)


@pytest.mark.parametrize(
    ("dm", "expected"),
    [(0.05, 3.2), (0.1, 4.1), (0.2, 4.2), (0.25, 4.0), (0.5, 4.5)],
)
def test_the_stability_completeness_matches_the_browser_page_at_every_bin_width(
    mags, dm, expected
):
    # These are the values the JavaScript implementation returns, checked in a
    # browser against the same catalogue. They agreed at the default width and
    # nowhere else: at dm 0.2 the package said 4.0 and the page said 4.2 from
    # b-stability curves that were identical row for row, because span/dm is
    # exactly 2.5 there and Python's round() sends a half down while
    # JavaScript's Math.round sends it up.
    assert mc_b_stability(mags, dm=dm) == pytest.approx(expected)


def test_the_averaging_window_rounds_the_same_way_magnitudes_are_binned(mags):
    # Halves up, so a window of 0.5 at dm 0.2 covers three bins rather than two.
    # Asserted through the public result rather than the private arithmetic: with
    # two bins the answer is 4.0, with three it is 4.2.
    assert mc_b_stability(mags, dm=0.2, span=0.5) == pytest.approx(4.2)
    assert mc_b_stability(mags, dm=0.2, span=0.4) == pytest.approx(4.0)
