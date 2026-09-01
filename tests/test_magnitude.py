"""Analytic checks of the energy relations and the Scordilis conversions.

Every expected value here is exact arithmetic on a published relation, so these tests
fail only if a coefficient has been mistyped.
"""

import numpy as np
import pytest

from tremor_lab.magnitude import (
    bath_mag,
    bath_ratio,
    energy_joules,
    mb_to_mw,
    ms_to_mw,
    to_mw,
)


def test_energy_exponent():
    assert np.log10(energy_joules(7.8)) == pytest.approx(16.5)


def test_energy_scales_by_ten_to_the_three_halves_per_magnitude_unit():
    assert energy_joules(6.0) / energy_joules(5.0) == pytest.approx(10**1.5)


def test_energy_vectorises():
    assert np.log10(energy_joules([5.0, 7.8])) == pytest.approx([12.3, 16.5])


def test_bath_expected_magnitude():
    assert bath_mag(7.8) == pytest.approx(6.65)
    assert bath_mag(7.8, delta_mb=1.1) == pytest.approx(6.7)


def test_bath_ratio_at_a_one_unit_anomaly():
    assert np.log10(bath_ratio(7.55, 7.7, 1.15)) == pytest.approx(1.5)


def test_bath_ratio_is_unity_at_the_expected_magnitude():
    assert bath_ratio(bath_mag(7.8), 7.8) == pytest.approx(1.0)


def test_scordilis_ms_lower_branch():
    assert ms_to_mw(6.0) == pytest.approx(6.09)


def test_scordilis_ms_upper_branch():
    assert ms_to_mw(7.0) == pytest.approx(7.01)


def test_scordilis_ms_branches_switch_across_the_uncalibrated_gap():
    assert ms_to_mw(6.1) == pytest.approx(6.157)
    assert ms_to_mw(6.2) == pytest.approx(6.218)


def test_scordilis_mb():
    assert mb_to_mw(5.0) == pytest.approx(5.28)


def test_to_mw_dispatches_on_the_scale_label():
    assert to_mw(6.0, "Ms") == pytest.approx(6.09)
    assert to_mw(5.0, "mb") == pytest.approx(5.28)
    assert to_mw(6.5, "mww") == pytest.approx(6.5)


def test_to_mw_passes_through_scales_with_no_published_relation():
    assert to_mw(4.2, "ML") == pytest.approx(4.2)
    assert to_mw(4.2, "") == pytest.approx(4.2)


def test_to_mw_vectorises_over_a_mixed_scale_catalogue():
    mag = [6.0, 5.0, 6.5, 4.2]
    mtype = [" Ms ", "mb", "Mw", "ml"]
    assert to_mw(mag, mtype) == pytest.approx([6.09, 5.28, 6.5, 4.2])
