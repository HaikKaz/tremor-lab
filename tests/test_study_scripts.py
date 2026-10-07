"""The study scripts in examples/ make claims, so the parts they rest on are tested.

`calibration_study.py` and `synthetic_incompleteness.py` take an hour between them and
produce the numbers behind two sections of the paper. A mistake in how they simulate or
summarise would not raise an error; it would produce a confident table. These tests run
the pieces on small inputs where the answer is known: the interval arithmetic, the
simulators, and the fit that recovers a detection function from the curve it implies.
"""

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, EXAMPLES / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


calibration = _load("calibration_study")
synthetic = _load("synthetic_incompleteness")


@pytest.mark.parametrize("module", [calibration, synthetic])
def test_wilson_interval_reproduces_known_values(module):
    # 50 of 100 and 0 of 20 are the textbook cases, 0.4038 to 0.5962 and 0 to 0.1611.
    lo, hi = module.wilson(50, 100)
    assert lo == pytest.approx(0.4038, abs=5e-4)
    assert hi == pytest.approx(0.5962, abs=5e-4)
    lo, hi = module.wilson(0, 20)
    assert lo == 0.0
    assert hi == pytest.approx(0.1611, abs=5e-4)
    assert all(math.isnan(v) for v in module.wilson(0, 0))


def test_a_null_sequence_has_the_requested_size_and_stays_in_the_window():
    config = calibration.CONFIGS[0]
    times = calibration.draw_sequence(config, np.random.default_rng(1))
    assert times.size == config["n"]
    assert np.all(np.diff(times) >= 0)
    assert times[0] > 0
    assert times[-1] <= calibration.WINDOW_DAYS


def test_the_second_sequence_scenario_puts_its_share_of_events_after_the_start():
    config = next(c for c in calibration.CONFIGS if c["second"] > 0)
    times = calibration.draw_sequence(config, np.random.default_rng(2))
    assert times.size == config["n"]
    # Events of the first sequence can fall after the second one begins as well, so
    # this is a floor: the second sequence's own share cannot be earlier than its start.
    late = np.sum(times >= config["second_start"])
    assert late >= round(config["n"] * config["second"])


def test_each_calibration_sequence_has_its_own_seed_and_does_not_depend_on_order():
    job = (0, 7, 9)
    first = calibration.one_sequence(job)
    again = calibration.one_sequence(job)
    other = calibration.one_sequence((0, 8, 9))
    assert first == again
    assert first["p_calibrated"] != other["p_calibrated"] or first["ks"] != other["ks"]


def test_the_asymptotic_p_value_is_more_generous_than_the_calibrated_one_on_a_null():
    # The reason the study exists: on a correct decay the textbook p-value sits near
    # 0.85 and the calibrated one near 0.5.
    rows = [calibration.one_sequence((0, i, 49)) for i in range(12)]
    asymptotic = np.mean([r["p_asymptotic"] for r in rows])
    calibrated = np.mean([r["p_calibrated"] for r in rows])
    assert asymptotic > calibrated + 0.15


def test_detection_function_thins_magnitudes_the_way_it_says():
    scenario = {
        "b": 1.0,
        "mu": 1.7,
        "sigma": 0.4,
        "n_generated": 2_000_000,
        "floor": 0.5,
        "record_from": 0.5,
        "decimals": 3,
    }
    mags = synthetic.catalogue(scenario, np.random.default_rng(3))
    # Half of the events at mu are recorded: the recorded count per 0.05 around mu is
    # half the generated count there, and by mu + 2.5 sigma nearly all are.
    generated = 2_000_000 * (10 ** (-(1.675 - 0.5)) - 10 ** (-(1.725 - 0.5)))
    recorded = np.sum((mags >= 1.675) & (mags < 1.725))
    assert recorded == pytest.approx(0.5 * generated, rel=0.06)
    full = 2_000_000 * (10 ** (-(2.7 - 0.5)) - 10 ** (-(2.75 - 0.5)))
    kept = np.sum((mags >= 2.7) & (mags < 2.75))
    assert kept == pytest.approx(0.995 * full, rel=0.1)


def test_a_complete_scenario_keeps_every_event_and_rounds_to_the_grid():
    scenario = synthetic.SCENARIOS[0]
    mags = synthetic.catalogue(scenario, np.random.default_rng(4))
    assert mags.size == scenario["n_generated"]
    assert mags.min() >= scenario["record_from"] - 1e-9
    assert np.allclose(mags * 10, np.round(mags * 10))


def test_the_expected_b_matches_a_simulation_under_the_same_detection_function():
    b, mu, sigma = 1.0, 1.7, 0.4
    scenario = {
        "b": b,
        "mu": mu,
        "sigma": sigma,
        "n_generated": 3_000_000,
        "floor": 0.5,
        "record_from": 0.5,
        "decimals": 6,
    }
    mags = synthetic.catalogue(scenario, np.random.default_rng(5))
    thresholds = np.array([1.5, 2.0, 2.4])
    expected = synthetic.expected_b(thresholds, b, mu, sigma)
    for threshold, want in zip(thresholds, expected, strict=True):
        sample = mags[mags >= threshold - synthetic.DM / 2]
        got = 1.0 / (math.log(10.0) * (sample.mean() - (threshold - synthetic.DM / 2)))
        assert got == pytest.approx(want, abs=0.015), threshold
    # The curve rises and levels off at the true b.
    assert expected[0] < expected[1] < expected[2] < b + 0.01


def test_the_detection_fit_recovers_the_function_that_made_the_curve():
    thresholds = np.round(np.arange(1.5, 3.1, 0.1), 2)
    curve = synthetic.expected_b(thresholds, 1.0, 1.7, 0.4)
    fit = synthetic.fit_detection(thresholds, curve, np.full(thresholds.size, 0.01))
    assert fit["b_true"] == pytest.approx(1.0, abs=0.02)
    assert fit["mu"] == pytest.approx(1.7, abs=0.05)
    assert fit["sigma"] == pytest.approx(0.4, abs=0.05)
    assert fit["weighted_ss"] < 1e-3
    assert not fit["at_bound"]


def test_the_oracle_threshold_is_the_first_grid_value_at_ninety_nine_per_cent():
    for scenario in synthetic.SCENARIOS:
        threshold = synthetic.oracle_threshold(scenario)
        if scenario["mu"] is None:
            assert threshold == scenario["record_from"]
            continue
        level = scenario["mu"] + 2.326 * scenario["sigma"]
        assert threshold >= level - 1e-9
        assert threshold - synthetic.DM < level + 1e-9


def test_the_chi_squared_size_check_counts_only_catalogues_with_a_candidate():
    result = synthetic.chi_squared_size(replicates=40)
    assert result["scenario"] == "complete"
    assert 0 < result["replicates"] <= 40
    lo, hi = result["rejected_ci95"]
    assert 0.0 <= lo <= result["rejected_at_0.05"] <= hi <= 1.0


def test_a_replicate_of_the_complete_scenario_finds_the_true_b_and_reports_everything():
    row = synthetic.one_replicate((0, 3, 49))
    assert row["n"] > 8000
    for label in (
        "maxcurvature",
        "goodness_of_fit",
        "b_stability",
        "plateau",
        "oracle",
    ):
        found = row[label]
        assert found is not None, label
        assert found["b"] == pytest.approx(1.0, abs=0.12), label
        assert found["sigma"] > 0
    assert len(row["curve"]) == len(synthetic.CURVE)
    assert row["first_candidate"]["threshold"] == pytest.approx(1.7)
