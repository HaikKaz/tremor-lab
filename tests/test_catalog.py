"""Catalogue reading, windowing, distance and magnitude homogenisation.

The two fixture files describe the same four in-window events in the two export
formats the tool has to cope with, so agreement between them is itself a test.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tremor_lab import constants
from tremor_lab.catalog import CATALOG_COLUMNS, haversine_km, read_catalog, window

DATA = Path(__file__).parent / "data"

MAINSHOCK = {"t": "2023-02-06 01:17:32", "lat": 37.0, "lon": 37.0, "mw": 7.8}

KOERI_COLUMNS = {
    "date": "Tarih",
    "time": "Saat",
    "lat": "Enlem",
    "lon": "Boylam",
    "mag": "Mag",
}
USGS_COLUMNS = {
    "datetime": "time",
    "lat": "latitude",
    "lon": "longitude",
    "mag": "mag",
}


def read_koeri(**kwargs):
    return read_catalog(
        DATA / "raw_catalog_koeri.csv", KOERI_COLUMNS, MAINSHOCK, **kwargs
    )


def test_one_degree_of_latitude_is_about_111_km():
    assert haversine_km(0, 0, 1, 0) == pytest.approx(111.19, abs=0.5)


def test_distance_from_a_point_to_itself_is_zero():
    assert haversine_km(37.0, 37.0, 37.0, 37.0) == pytest.approx(0.0)


def test_one_degree_of_longitude_shrinks_with_the_cosine_of_latitude():
    # Without the cos(lat) factor this would also return 111.19, so this is the
    # test that holds the longitude half of the formula in place.
    assert haversine_km(37.0, 37.0, 37.0, 38.0) == pytest.approx(88.80, abs=0.05)


def test_distance_is_symmetric():
    there = haversine_km(37.0, 36.0, 38.5, 38.2)
    back = haversine_km(38.5, 38.2, 37.0, 36.0)
    assert there == pytest.approx(back)


def test_distance_broadcasts_a_point_against_an_array():
    d = haversine_km(0.0, 0.0, np.array([0.0, 1.0]), np.array([0.0, 0.0]))
    assert d == pytest.approx([0.0, 111.19], abs=0.5)


def test_sphere_radius_can_be_overridden(monkeypatch):
    assert haversine_km(0, 0, 1, 0, radius_km=1.0) == pytest.approx(np.radians(1.0))
    monkeypatch.setattr(constants, "EARTH_RADIUS_KM", 1.0)
    assert haversine_km(0, 0, 1, 0) == pytest.approx(np.radians(1.0))


def test_reader_returns_the_documented_columns():
    assert tuple(read_koeri(mag_type_col="Tip").columns) == CATALOG_COLUMNS


def test_window_keeps_only_events_after_the_mainshock_and_inside_the_window():
    cat = read_koeri(mag_type_col="Tip", window_days=180)
    # 1 day before, and exactly at, the origin time are both excluded; so is day 181.
    assert cat["dt_days"].to_list() == pytest.approx([0.5, 1.0, 2.0, 30.0])


def test_window_length_is_configurable():
    assert len(read_koeri(mag_type_col="Tip", window_days=1.0)) == 2
    assert len(read_koeri(mag_type_col="Tip", window_days=200.0)) == 5


def test_default_window_length_follows_the_constant(monkeypatch):
    monkeypatch.setattr(constants, "WINDOW_DAYS", 1.0)
    assert len(read_koeri(mag_type_col="Tip")) == 2


def test_distances_are_measured_from_the_mainshock():
    cat = read_koeri(mag_type_col="Tip")
    assert cat["dist_km"].to_list() == pytest.approx([0.0, 111.19, 0.0, 0.0], abs=0.5)


def test_magnitudes_are_homogenised_when_a_type_column_is_given():
    cat = read_koeri(mag_type_col="Tip")
    # ML passes through; Ms 4.2 and mb 5.5 convert by Scordilis; Mw passes through.
    assert cat["mw"].to_list() == pytest.approx([4.5, 4.884, 5.705, 3.9])


def test_magnitudes_are_taken_as_reported_without_a_type_column():
    assert read_koeri()["mw"].to_list() == pytest.approx([4.5, 4.2, 5.5, 3.9])


def test_unparseable_rows_are_dropped_and_counted():
    cat = read_koeri(mag_type_col="Tip")
    assert cat.attrs["rows_read"] == 10
    assert cat.attrs["rows_parsed"] == 7
    assert cat.attrs["events_in_window"] == 4


def test_the_two_export_formats_give_the_same_window():
    koeri = read_koeri(mag_type_col="Tip")
    usgs = read_catalog(
        DATA / "raw_catalog_usgs.csv", USGS_COLUMNS, MAINSHOCK, mag_type_col="magType"
    )
    for column in ("mw", "dt_days", "dist_km"):
        assert usgs[column].to_list() == pytest.approx(koeri[column].to_list())


def test_a_utc_offset_is_dropped_rather_than_shifted():
    usgs = read_catalog(DATA / "raw_catalog_usgs.csv", USGS_COLUMNS, MAINSHOCK)
    assert usgs["t"].iloc[0] == pd.Timestamp("2023-02-06 13:17:32")


def test_dayfirst_resolves_an_ambiguous_european_date(tmp_path):
    path = tmp_path / "ambiguous.csv"
    path.write_text("date,time,lat,lon,mag\n06.02.2023,13:17:32,37.0,37.0,4.5\n")
    columns = {"date": "date", "time": "time", "lat": "lat", "lon": "lon", "mag": "mag"}
    # Read month-first, 06.02.2023 silently becomes 2 June: the event stays inside the
    # window but lands 116.5 days after the mainshock instead of half a day after it.
    month_first = read_catalog(path, columns, MAINSHOCK)
    assert month_first["dt_days"].to_list() == pytest.approx([116.5])
    day_first = read_catalog(path, columns, MAINSHOCK, dayfirst=True)
    assert day_first["dt_days"].to_list() == pytest.approx([0.5])


def test_window_can_be_called_directly_on_a_parsed_frame():
    frame = pd.DataFrame(
        {
            "t": pd.to_datetime(["2023-02-06 13:17:32", "2023-02-07 01:17:32"]),
            "lat": [37.0, 38.0],
            "lon": [37.0, 37.0],
            "mw": [4.5, 4.2],
        }
    )
    out = window(frame, MAINSHOCK, window_days=180)
    assert out["dt_days"].to_list() == pytest.approx([0.5, 1.0])


def test_conversions_are_counted_so_a_report_can_state_what_happened():
    cat = read_koeri(mag_type_col="Tip")
    # the fixture's four in-window events are ML, Ms, mb and Mw
    assert cat.attrs["scale_counts"]["ms"] == 1
    assert cat.attrs["scale_counts"]["mb"] == 1
    assert cat.attrs["scale_counts"]["mw"] == 1


def test_a_column_holding_no_scale_labels_converts_nothing():
    # KOERI's "Type" column holds "Earthquake"; naming it as the scale column
    # must not produce a report claiming a conversion took place.
    path = DATA / "raw_catalog_koeri.csv"
    frame = pd.read_csv(path)
    frame["Kind"] = "Earthquake"
    scratch = DATA.parent / "_kind.csv"
    frame.to_csv(scratch, index=False)
    try:
        cat = read_catalog(scratch, KOERI_COLUMNS, MAINSHOCK, mag_type_col="Kind")
        counts = cat.attrs["scale_counts"]
        assert counts["ms"] == 0 and counts["mb"] == 0
        assert counts["unconverted"] == len(frame)
        assert cat["mw"].to_list() == pytest.approx([4.5, 4.2, 5.5, 3.9])
    finally:
        scratch.unlink()


def test_no_distance_limit_is_applied_unless_asked_for():
    cat = read_koeri(mag_type_col="Tip")
    assert cat.attrs["radius_km"] is None
    assert cat.attrs["removed_by_radius"] == 0
    assert len(cat) == 4


def test_a_distance_limit_removes_only_events_beyond_it_and_says_how_many():
    # the fixture's four in-window events sit at 0, 111.2, 0 and 0 km
    near = read_koeri(mag_type_col="Tip", radius_km=50.0)
    assert len(near) == 3
    assert near.attrs["radius_km"] == 50.0
    assert near.attrs["removed_by_radius"] == 1
    assert near["dist_km"].max() <= 50.0
    assert near.attrs["farthest_km"] == pytest.approx(0.0)


def test_a_generous_distance_limit_removes_nothing():
    cat = read_koeri(mag_type_col="Tip", radius_km=1000.0)
    assert len(cat) == 4
    assert cat.attrs["removed_by_radius"] == 0
    assert cat.attrs["farthest_km"] == pytest.approx(111.19, abs=0.5)
