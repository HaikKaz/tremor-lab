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


def test_a_column_name_that_is_not_in_the_file_is_named_in_the_error():
    # A bare KeyError('Latitude') tells a reader nothing about which of their
    # choices was wrong, or what the alternatives were.
    with pytest.raises(KeyError, match="not in the file"):
        read_catalog(
            DATA / "raw_catalog_koeri.csv",
            {**KOERI_COLUMNS, "lat": "Latitude"},
            MAINSHOCK,
        )


def test_a_missing_magnitude_column_is_reported_as_such():
    with pytest.raises(KeyError, match="no magnitude column"):
        read_catalog(
            DATA / "raw_catalog_koeri.csv",
            {"date": "Tarih", "time": "Saat"},
            MAINSHOCK,
        )


def test_a_missing_time_column_is_reported_as_such():
    with pytest.raises(KeyError, match="no time column"):
        read_catalog(DATA / "raw_catalog_koeri.csv", {"mag": "Mag"}, MAINSHOCK)


def test_a_catalogue_without_coordinates_is_still_usable(tmp_path):
    path = tmp_path / "nocoords.csv"
    path.write_text(
        "Tarih,Saat,Mag\n2023.02.06,13:17:32,4.5\n2023.02.07,01:17:32,4.2\n"
    )
    cat = read_catalog(path, {"date": "Tarih", "time": "Saat", "mag": "Mag"}, MAINSHOCK)
    assert len(cat) == 2
    assert cat["dist_km"].isna().all()
    assert cat["dt_days"].to_list() == pytest.approx([0.5, 1.0])


def test_a_distance_limit_is_refused_when_no_distance_can_be_measured(tmp_path):
    # Comparisons against NaN are false, so without this guard the limit would
    # remove every event and the report would blame the limit for it.
    path = tmp_path / "nocoords.csv"
    path.write_text("Tarih,Saat,Mag\n2023.02.06,13:17:32,4.5\n")
    with pytest.raises(ValueError, match="no distance could be computed"):
        read_catalog(
            path,
            {"date": "Tarih", "time": "Saat", "mag": "Mag"},
            MAINSHOCK,
            radius_km=500.0,
        )


def test_a_distance_limit_is_refused_when_the_mainshock_has_no_epicentre():
    with pytest.raises(ValueError, match="no distance could be computed"):
        read_catalog(
            DATA / "raw_catalog_koeri.csv",
            KOERI_COLUMNS,
            {**MAINSHOCK, "lat": float("nan"), "lon": float("nan")},
            radius_km=500.0,
        )


def test_dayfirst_leaves_a_year_first_date_exactly_as_written():
    # pandas applies dayfirst to every date, not only to the ambiguous ones: it
    # once read the fixture's 2023.02.06 as "%Y.%d.%m" and moved 6 February to
    # 2 June. A year-first date has no ambiguity to resolve, so the flag must not
    # touch it - and the reader must say out loud that it did not.
    with pytest.warns(UserWarning, match="state the year first"):
        cat = read_koeri(mag_type_col="Tip", dayfirst=True)
    assert cat["dt_days"].to_list() == pytest.approx([0.5, 1.0, 2.0, 30.0])
    assert cat.attrs["rows_parsed"] == 7


def test_dayfirst_leaves_an_iso_stamp_exactly_as_written():
    with pytest.warns(UserWarning, match="state the year first"):
        usgs = read_catalog(
            DATA / "raw_catalog_usgs.csv", USGS_COLUMNS, MAINSHOCK, dayfirst=True
        )
    assert usgs["dt_days"].to_list() == pytest.approx([0.5, 1.0, 2.0, 30.0])
    assert usgs["t"].iloc[0] == pd.Timestamp("2023-02-06 13:17:32")


def test_dayfirst_still_applies_to_the_rows_that_are_genuinely_ambiguous(tmp_path):
    # One file, both shapes: the year-first row is read as written and the
    # day-first row is read day-first, so neither correction costs the other.
    path = tmp_path / "mixed.csv"
    path.write_text(
        "date,time,lat,lon,mag\n"
        "2023.02.06,13:17:32,37.0,37.0,4.5\n"
        "07.02.2023,01:17:32,37.0,37.0,4.2\n"
    )
    columns = {"date": "date", "time": "time", "lat": "lat", "lon": "lon", "mag": "mag"}
    with pytest.warns(UserWarning, match="1 of 2 rows"):
        cat = read_catalog(path, columns, MAINSHOCK, dayfirst=True)
    assert cat["dt_days"].to_list() == pytest.approx([0.5, 1.0])


def test_the_row_count_and_the_window_agree_on_what_a_usable_row_is(tmp_path):
    # Coordinates written with a decimal comma parse as nothing at all. The
    # report used to call every such row unreadable and then analyse all of them:
    # "0 usable" and "2 events in window" in the same block of text.
    path = tmp_path / "decimal_comma.csv"
    path.write_text(
        "Tarih,Saat,Enlem,Boylam,Mag\n"
        '2023.02.06,13:17:32,"37,0","37,0",4.5\n'
        '2023.02.07,01:17:32,"37,5","37,5",4.2\n'
    )
    cat = read_catalog(path, KOERI_COLUMNS, MAINSHOCK)
    assert cat.attrs["rows_parsed"] == 2
    assert cat.attrs["events_in_window"] == 2
    assert len(cat) == 2


def test_a_catalogue_with_latitudes_but_no_longitudes_keeps_its_events(tmp_path):
    # Half a coordinate costs the distance and nothing else. Pooling the two
    # columns in one test meant a present latitude made the absent longitude
    # compulsory, so a file that was fine except for distances returned nothing.
    path = tmp_path / "halfcoords.csv"
    path.write_text(
        "Tarih,Saat,Enlem,Boylam,Mag\n"
        "2023.02.06,13:17:32,37.0,,4.5\n"
        "2023.02.07,01:17:32,38.0,,4.2\n"
    )
    cat = read_catalog(path, KOERI_COLUMNS, MAINSHOCK)
    assert cat["dt_days"].to_list() == pytest.approx([0.5, 1.0])
    assert cat.attrs["rows_parsed"] == 2
    assert cat["dist_km"].isna().all()


def test_an_unmeasurable_distance_is_reported_as_none_rather_than_nan(tmp_path):
    # The run report printed "farthest kept nan km" as though nan were a distance.
    path = tmp_path / "nocoords.csv"
    path.write_text(
        "Tarih,Saat,Mag\n2023.02.06,13:17:32,4.5\n2023.02.07,01:17:32,4.2\n"
    )
    cat = read_catalog(path, {"date": "Tarih", "time": "Saat", "mag": "Mag"}, MAINSHOCK)
    assert len(cat) == 2
    assert cat.attrs["farthest_km"] is None


def test_a_time_column_that_mixes_utc_offsets_is_refused_by_name(tmp_path):
    # Mixed offsets used to leave pandas with a column of plain objects rather
    # than timestamps, and the run died much later with "unsupported operand
    # type(s) for -: 'numpy.ndarray' and 'Timestamp'", naming neither the column
    # nor the clock. The offsets are now read from the text instead, because
    # pandas warns that it will raise on such a column in a later version, and a
    # check resting on the current behaviour would start failing on somebody
    # else's machine with a pandas error in place of this explanation.
    path = tmp_path / "two_clocks.csv"
    path.write_text(
        "time,latitude,longitude,mag\n"
        "2023-02-06T13:17:32+00:00,37.0,37.0,4.5\n"
        "2023-02-07T04:17:32+03:00,37.0,37.0,4.2\n"
    )
    with pytest.raises(ValueError, match="mixes UTC offsets"):
        read_catalog(path, USGS_COLUMNS, MAINSHOCK)


def test_one_blank_date_costs_its_own_row_and_no_others(tmp_path):
    # A YYYYMMDD column is whole numbers until a cell is left blank, when pandas
    # reads it as decimals; the reader then offered "20230206.0" to the date
    # parser and lost every row in the file rather than the one that was empty.
    path = tmp_path / "compact.csv"
    path.write_text(
        "date,time,lat,lon,mag\n"
        "20230206,13:17:32,37.0,37.0,4.5\n"
        ",01:17:32,37.0,37.0,4.2\n"
        "20230208,01:17:32,37.0,37.0,5.5\n"
    )
    columns = {"date": "date", "time": "time", "lat": "lat", "lon": "lon", "mag": "mag"}
    cat = read_catalog(path, columns, MAINSHOCK)
    assert cat.attrs["rows_read"] == 3
    assert cat.attrs["rows_parsed"] == 2
    assert cat["dt_days"].to_list() == pytest.approx([0.5, 2.0])


def test_one_blank_stamp_costs_its_own_row_in_a_single_column_file(tmp_path):
    path = tmp_path / "compact_single.csv"
    path.write_text(
        "time,latitude,longitude,mag\n"
        "20230206,37.0,37.0,4.5\n"
        ",37.0,37.0,4.2\n"
        "20230208,37.0,37.0,5.5\n"
    )
    cat = read_catalog(path, USGS_COLUMNS, MAINSHOCK)
    assert cat.attrs["rows_parsed"] == 2
    # Both dates are readable; midnight on the 6th is simply before the mainshock.
    assert cat["dt_days"].to_list() == pytest.approx([1.947], abs=0.01)


def test_the_scale_column_is_read_whether_it_is_a_role_or_an_argument():
    # The CLI passes the whole [catalog.columns] table, so the scale column can
    # arrive as a role. It used to be checked against the file and then never
    # read, and the magnitudes came back unconverted with nothing said about it.
    by_role = read_catalog(
        DATA / "raw_catalog_koeri.csv",
        {**KOERI_COLUMNS, "mag_type": "Tip"},
        MAINSHOCK,
    )
    assert by_role["mw"].to_list() == pytest.approx([4.5, 4.884, 5.705, 3.9])
    assert (
        by_role.attrs["scale_counts"]
        == read_koeri(mag_type_col="Tip").attrs["scale_counts"]
    )


def test_the_conversions_can_be_counted_over_the_analysed_events_alone():
    # scale_counts describes the download; a report that says "converted" about
    # the sequence needs the same tally over the events that were kept.
    cat = read_koeri(mag_type_col="Tip")
    assert cat.attrs["scale_counts"]["unconverted"] == 7
    assert cat.attrs["scale_counts_in_window"] == {
        "ms": 1,
        "mb": 1,
        "mw": 1,
        "unconverted": 1,
    }


def test_a_distance_limit_does_not_blame_the_coordinates_for_a_bad_magnitude(tmp_path):
    # Every row was already gone - the magnitudes read "M4.5" - so there was
    # nothing for the limit to remove. The guard used to fire anyway and send the
    # author off to check coordinates and an epicentre that were both correct.
    path = tmp_path / "text_magnitudes.csv"
    path.write_text(
        "Tarih,Saat,Enlem,Boylam,Mag\n"
        "2023.02.06,13:17:32,37.0,37.0,M4.5\n"
        "2023.02.07,01:17:32,38.0,37.0,M4.2\n"
    )
    cat = read_catalog(path, KOERI_COLUMNS, MAINSHOCK, radius_km=500.0)
    assert len(cat) == 0
    assert cat.attrs["rows_read"] == 2
    assert cat.attrs["rows_parsed"] == 0


# --------------------------------------------------- regressions found by review
# Both of these were introduced by a fix and caught by the pass that checked it.


def test_a_uniform_offset_is_dropped_even_when_a_row_is_unreadable():
    """The offset must come off however the column reached pandas.

    A column of tz-aware stamps arrives as a datetime column and the offset is
    stripped. A column of tz-aware stamps with one unreadable cell in it arrives
    as plain objects, and rebuilding those gives a tz-AWARE column back - which
    was returned as it came. The run then died sixty lines later on "Cannot
    subtract tz-naive and tz-aware datetime-like objects", which is the exact
    failure this code exists to prevent.
    """
    import warnings

    from tremor_lab.catalog import _one_clock, _read_dates

    raw = pd.Series(
        ["2023-02-06T13:17:32+03:00", "not-a-date", "2023-02-07T13:17:32+03:00"]
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parsed = _read_dates(raw, True, "'time'")
    result = _one_clock(parsed, "'time'")
    assert not isinstance(result.dtype, pd.DatetimeTZDtype)
    # The subtraction that used to raise.
    assert (result - pd.Timestamp("2023-02-06 01:17:32")).iloc[0].days == 0


def test_a_year_first_date_is_read_as_written_however_it_separates_its_parts():
    """Whitespace is a separator too.

    The pattern once required punctuation straight after the year, so
    "2023 02 06" was not recognised as year-first, dayfirst was forwarded to
    pandas, and 6 February silently became 2 June - with none of the warning
    that this recognition exists to raise.
    """
    import warnings

    from tremor_lab.catalog import _read_dates

    for text in ("2023 02 06 13:17:32", "2023-02-06 13:17:32", "2023.02.06 13:17:32"):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            parsed = _read_dates(pd.Series([text]), True, "'time'")
        assert str(parsed.iloc[0]).startswith("2023-02-06"), text
        assert caught, f"no warning that dayfirst did nothing for {text!r}"
