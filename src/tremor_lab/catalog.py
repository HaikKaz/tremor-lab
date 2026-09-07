"""Reading earthquake catalogues and reducing them to an aftershock window.

The analysis-ready table has the columns named in `CATALOG_COLUMNS`; `analyze_case`
consumes exactly that shape, so the reader and the analysis cannot drift apart.
"""

import re
import warnings
from collections.abc import Mapping

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray

from tremor_lab import constants
from tremor_lab.magnitude import to_mw

CATALOG_COLUMNS = ("t", "lat", "lon", "mw", "dt_days", "dist_km")

# A date whose first field is a four-digit year: 2023.02.06, 2023-02-06T13:17:32Z,
# 2023/02/06. Written that way the day and the month cannot be confused, which is
# why the day-first flag must never be applied to one. See `_read_dates`.
# A stamp that opens with a four-digit year, however it separates the parts.
# The separator was once required to be punctuation, which let "2023 02 06"
# through to be read day-first and silently turned into 6 February's data for
# 2 June, with none of the warning this pattern exists to trigger.
_YEAR_FIRST = re.compile(r"\s*\d{4}\s*(?:[-/.]|\s)\s*\d{1,2}")


def haversine_km(
    lat1: ArrayLike,
    lon1: ArrayLike,
    lat2: ArrayLike,
    lon2: ArrayLike,
    radius_km: float | None = None,
) -> NDArray[np.float64]:
    """
    Great-circle distance between two points, or between a point and an array.

    Parameters
    ----------
    lat1, lon1, lat2, lon2 : array_like
        Latitudes and longitudes in decimal degrees.
    radius_km : float, optional
        Sphere radius. Defaults to `constants.EARTH_RADIUS_KM` (6371 km).

    Returns
    -------
    ndarray
        Distance in kilometres.
    """
    radius_km = constants.EARTH_RADIUS_KM if radius_km is None else radius_km
    lat1, lon1, lat2, lon2 = (np.asarray(v, float) for v in (lat1, lon1, lat2, lon2))
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = (
        np.sin(dlat / 2) ** 2
        + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon / 2) ** 2
    )
    return 2 * radius_km * np.arcsin(np.sqrt(a))


def _usable_rows(catalog: pd.DataFrame) -> pd.Series:
    """Which rows can be analysed at all: the one rule, asked once.

    This is the only place that decides what "usable" means, because there were
    once two places. The window dropped a row when a coordinate was missing and
    the file carried coordinates *somewhere*; the row counter dropped it whenever
    a coordinate *role* had been named in the settings, even if that column held
    nothing readable. On a file whose coordinates were all unreadable the printed
    report then said "0 usable" on one line and "120 events in window" three lines
    below it, and a referee could not tell which number was the real one.

    Each coordinate is judged on its own column, too. Pooling them meant that a
    catalogue with latitudes but no longitudes had every single row deleted -
    worse than a catalogue with no coordinates at all, which is handled fine.
    A catalogue with no longitude column at all only loses the distance. A
    catalogue that has one, but leaves it blank on some rows, loses those
    rows: a coordinate that is present for most events and absent for a few
    is a defect in the file, not a catalogue without coordinates.
    """
    required = ["t", "mw"]
    for role in ("lat", "lon"):
        if role in catalog and catalog[role].notna().any():
            required.append(role)
    return catalog[required].notna().all(axis=1)


def window(
    catalog: pd.DataFrame,
    mainshock: Mapping,
    window_days: float | None = None,
    radius_km: float | None = None,
) -> pd.DataFrame:
    """
    Restrict a catalogue to the post-mainshock window and attach the offsets.

    Events at or before the mainshock origin time are excluded, so `dt_days` is
    strictly positive and the Omori fit never sees a zero elapsed time.

    Parameters
    ----------
    catalog : DataFrame
        Must carry t (datetime) and mw. Coordinates are optional: without them
        the distance to the mainshock is simply not computed, and a distance
        limit cannot be applied.
    mainshock : mapping
        Keys t, lat, lon, mw. The origin time must be on the same clock as the
        catalogue; no timezone conversion is applied.
    window_days : float, optional
        Length of the aftershock window in days. Defaults to
        `constants.WINDOW_DAYS` (180).
    radius_km : float, optional
        Keep only events within this great-circle distance of the mainshock
        epicentre. No limit by default: a spatial cut is a choice that changes
        the result, so it is never applied unless asked for, and the number of
        events it removed is recorded in ``attrs``.

    Returns
    -------
    DataFrame
        Columns as in `CATALOG_COLUMNS`, index reset. ``attrs`` also carries
        ``kept_rows``, the index labels the kept events had in `catalog`, so a
        caller can describe the events that were analysed rather than the file
        they came out of.
    """
    window_days = constants.WINDOW_DAYS if window_days is None else window_days
    # A row is kept or dropped by `_usable_rows` and by nothing else, so the count
    # a report prints and the rows an analysis sees can never tell two stories.
    out = catalog.loc[_usable_rows(catalog)].copy()
    for role in ("lat", "lon"):
        if role not in out:
            out[role] = np.nan
    origin = pd.Timestamp(mainshock["t"])
    out["dt_days"] = (out["t"] - origin).dt.total_seconds() / 86400.0
    out["dist_km"] = haversine_km(
        mainshock["lat"], mainshock["lon"], out["lat"].to_numpy(), out["lon"].to_numpy()
    )
    inside = (out["dt_days"] > 0) & (out["dt_days"] <= window_days)
    in_time = int(inside.sum())
    if radius_km is not None:
        # A comparison against NaN is false, so an unmeasurable distance would
        # silently remove every event and blame the limit for it.
        #
        # The `len(out)` half matters as much: on an empty frame "no distance is
        # finite" is trivially true, and the message below then sent the author to
        # check coordinates and an epicentre that were both perfectly good, when
        # the real fault was that every row had an unreadable time or magnitude.
        # With nothing left to filter, the row counts say what actually happened.
        if len(out) and not np.isfinite(out["dist_km"]).any():
            raise ValueError(
                "a distance limit of "
                f"{radius_km} km was asked for, but no distance could be "
                "computed: the catalogue has no usable coordinates, or the "
                "mainshock latitude and longitude are missing"
            )
        inside &= out["dist_km"] <= radius_km
    kept = out.loc[inside, list(CATALOG_COLUMNS)]
    farthest = float(kept["dist_km"].max()) if len(kept) else np.nan
    attrs = {
        "radius_km": radius_km,
        "removed_by_radius": in_time - len(kept),
        # None rather than NaN when nothing could be measured. Events with no
        # usable coordinates have no farthest one, and the run report printed the
        # bare word "nan" next to a kilometre sign as though it were a measurement.
        "farthest_km": farthest if np.isfinite(farthest) else None,
        "kept_rows": kept.index,
    }
    kept = kept.reset_index(drop=True)
    kept.attrs = attrs
    return kept


def read_catalog(
    path,
    columns: Mapping[str, str],
    mainshock: Mapping,
    window_days: float | None = None,
    mag_type_col: str | None = None,
    dayfirst: bool = False,
    radius_km: float | None = None,
) -> pd.DataFrame:
    """
    Read a catalogue file and reduce it to the analysis-ready aftershock window.

    Two catalogue shapes are supported: a single ISO datetime column (USGS exports,
    reported in UTC), given as ``columns["datetime"]``; or separate date and time
    columns (KOERI exports, reported in local time), given as ``columns["date"]`` and
    ``columns["time"]``. Both also need ``columns["lat"]``, ``columns["lon"]`` and
    ``columns["mag"]``.

    Rows whose time, position or magnitude cannot be parsed are dropped. Counts of rows
    read, rows parsed and events kept are attached to the result's ``attrs``, so a
    caller can report how much of the file survived. ``rows_parsed`` is decided by the
    same rule as the window itself, so the two counts cannot contradict each other.

    Parameters
    ----------
    path : str or path-like
        CSV file.
    columns : mapping of str to str
        Column names for the roles described above. ``columns["mag_type"]`` is
        accepted as an alternative to the `mag_type_col` argument below.
    mainshock : mapping
        Keys t, lat, lon, mw, with t on the same clock as the catalogue.
    window_days : float, optional
        Aftershock window length in days. Defaults to `constants.WINDOW_DAYS`.
    mag_type_col : str, optional
        Column of magnitude-scale labels. When given, magnitudes are homogenised to Mw
        by `tremor_lab.magnitude.to_mw`; when omitted they are taken as already on the
        moment scale.
    dayfirst : bool, optional
        Read day-first dates such as 06.02.2023 as 6 February rather than 2 June.
        False by default. It applies only to the rows where the day and the month
        really could be swapped: a year-first stamp such as 2023.02.06, or an ISO
        one such as 2023-02-06T13:17:32Z, is read as written whatever this says,
        and a warning states how many rows the flag did not touch.
    radius_km : float, optional
        Great-circle distance limit from the mainshock epicentre. No limit by
        default; see `window`.

    Returns
    -------
    DataFrame
        Columns as in `CATALOG_COLUMNS`.
    """
    frame = pd.read_csv(path)
    # The command line hands its whole [catalog.columns] table over, so the scale
    # column arrives as columns["mag_type"]; a library caller who wrote it there
    # rather than in `mag_type_col` used to have the role checked against the file
    # and then never read, and got un-homogenised magnitudes without being told.
    mag_type_col = mag_type_col or columns.get("mag_type")
    _check_columns(frame, columns, mag_type_col)
    parsed = pd.DataFrame({"t": _parse_times(frame, columns, dayfirst)})
    for role in ("lat", "lon"):
        parsed[role] = (
            pd.to_numeric(frame[columns[role]], errors="coerce")
            if role in columns
            else np.nan
        )
    mag = pd.to_numeric(frame[columns["mag"]], errors="coerce")
    labels = None
    scales: dict[str, int] = {}
    in_window_scales: dict[str, int] = {}
    if mag_type_col:
        labels = frame[mag_type_col].to_numpy(dtype=str)
        parsed["mw"] = to_mw(mag.to_numpy(), labels)
        scales = _scale_counts(labels)
    else:
        parsed["mw"] = mag

    out = window(parsed, mainshock, window_days, radius_km)
    if labels is not None:
        # The same tally restricted to the events that were actually analysed.
        # `scale_counts` describes the download, which is a different sentence.
        in_window_scales = _scale_counts(labels[np.asarray(out.attrs["kept_rows"])])
    out.attrs = {
        **out.attrs,
        "rows_read": len(frame),
        "rows_parsed": int(_usable_rows(parsed).sum()),
        "events_in_window": len(out),
        "scale_counts": scales,
        "scale_counts_in_window": in_window_scales,
    }
    return out


def _scale_counts(labels) -> dict[str, int]:
    """How many magnitudes each conversion branch of `to_mw` actually took.

    A scale column is not always a scale column: KOERI's "Type" holds
    "Earthquake". Counting the branches lets a caller report what happened
    instead of asserting that a conversion occurred.

    Called on the whole file it describes the download, including rows later
    dropped as unreadable or left outside the window; that is what lands in
    ``attrs["scale_counts"]``. ``attrs["scale_counts_in_window"]`` is the same
    tally over the analysed events only, for a report that wants to say
    "converted" about the sequence rather than about the file.
    """
    counts = {"ms": 0, "mb": 0, "mw": 0, "unconverted": 0}
    for raw in labels:
        label = str(raw).strip().lower()
        if label.startswith("mw"):
            counts["mw"] += 1
        elif label == "ms":
            counts["ms"] += 1
        elif label == "mb":
            counts["mb"] += 1
        else:
            counts["unconverted"] += 1
    return counts


def _check_columns(
    frame: pd.DataFrame, columns: Mapping[str, str], mag_type_col: str | None
) -> None:
    """Reject a column mapping before it produces a bare KeyError.

    Naming a column that is not in the file, or leaving out a role the reader
    needs, otherwise surfaces as `KeyError: 'lat'`, which tells a reader nothing
    about which of their choices was wrong.
    """
    if "mag" not in columns:
        raise KeyError(
            "no magnitude column was given; set columns['mag'] to one of "
            f"{list(frame.columns)}"
        )
    if not ({"datetime", "date"} & set(columns)):
        raise KeyError(
            "no time column was given; set columns['datetime'] for a single "
            "ISO stamp, or columns['date'] and columns['time'] for separate "
            f"date and clock columns, from {list(frame.columns)}"
        )
    known = {"datetime", "date", "time", "lat", "lon", "mag", "mag_type", "dt_days"}
    unknown = sorted(set(columns) - known)
    if unknown:
        raise KeyError(
            f"unknown column roles {unknown}; a mistyped role is dropped in "
            f"silence and its column never read. The roles are {sorted(known)}"
        )
    named = {role: name for role, name in columns.items() if role != "dt_days"}
    if mag_type_col:
        named["mag_type"] = mag_type_col
    missing = {r: n for r, n in named.items() if n and n not in frame.columns}
    if missing:
        raise KeyError(
            f"these column names are not in the file: {missing}. "
            f"The file has {list(frame.columns)}"
        )


def _as_text(column: pd.Series) -> pd.Series:
    """The date or time column as the text a person wrote in it.

    A compact YYYYMMDD column is read by pandas as whole numbers, and as decimals
    the moment one cell is left blank. Handing that straight to `astype(str)`
    turned 20230206 into "20230206.0", which no date parser accepts, so one empty
    cell made every row in the file unreadable instead of only its own.

    Missing cells come back as an empty string rather than the word "nan", so the
    caller can tell an absent date from a misspelt one.
    """
    if pd.api.types.is_numeric_dtype(column):
        text = column.map(_number_as_written)
    else:
        text = column.astype(str).str.strip()
    return text.mask(column.isna(), "")


def _number_as_written(value) -> str:
    """Undo the decimal point pandas adds to a numeric date column."""
    if pd.isna(value) or not np.isfinite(value):
        return ""
    number = float(value)
    return f"{number:.0f}" if number == int(number) else str(value)


def _parse_times(
    frame: pd.DataFrame, columns: Mapping[str, str], dayfirst: bool
) -> pd.Series:
    """Combine the configured date and time columns into one naive datetime series."""
    if "datetime" in columns:
        source = repr(columns["datetime"])
        text = _as_text(frame[columns["datetime"]])
        # A blank cell must make its own row unreadable. Left as an empty string it
        # is not empty for long: pandas reads a bare "13:17:32" as today's date.
        raw = text.where(text != "")
    else:
        source = f"{columns['date']!r} and {columns['time']!r}"
        date_text = _as_text(frame[columns["date"]])
        time_text = _as_text(frame[columns["time"]])
        complete = (date_text != "") & (time_text != "")
        raw = (date_text + " " + time_text).where(complete)
    _refuse_mixed_offsets(raw, source)
    return _one_clock(_read_dates(raw, dayfirst, source), source)


_UTC_OFFSET = re.compile(r"(Z|[+-]\d{2}:?\d{2})\s*$")


def _refuse_mixed_offsets(raw: pd.Series, source: str) -> None:
    """Refuse a time column carrying more than one UTC offset, read from the text.

    A UTC offset cannot be subtracted from a naive mainshock time, and when the
    offsets differ there is no single clock to drop. This is decided from the
    strings rather than from what pandas makes of them, for two reasons. Pandas
    currently answers a mixed-offset column with a column of plain objects and a
    FutureWarning saying it will raise instead in a later version, so a check
    built on that behaviour would start failing on somebody else's machine with a
    pandas error in place of this explanation. And the offsets are plainly there
    in the text, so there is no reason to infer them.
    """
    found = raw.dropna().astype(str).str.extract(_UTC_OFFSET.pattern, expand=False)
    clocks = sorted({o.strip() for o in found.dropna().unique()})
    if len(clocks) > 1:
        raise ValueError(
            f"the time column ({source}) mixes UTC offsets: {', '.join(clocks)}. "
            "The mainshock origin time is one instant on one clock, and this "
            "reader will not guess which of these it is on; put the catalogue on "
            "a single UTC offset, or a single local time, before reading it"
        )


def _read_dates(raw: pd.Series, dayfirst: bool, source: str) -> pd.Series:
    """Parse the assembled text, applying day-first only where it can apply.

    pandas' `dayfirst` is not the "only when the date is ambiguous" switch it is
    widely taken for, and this package's own documentation used to describe it that
    way. pandas guesses one format from the first row and swaps the day and month
    slots of *every* date, year-first stamps included: reading a KOERI export whose
    dates read 2023.02.06 with dayfirst=True made pandas settle on "%Y.%d.%m", so
    6 February became 2 June, every day past the 12th was dropped as unreadable,
    and a 3469-event sequence quietly became 1419 events with a different Omori
    decay. No error, no warning.

    So the flag is applied only to the rows where the day and the month really
    could be confused, and a row that states its year first is read as written.
    Refusing the combination outright was the alternative; it was rejected because
    the file in front of the user parses perfectly, and refusing to read a good
    file teaches nobody anything. A warning is raised instead, because a setting
    that silently does nothing is its own kind of lie.
    """
    if not dayfirst:
        return pd.to_datetime(raw, errors="coerce")
    year_first = raw.str.match(_YEAR_FIRST, na=False)
    if not year_first.any():
        return pd.to_datetime(raw, errors="coerce", dayfirst=True)
    example = raw[year_first].iloc[0]
    unaffected = int(year_first.sum())
    warnings.warn(
        f"dayfirst=True was asked for, but {unaffected} of {len(raw)} rows of "
        f"{source} state the year first, as in {example!r}. In that order the day "
        "and the month cannot be confused, so those rows were read exactly as "
        "written"
        + (
            "; the setting changed nothing"
            if unaffected == len(raw)
            else " and the flag was applied only to the rows that did not"
        ),
        stacklevel=2,
    )
    if year_first.all():
        return pd.to_datetime(raw, errors="coerce")
    as_written = pd.to_datetime(raw.where(year_first), errors="coerce")
    day_first = pd.to_datetime(raw.where(~year_first), errors="coerce", dayfirst=True)
    return as_written.where(year_first, day_first)


def _one_clock(t: pd.Series, source: str) -> pd.Series:
    """Drop a single UTC offset, and refuse a column that carries several.

    A UTC offset cannot be subtracted from a naive mainshock time. Both supported
    catalogues sit on one clock, so a uniform offset is dropped rather than
    shifted. When the offsets differ there is no such clock: pandas cannot build a
    datetime column at all and hands back a column of plain objects, which the old
    guard did not recognise, and the run died about sixty lines later with
    "unsupported operand type(s) for -: 'numpy.ndarray' and 'Timestamp'" - a
    message naming neither a column nor a timezone. Say what is wrong instead.
    """
    if isinstance(t.dtype, pd.DatetimeTZDtype):
        return t.dt.tz_localize(None)
    if pd.api.types.is_datetime64_any_dtype(t):
        return t
    clocks = sorted({_clock_name(v) for v in t if isinstance(v, pd.Timestamp)})
    if len(clocks) < 2:
        # Rebuilding a column of objects that all carry the same offset gives a
        # tz-aware column back, so it goes through the same stripping as a column
        # that arrived aware. Returning it as it came left the offset on, and the
        # run died on "Cannot subtract tz-naive and tz-aware datetime-like
        # objects" - the very error the rest of this function exists to prevent.
        rebuilt = pd.to_datetime(t, errors="coerce")
        if isinstance(rebuilt.dtype, pd.DatetimeTZDtype):
            return rebuilt.dt.tz_localize(None)
        return rebuilt
    raise ValueError(
        f"the time column ({source}) mixes clocks: {', '.join(clocks)}. The "
        "mainshock origin time is one instant on one clock, and this reader will "
        "not guess which of these it is on; put the catalogue on a single UTC "
        "offset, or a single local time, before reading it"
    )


def _clock_name(stamp: pd.Timestamp) -> str:
    """Name the clock a timestamp is on, for the message above."""
    offset = stamp.utcoffset()
    if offset is None:
        return "no offset"
    minutes = int(offset.total_seconds() // 60)
    sign = "+" if minutes >= 0 else "-"
    return f"UTC{sign}{abs(minutes) // 60:02d}:{abs(minutes) % 60:02d}"
