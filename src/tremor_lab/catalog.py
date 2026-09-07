"""Reading earthquake catalogues and reducing them to an aftershock window.

The analysis-ready table has the columns named in `CATALOG_COLUMNS`; `analyze_case`
consumes exactly that shape, so the reader and the analysis cannot drift apart.
"""

from collections.abc import Mapping

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray

from tremor_lab import constants
from tremor_lab.magnitude import to_mw

CATALOG_COLUMNS = ("t", "lat", "lon", "mw", "dt_days", "dist_km")


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
        Must carry t (datetime), lat, lon and mw.
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
        Columns as in `CATALOG_COLUMNS`, index reset.
    """
    window_days = constants.WINDOW_DAYS if window_days is None else window_days
    out = catalog.dropna(subset=["t", "lat", "lon", "mw"]).copy()
    origin = pd.Timestamp(mainshock["t"])
    out["dt_days"] = (out["t"] - origin).dt.total_seconds() / 86400.0
    out["dist_km"] = haversine_km(
        mainshock["lat"], mainshock["lon"], out["lat"].to_numpy(), out["lon"].to_numpy()
    )
    inside = (out["dt_days"] > 0) & (out["dt_days"] <= window_days)
    in_time = int(inside.sum())
    if radius_km is not None:
        inside &= out["dist_km"] <= radius_km
    kept = out.loc[inside, list(CATALOG_COLUMNS)].reset_index(drop=True)
    kept.attrs = {
        "radius_km": radius_km,
        "removed_by_radius": in_time - len(kept),
        "farthest_km": float(kept["dist_km"].max()) if len(kept) else None,
    }
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
    caller can report how much of the file survived.

    Parameters
    ----------
    path : str or path-like
        CSV file.
    columns : mapping of str to str
        Column names for the roles described above.
    mainshock : mapping
        Keys t, lat, lon, mw, with t on the same clock as the catalogue.
    window_days : float, optional
        Aftershock window length in days. Defaults to `constants.WINDOW_DAYS`.
    mag_type_col : str, optional
        Column of magnitude-scale labels. When given, magnitudes are homogenised to Mw
        by `tremor_lab.magnitude.to_mw`; when omitted they are taken as already on the
        moment scale.
    dayfirst : bool, optional
        Interpret ambiguous dates as day-first. False by default. Set it for European
        exports such as 06.02.2023: read as month-first, that date silently becomes
        2 June rather than 6 February.
    radius_km : float, optional
        Great-circle distance limit from the mainshock epicentre. No limit by
        default; see `window`.

    Returns
    -------
    DataFrame
        Columns as in `CATALOG_COLUMNS`.
    """
    frame = pd.read_csv(path)
    parsed = pd.DataFrame(
        {
            "t": _parse_times(frame, columns, dayfirst),
            "lat": pd.to_numeric(frame[columns["lat"]], errors="coerce"),
            "lon": pd.to_numeric(frame[columns["lon"]], errors="coerce"),
        }
    )
    mag = pd.to_numeric(frame[columns["mag"]], errors="coerce")
    scales: dict[str, int] = {}
    if mag_type_col:
        labels = frame[mag_type_col].to_numpy(dtype=str)
        parsed["mw"] = to_mw(mag.to_numpy(), labels)
        scales = _scale_counts(labels)
    else:
        parsed["mw"] = mag

    out = window(parsed, mainshock, window_days, radius_km)
    out.attrs = {
        **out.attrs,
        "rows_read": len(frame),
        "rows_parsed": int(parsed.notna().all(axis=1).sum()),
        "events_in_window": len(out),
        "scale_counts": scales,
    }
    return out


def _scale_counts(labels) -> dict[str, int]:
    """How many magnitudes each conversion branch of `to_mw` actually took.

    A scale column is not always a scale column: KOERI's "Type" holds
    "Earthquake". Counting the branches lets a caller report what happened
    instead of asserting that a conversion occurred.
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


def _parse_times(
    frame: pd.DataFrame, columns: Mapping[str, str], dayfirst: bool
) -> pd.Series:
    """Combine the configured date and time columns into one naive datetime series."""
    if "datetime" in columns:
        raw = frame[columns["datetime"]].astype(str)
    else:
        raw = (
            frame[columns["date"]].astype(str).str.strip()
            + " "
            + frame[columns["time"]].astype(str).str.strip()
        )
    t = pd.to_datetime(raw, errors="coerce", dayfirst=dayfirst)
    # A UTC offset cannot be subtracted from a naive mainshock time. Both catalogues
    # sit on a single clock, so the offset is dropped rather than shifted.
    if isinstance(t.dtype, pd.DatetimeTZDtype):
        t = t.dt.tz_localize(None)
    return t
