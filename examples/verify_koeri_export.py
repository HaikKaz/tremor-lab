"""Check a raw KOERI export against the bundled fixture, then look below M 3.0.

    python examples/verify_koeri_export.py slice1.txt slice2.txt ...

WHY THIS EXISTS
---------------
`tests/data/kahramanmaras_180d.csv` is derived: it carries elapsed days, magnitudes
and coordinates, but the query that produced it was never written down. The
completeness magnitude it yields rests on four bins sitting immediately above its
M 3.0 truncation floor, so the open question is what the catalogue looks like below
that floor, and answering it needs the original, untruncated export.

An untruncated export is only useful if it selects the same events as the fixture
above M 3.0. If the spatial or temporal selection differs, any completeness estimate
from the events below M 3.0 describes a different catalogue and cannot be compared
with anything already reported. So this script does not merely count rows. It
matches the export against the fixture event by event and reports how well the two
agree, before reporting anything about the deeper population.

WHAT A KOERI EXPORT LOOKS LIKE
------------------------------
Tab-separated, with a one-line header carrying these columns:

    No, Event ID, Date, Origin Time, Latitude, Longitude, Depth(km),
    xM, MD, ML, Mw, Ms, Mb, Type, Location

Date is `2023.02.06`, Origin Time is `04:17:32.20`, and `xM` is the largest
reported magnitude, which is the column the fixture uses. Place names carry Turkish
characters in a non-UTF-8 encoding.

**Times are local (UTC+3), while the mainshock time in `examples/kahramanmaras.toml`
is UTC**, so a naive match is three hours out. The script tries several offsets and
reports which one aligns rather than assuming; on the fixture the answer is +3 h,
which is how that convention was established in the first place.

THE QUERY TO RUN
----------------
KOERI caps a query at 50,000 records and returns them newest first, so a wide date
range silently comes back covering only the most recent months. Ask in slices, each
well under the cap, and pass them all to this script. The fixture's own footprint is

    lat 33.7 to 40.6,  lon 33.0 to 41.5

so a box of lat 33.5 to 41.0 and lon 32.5 to 42.0 contains it with margin. The
window runs from 6 February 2023 to 5 August 2023.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "data" / "kahramanmaras_180d.csv"

ORIGIN_UTC = pd.Timestamp("2023-02-06 01:17:32")
WINDOW_DAYS = 180.0
FIXTURE_N = 3469
FIXTURE_FLOOR = 3.0
KOERI_CAP = 50000

# Matching tolerances. Coordinates are published to four decimals and magnitudes to
# one, so these are loose enough for rounding and tight enough to mean something.
TOL_DAYS = 2.0 / 86400.0
TOL_DEG = 0.002
TOL_MAG = 0.051

REQUIRED = ("Date", "Origin Time", "Latitude", "Longitude", "xM")


def read_export(paths: list[Path]) -> pd.DataFrame:
    """Read and concatenate KOERI exports, dropping duplicates by event id."""
    frames = []
    for p in paths:
        # KOERI does not export UTF-8. cp1254 is the Turkish code page; latin-1
        # accepts any byte, so it is the backstop that cannot fail.
        df = None
        for enc in ("utf-8", "cp1254", "latin-1"):
            try:
                df = pd.read_csv(p, sep="\t", dtype=str, encoding=enc)
                break
            except UnicodeDecodeError:
                continue
        if df is None:
            raise SystemExit(f"{p.name}: could not decode in any known encoding")
        df = df.rename(columns=lambda c: c.strip())
        missing = sorted(set(REQUIRED) - set(df.columns))
        if missing:
            raise SystemExit(f"{p.name}: missing expected columns {missing}")
        if len(df) >= KOERI_CAP:
            print(
                f"  WARNING  {p.name} has {len(df)} rows, at or above the\n"
                f"           {KOERI_CAP} cap. It is probably truncated;"
                f" split this slice."
            )
        frames.append(df)

    raw = pd.concat(frames, ignore_index=True)
    before = len(raw)
    if "Event ID" in raw.columns:
        raw = raw.drop_duplicates(subset="Event ID")
    if before != len(raw):
        print(f"  dropped {before - len(raw)} duplicate rows across slices")

    out = pd.DataFrame()
    out["t"] = pd.to_datetime(
        raw["Date"].str.strip() + " " + raw["Origin Time"].str.strip(),
        format="%Y.%m.%d %H:%M:%S.%f",
        errors="coerce",
    )
    out["lat"] = pd.to_numeric(raw["Latitude"], errors="coerce")
    out["lon"] = pd.to_numeric(raw["Longitude"], errors="coerce")
    out["mag"] = pd.to_numeric(raw["xM"], errors="coerce")
    bad = int(out.isna().any(axis=1).sum())
    if bad:
        print(f"  dropped {bad} rows that would not parse")
    return out.dropna().sort_values("t").reset_index(drop=True)


def elapsed_days(export: pd.DataFrame, offset_h: float) -> pd.Series:
    """Days since the mainshock, reading the export clock as UTC plus offset_h."""
    shift = ORIGIN_UTC + pd.Timedelta(offset_h, unit="h")
    return (export["t"] - shift).dt.total_seconds() / 86400.0


def match_count(export: pd.DataFrame, fixture: pd.DataFrame, offset_h: float) -> int:
    """How many fixture events the export reproduces at this clock offset."""
    cand = export.assign(dt_days=elapsed_days(export, offset_h))
    cand = cand[
        (cand.dt_days >= -TOL_DAYS)
        & (cand.dt_days <= WINDOW_DAYS)
        & (cand.mag >= FIXTURE_FLOOR - TOL_MAG)
    ]
    if cand.empty:
        return 0
    # Match on time first, the most discriminating field, then confirm the rest.
    times = cand.dt_days.to_numpy()
    order = np.argsort(times)
    sorted_times = times[order]
    mags = cand.mag.to_numpy()[order]
    lats = cand.lat.to_numpy()[order]
    lons = cand.lon.to_numpy()[order]

    hits = 0
    for ft, fm, fla, flo in fixture[["dt_days", "mw", "lat", "lon"]].to_numpy():
        lo = np.searchsorted(sorted_times, ft - TOL_DAYS, "left")
        hi = np.searchsorted(sorted_times, ft + TOL_DAYS, "right")
        if hi <= lo:
            continue
        near = (
            (np.abs(mags[lo:hi] - fm) <= TOL_MAG)
            & (np.abs(lats[lo:hi] - fla) <= TOL_DEG)
            & (np.abs(lons[lo:hi] - flo) <= TOL_DEG)
        )
        if near.any():
            hits += 1
    return hits


def report_depth(export: pd.DataFrame, offset_h: float) -> None:
    """What the export holds below the fixture's M 3.0 floor."""
    win = export.assign(dt_days=elapsed_days(export, offset_h))
    win = win[(win.dt_days >= 0) & (win.dt_days <= WINDOW_DAYS)]
    deep = win[win.mag < FIXTURE_FLOOR - TOL_MAG]
    print(
        f"\nin the 180-day window: {len(win)} events, of which {len(deep)} below M 3.0"
    )
    if not len(deep):
        return
    print("  incremental counts by 0.1 bin, from the new floor upward:")
    lo = float(np.floor(win.mag.min() * 10) / 10)
    for edge in np.round(np.arange(lo, lo + 0.8, 0.1), 1):
        n = int(((win.mag >= edge - 0.05) & (win.mag < edge + 0.05)).sum())
        print(f"    M {edge:.1f}  {n}")
    print("\n  Now rerun the analysis on this file. If maximum curvature still")
    print("  returns 3.4 the truncation was not driving it; if it moves,")
    print("  the caveat in section 9e becomes a measured result.")


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    paths = [Path(a) for a in argv]
    for p in paths:
        if not p.exists():
            raise SystemExit(f"no such file: {p}")

    fixture = pd.read_csv(FIXTURE)
    print(
        f"fixture: {len(fixture)} events, M {fixture.mw.min()} to {fixture.mw.max()}\n"
    )

    print(f"reading {len(paths)} export file(s)")
    export = read_export(paths)
    print(f"  {len(export)} events, {export.t.min()} to {export.t.max()}")
    print(f"  magnitudes {export.mag.min()} to {export.mag.max()}")

    if export.mag.min() >= FIXTURE_FLOOR - TOL_MAG:
        print(
            "\n  WARNING  this export goes no deeper than the fixture does, so it\n"
            "           cannot answer the question below M 3.0. Re-query at\n"
            "           minimum magnitude 0."
        )

    # A correct export in local time starts minutes after the mainshock, so
    # "starts later than the origin" is normal and must not trip this. What is
    # fatal is missing the beginning of the sequence, which is what the record cap
    # produces: it comes back years adrift.
    if export.t.max() < ORIGIN_UTC or export.t.min() > ORIGIN_UTC + pd.Timedelta(
        1, unit="D"
    ):
        print(
            f"\n  STOP  this export spans {export.t.min()} to {export.t.max()},\n"
            f"        which misses the first day after the mainshock at"
            f" {ORIGIN_UTC} UTC.\n"
            "        The start of the sequence is where completeness matters most.\n"
            "        This is what the record cap returned newest-first looks like:\n"
            "        ask for a narrower date range."
        )
        return 1

    print("\nclock offset (KOERI publishes local time, the mainshock above is UTC):")
    best, best_hits = 0.0, -1
    for off in (0.0, 3.0, -3.0):
        hits = match_count(export, fixture, off)
        flag = ""
        if hits > best_hits:
            best, best_hits = off, hits
            flag = "  <-- best"
        print(f"  offset {off:+.0f} h : {hits:5d} / {FIXTURE_N} matched{flag}")

    pct = 100.0 * best_hits / FIXTURE_N
    print(
        f"\nbest match: {best_hits} of {FIXTURE_N} ({pct:.1f}%) at offset {best:+.0f} h"
    )
    if pct >= 99.0:
        print("  Good. The export selects essentially the same events as the")
        print("  fixture, so what it says below M 3.0 describes that catalogue.")
    elif pct >= 90.0:
        print("  Usable but imperfect. Record the exact query you ran alongside")
        print("  any completeness number from it, and say the match was inexact.")
    else:
        print("  NOT usable as a drop-in. The selection differs from the")
        print("  fixture's, so a completeness estimate from it describes a")
        print("  different catalogue. Widen the box (the fixture spans")
        print("  lat 33.7-40.6, lon 33.0-41.5) and check the date slices cover")
        print("  2023-02-06 to 2023-08-05 without hitting the record cap.")

    report_depth(export, best)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
