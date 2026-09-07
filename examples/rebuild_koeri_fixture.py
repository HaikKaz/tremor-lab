"""Rebuild the KOERI-shaped fixture from the windowed one, and check it.

    python examples/rebuild_koeri_fixture.py

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
`tests/data/kahramanmaras_180d.csv` is a *derived* file: it carries elapsed days
and magnitudes, already windowed and already on the moment scale. Reading it back
exercises none of the catalogue pipeline, so the headline reproduction was
demonstrated on a path that skipped time parsing, magnitude homogenisation and
windowing entirely.

This script reconstructs the timestamps those elapsed days imply, and writes
`tests/data/kahramanmaras_180d_koeri.csv` in the shape a KOERI export actually
has: separate date and time columns, local clock, whole seconds, a magnitude
column and a scale column. `tests/test_regression.py` then runs the whole
pipeline over it and reproduces the locked values.

**It is a pipeline test, not provenance.** The timestamps are derived from the
elapsed days, so this file cannot corroborate the elapsed days themselves. The
original KOERI download is the only thing that can, and it should be archived
with the release. Do not present this file as the source data.

That it reproduces the locked values at whole-second precision is itself worth
knowing: the resolution an agency publishes is sufficient for these estimators.
"""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WINDOWED = ROOT / "tests" / "data" / "kahramanmaras_180d.csv"
KOERI = ROOT / "tests" / "data" / "kahramanmaras_180d_koeri.csv"

ORIGIN = "2023-02-06 01:17:32"
MAINSHOCK = {"t": ORIGIN, "lat": 37.1578, "lon": 36.8092, "mw": 7.8}
COLUMNS = {
    "date": "Tarih",
    "time": "Saat",
    "lat": "Enlem",
    "lon": "Boylam",
    "mag": "xM (Biggest Mag)",
}


# Events deliberately outside the window. Without these the row count and the
# in-window count would be the same number, and asserting one against the other
# would test nothing: every row would be inside the window by construction.
OUTSIDE_DAYS = [-30.0, -7.5, -1.0, -0.25, -0.001, 0.0, 180.0001, 181.0, 200.0, 365.0]


def build() -> pd.DataFrame:
    """Turn elapsed days back into the wall-clock stamps they imply."""
    windowed = pd.read_csv(WINDOWED)
    elapsed = list(windowed["dt_days"]) + OUTSIDE_DAYS
    mags = list(windowed["mw"]) + [4.0] * len(OUTSIDE_DAYS)
    lats = list(windowed["lat"]) + [37.2] * len(OUTSIDE_DAYS)
    lons = list(windowed["lon"]) + [37.0] * len(OUTSIDE_DAYS)
    stamps = pd.Timestamp(ORIGIN) + pd.to_timedelta(pd.Series(elapsed), unit="D")
    frame = pd.DataFrame(
        {
            "Tarih": stamps.dt.strftime("%Y.%m.%d"),
            "Saat": stamps.dt.strftime("%H:%M:%S"),
            "Enlem": lats,
            "Boylam": lons,
            "Depth(km)": 10.0,
            "xM (Biggest Mag)": mags,
            "Tip": "Mw",
        }
    )
    return frame.sort_values(["Tarih", "Saat"], kind="stable").reset_index(drop=True)


def main() -> int:
    from tremor_lab import analyze_case, read_catalog

    frame = build()
    KOERI.write_text(frame.to_csv(index=False), encoding="utf-8")
    print(
        f"wrote {KOERI.name}: {len(frame)} rows, KOERI shape, whole seconds, "
        f"of which {len(OUTSIDE_DAYS)} lie outside the window and must be excluded"
    )

    catalog = read_catalog(
        KOERI, COLUMNS, MAINSHOCK, window_days=180, mag_type_col="Tip"
    )
    result = analyze_case(
        catalog, MAINSHOCK, mc_threshold=3.5, window_days=180, n_boot=0
    )
    b, omori = result["b_value"], result["omori"]

    print("\nthrough the full pipeline: parse, homogenise, window, estimate")
    checks = [
        ("rows in the file", len(frame), 3469 + len(OUTSIDE_DAYS), 0),
        ("events in window", result["n_events"], 3469, 0),
        ("Mc", result["mc"], 3.4, 0),
        ("events above M 3.5", b.n, 1529, 0),
        ("b-value", round(b.b, 3), 0.844, 0),
        ("Omori p", round(omori.p, 3), 1.161, 0),
        ("Omori c", round(omori.c, 3), 0.497, 0),
        ("Omori k", round(omori.k), 359, 1),
    ]
    ok = True
    for name, got, want, tol in checks:
        matches = abs(got - want) <= tol
        ok &= matches
        verdict = "ok" if matches else "DIFFERS"
        print(f"  {name:<22}{got!s:>10}  want {want}   {verdict}")
    print("\n" + ("Reproduced from a raw-shaped catalogue." if ok else "DIFFERS."))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
