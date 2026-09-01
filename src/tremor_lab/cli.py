"""Command-line entry point: run one analysis from a settings file.

    tremor-lab run settings.toml

The settings file describes the catalogue, its column names, the mainshock, and
optionally any constant to override. It exists so that an analysis is a file that can
be archived alongside the paper rather than a sequence of arguments someone has to
remember.
"""

import argparse
import sys
import tomllib
from pathlib import Path
from typing import Any

import pandas as pd

from tremor_lab import __version__, constants
from tremor_lab.analysis import analyze_case
from tremor_lab.catalog import read_catalog


def main(argv: list[str] | None = None) -> int:
    """Run the command-line tool. Returns the process exit status."""
    args = _build_parser().parse_args(argv)
    try:
        settings = _load_settings(args.settings)
        result, source = _run(settings, args.settings.parent)
    except (OSError, KeyError, TypeError, ValueError, tomllib.TOMLDecodeError) as e:
        # KeyError formats itself with repr, which quotes the message and doubles
        # every backslash in a Windows path.
        message = e.args[0] if isinstance(e, KeyError) and e.args else e
        print(f"tremor-lab: {message}", file=sys.stderr)
        return 2
    print(_report(result, source))
    if args.fmd_out and result["fmd"] is not None:
        _write_fmd(result, args.fmd_out)
        print(f"\nfrequency-magnitude table written to {args.fmd_out}")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tremor-lab", description="Aftershock-sequence statistics."
    )
    parser.add_argument(
        "--version", action="version", version=f"Tremor Lab {__version__}"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="analyse one catalogue from a settings file")
    run.add_argument("settings", type=Path, help="TOML settings file")
    run.add_argument(
        "--fmd-out",
        type=Path,
        metavar="CSV",
        help="also write the frequency-magnitude table for the Gutenberg-Richter plot",
    )
    return parser


def _load_settings(path: Path) -> dict[str, Any]:
    with open(path, "rb") as handle:
        settings = tomllib.load(handle)
    for name, value in settings.get("constants", {}).items():
        if not hasattr(constants, name):
            raise KeyError(f"{path}: unknown constant {name!r} in [constants]")
        setattr(constants, name, value)
    for section in ("catalog", "mainshock"):
        if section not in settings:
            raise KeyError(f"{path}: missing [{section}] section")
    return settings


def _run(settings: dict[str, Any], base: Path) -> tuple[dict[str, Any], Path]:
    """Load the configured catalogue and analyse it."""
    catalog_settings = dict(settings["catalog"])
    analysis = dict(settings.get("analysis", {}))
    mainshock = settings["mainshock"]
    columns = catalog_settings.pop("columns", {})
    path = base / catalog_settings.pop("path")
    dayfirst = catalog_settings.pop("dayfirst", False)
    if catalog_settings:
        raise KeyError(f"unknown keys in [catalog]: {sorted(catalog_settings)}")
    threshold = analysis.pop("threshold", None)
    window_days = analysis.pop("window_days", None)

    if "dt_days" in columns:
        # Already windowed: elapsed times are in the file, so there is nothing to parse
        # and no mainshock position to measure against.
        frame = pd.read_csv(path)
        catalog = frame.rename(
            columns={columns["dt_days"]: "dt_days", columns["mag"]: "mw"}
        )
    else:
        catalog = read_catalog(
            path,
            columns,
            mainshock,
            window_days=window_days,
            mag_type_col=columns.get("mag_type"),
            dayfirst=dayfirst,
        )
    result = analyze_case(
        catalog,
        mainshock,
        mc_threshold=threshold,
        window_days=window_days,
        **analysis,
    )
    result.update({k: v for k, v in catalog.attrs.items() if k.startswith("rows_")})
    scale_col = columns.get("mag_type")
    result["magnitude_note"] = (
        f"column {columns.get('mag')!r}, scale column {scale_col!r}, converted to "
        f"Mw by the Scordilis relations where the scale is Ms or mb"
        if scale_col
        else f"column {columns.get('mag')!r}, used as published; no scale column "
        f"given, so no homogenisation was applied"
    )
    return result, path


def _report(result: dict[str, Any], source: Path) -> str:
    """Render the analysis as the block of text a reader would paste into a table."""
    read = result.get("rows_read")
    lines = [
        f"Tremor Lab {__version__}",
        f"catalogue            {source.name}",
    ]
    if read is not None and result.get("rows_parsed") is not None:
        dropped = read - result["rows_parsed"]
        lines.append(
            f"rows                 {read} read, {result['rows_parsed']} usable"
            + (f", {dropped} unreadable and left out" if dropped else "")
        )
    if result.get("n_unusable"):
        lines.append(f"incomplete rows      {result['n_unusable']} dropped")
    lines += [
        f"events in window     {result['n_events']}",
        f"completeness Mc      {_show(result['mc'])}",
        f"threshold used       {_show(result['threshold'])}"
        f"  ({_show(result['n_above'])} events at or above it)",
    ]
    b, omori, spread = result["b_value"], result["omori"], result["omori_bootstrap"]
    if b is None or omori is None:
        lines.append(f"not estimated        {result['note']}")
    else:
        lines += [
            f"b-value              {b.b:.3f} +/- {b.sigma:.3f}  on {b.n} events",
            f"Omori p              {omori.p:.3f}"
            + (f" +/- {spread.p_std:.3f}" if spread else ""),
            f"Omori c              {omori.c:.3f}"
            + (f" +/- {spread.c_std:.3f}" if spread else ""),
            f"Omori k              {omori.k:.1f}",
        ]
        if result.get("omori_warning"):
            lines.append(f"CAUTION              {result['omori_warning']}")
    lines += [
        f"mainshock energy     {result['mainshock_energy_j']:.3e} J",
        f"Bath expectation     M {result['bath_mag']:.2f}",
    ]
    if result.get("magnitude_note"):
        lines.append(f"magnitudes           {result['magnitude_note']}")
    return "\n".join(lines)


def _show(value) -> str:
    return "not estimated" if value is None else f"{value}"


def _write_fmd(result: dict[str, Any], path: Path) -> None:
    edges, inc, cum = result["fmd"]
    pd.DataFrame({"magnitude": edges, "incremental": inc, "cumulative": cum}).to_csv(
        path, index=False
    )
