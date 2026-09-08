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
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd

from tremor_lab import __version__, constants
from tremor_lab.analysis import analyze_case
from tremor_lab.catalog import read_catalog

# The report is written in two columns: a label, then the value. Continuation lines
# line up under the value rather than under the label, so a block that runs over one
# line still reads as one entry.
_INDENT = " " * 21


def main(argv: list[str] | None = None) -> int:
    """Run the command-line tool. Returns the process exit status."""
    args = _build_parser().parse_args(argv)
    previous: dict[str, Any] = {}
    try:
        settings = _load_settings(args.settings)
        previous = _apply_constants(settings.get("constants", {}))
        result, source = _run(settings, args.settings.parent)
        # Built while the overrides are still in force, because the last line of
        # the report says which of them were changed from the published values.
        report = _report(result, source)
    except (OSError, KeyError, TypeError, ValueError, tomllib.TOMLDecodeError) as e:
        return _refuse(e)
    finally:
        for name, value in previous.items():
            setattr(constants, name, value)
    print(report)
    if not args.fmd_out:
        return 0
    # Everything below is the second half of what was asked for. It used to sit
    # outside the error handler above, so an output folder that did not exist, or a
    # CSV still open in Excel, ended the run with a Python traceback and exit 1 while
    # every other failure printed one line and exited 2. And when the window held no
    # events the table was skipped in silence at exit 0, which left whatever an
    # earlier run had put at that path in place for a plotting step to pick up as
    # though it belonged to this sequence.
    if result["fmd"] is None:
        return _refuse(
            f"no frequency-magnitude table was written to {args.fmd_out}: the window "
            f"holds no events, so there is no distribution to tabulate. Anything "
            f"already at that path is left over from an earlier run"
        )
    try:
        _write_fmd(result, args.fmd_out)
    except (OSError, ValueError) as e:
        return _refuse(
            f"the frequency-magnitude table could not be written to {args.fmd_out}: {e}"
        )
    print(f"\nfrequency-magnitude table written to {args.fmd_out}")
    return 0


def _refuse(problem) -> int:
    """Report one problem the way every failure in this tool is reported, and exit 2."""
    # KeyError formats itself with repr, which quotes the message and doubles
    # every backslash in a Windows path.
    message = (
        problem.args[0] if isinstance(problem, KeyError) and problem.args else problem
    )
    print(f"tremor-lab: {message}", file=sys.stderr)
    return 2


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
    for section in ("catalog", "mainshock"):
        if section not in settings:
            raise KeyError(f"{path}: missing [{section}] section")
    unknown = sorted(set(settings) - {"catalog", "mainshock", "analysis", "constants"})
    if unknown:
        raise KeyError(
            f"{path}: unknown section(s) {unknown}. A mistyped section name would "
            f"otherwise be ignored and the run would report different numbers"
        )
    for name, value in settings.get("constants", {}).items():
        _check_constant(path, name, value)
    for name, value in settings.get("analysis", {}).items():
        _check_analysis_value(path, name, value)
    return settings


def _apply_constants(overrides: dict[str, Any]) -> dict[str, Any]:
    """Put the overrides in force, and hand back what was there before.

    Nothing here decides whether they are allowed; `_load_settings` has already
    refused anything that is not a published constant or not a number. This only
    swaps the values in, so that a caller can swap them back.
    """
    previous = {name: getattr(constants, name) for name in overrides}
    for name, value in overrides.items():
        setattr(constants, name, value)
    return previous


def _check_constant(path: Path, name: str, value: Any) -> None:
    """Refuse an entry in [constants] that is not a published constant, or not a number.

    The gate here was once `hasattr`, which admits every attribute of the module and
    not only the published constants. `_PUBLISHED` is one of them: it is the snapshot
    the report compares against, so a settings file could set `_PUBLISHED = {}`,
    change a constant, and have the report state that the published defaults were
    used. Only the published names may be set.

    Nothing checked the value either. `DELTA_MB = true` ran to completion and printed
    a Bath expectation of M 6.80 instead of M 6.65, because Python's True is
    arithmetically 1; `DM = "0.1"` failed much later with a numpy casting message that
    named neither the constant nor the file. Every published constant is a number, so
    anything else is refused here, by name, before it can reach a formula.
    """
    if name not in constants._PUBLISHED:
        raise KeyError(
            f"{path}: unknown constant {name!r} in [constants]. The constants that "
            f"can be set here are {sorted(constants._PUBLISHED)}"
        )
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(
            f"{path}: [constants] {name} must be a number, but {value!r} was given. "
            f"Its published value is {constants._PUBLISHED[name]}"
        )


_ANALYSIS_NUMBERS = {
    "threshold",
    "window_days",
    "dm",
    "mc_correction",
    "n_boot",
    "seed",
    "min_events_for_mc",
    "min_events_for_fit",
    "n_fit_simulations",
}


def _check_analysis_value(path: Path, name: str, value: Any) -> None:
    """Refuse an [analysis] entry that is not the number it has to be.

    These are the same quantities as the [constants] overrides, reached by a
    different route, and they had no check at all. `mc_correction = true` ran to
    exit 0 and printed a completeness magnitude of 4.2 instead of 3.4, because
    Python's True is arithmetically 1; `dm = "0.1"` died inside numpy with a
    message about casting rules that named neither the setting nor the file.

    An unknown name is left alone here: `analyze_case` refuses it by name, and its
    suggestion of the nearest real keyword is better than anything this could say.
    """
    if name in _ANALYSIS_NUMBERS and (
        isinstance(value, bool) or not isinstance(value, (int, float))
    ):
        raise TypeError(
            f"{path}: [analysis] {name} must be a number, but {value!r} was given"
        )


def _run(settings: dict[str, Any], base: Path) -> tuple[dict[str, Any], Path]:
    """Load the configured catalogue and analyse it."""
    catalog_settings = dict(settings["catalog"])
    analysis = dict(settings.get("analysis", {}))
    mainshock = settings["mainshock"]
    columns = catalog_settings.pop("columns", {})
    path = base / catalog_settings.pop("path")
    dayfirst = catalog_settings.pop("dayfirst", False)
    radius_km = catalog_settings.pop("radius_km", None)
    if catalog_settings:
        raise KeyError(f"unknown keys in [catalog]: {sorted(catalog_settings)}")
    threshold = analysis.pop("threshold", None)
    window_days = analysis.pop("window_days", None)
    already_windowed = "dt_days" in columns

    if already_windowed:
        if radius_km is not None:
            raise KeyError(
                "radius_km cannot be applied to a catalogue configured with "
                "dt_days: that file carries elapsed days, not positions to "
                "measure from. Remove radius_km, or configure the raw export "
                "with date, time, lat and lon columns"
            )
        if columns.get("mag_type"):
            # No magnitude is converted on this route: `to_mw` is reached only
            # through `read_catalog`. Accepting the entry meant a file whose scale
            # column really did say Ms was analysed on raw Ms values while the
            # report named the column as though a Scordilis conversion had been
            # applied to them.
            raise KeyError(
                "mag_type cannot be applied to a catalogue configured with "
                "dt_days: that file is taken as it stands and no magnitude is "
                "converted to Mw on this route. Remove mag_type, or configure the "
                "raw export with date, time, lat and lon columns, where the scale "
                "column is read and the conversion is counted and reported"
            )
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
            radius_km=radius_km,
        )
    result = analyze_case(
        catalog,
        mainshock,
        mc_threshold=threshold,
        window_days=window_days,
        **analysis,
    )
    carried = (
        "rows_read",
        "rows_parsed",
        "radius_km",
        "removed_by_radius",
        "farthest_km",
    )
    result.update({k: v for k, v in catalog.attrs.items() if k in carried})
    result["seed"] = analysis.get("seed", 0)
    result["settings_used"] = _choices(
        analysis,
        window_days,
        threshold,
        dayfirst=dayfirst,
        already_windowed=already_windowed,
    )
    result["magnitude_note"] = _magnitude_note(
        columns,
        catalog.attrs.get("scale_counts"),
        already_windowed=already_windowed,
    )
    return result, path


def _choices(
    analysis: Mapping[str, Any],
    window_days,
    threshold,
    *,
    dayfirst: bool,
    already_windowed: bool,
) -> dict[str, Any]:
    """The value actually used for every choice that moves a number.

    Each of these can arrive by two routes that leave no trace in one another:
    `[analysis] mc_correction = 0.0` is a keyword argument passed to the estimator,
    while `[constants] MC_CORRECTION = 0.0` reassigns the module. The report used to
    list only the second. A run over a 30-day window with 0.5-wide bins and no
    completeness correction therefore ended with exactly the same closing line as the
    reference run, and a referee holding the two reports had no way to tell which
    choices had produced which numbers.

    Resolving each choice here - from the keyword where one was given, and from the
    constant otherwise - is what lets the report state the value that was used rather
    than the route it came by.
    """
    return {
        "window_days": constants.WINDOW_DAYS if window_days is None else window_days,
        "threshold": threshold,
        "dm": analysis.get("dm", constants.DM),
        "mc_correction": analysis.get("mc_correction", constants.MC_CORRECTION),
        "n_boot": analysis.get("n_boot", constants.N_BOOT),
        "n_fit_simulations": analysis.get(
            "n_fit_simulations", constants.N_FIT_SIMULATIONS
        ),
        "seed": analysis.get("seed", 0),
        # Not an [analysis] setting, but it decides how every timestamp is read
        # and so which events fall inside the window at all. A block that leaves
        # it out is not a complete account of a run on a raw catalogue.
        "dayfirst": dayfirst,
        "already_windowed": already_windowed,
    }


def _magnitude_note(
    columns: Mapping[str, str], scales, *, already_windowed: bool
) -> str:
    """State what happened to the magnitudes, from counts rather than intent."""
    source = f"column {columns.get('mag')!r}"
    if already_windowed:
        # A file of elapsed days and magnitudes is read straight through; no
        # conversion happens on this route, so the note must not name a scale
        # column, which would read as a homogenisation that never took place.
        return (
            f"{source}, used as published; a catalogue given as elapsed days is read "
            f"as it stands and no scale conversion is performed"
        )
    scale_col = columns.get("mag_type")
    if not scale_col:
        return f"{source}, used as published; no scale column given"
    # An empty or missing record of what was converted means nothing was, which
    # the next branch already says correctly; a separate branch for it only
    # created a second sentence for the same outcome.
    scales = scales or {}
    converted = scales.get("ms", 0) + scales.get("mb", 0)
    if not converted:
        return (
            f"{source}, scale column {scale_col!r} held no Ms or mb labels, so "
            f"nothing was converted and every magnitude was used as published"
        )
    return (
        f"{source}, scale column {scale_col!r}: {scales.get('ms', 0)} Ms and "
        f"{scales.get('mb', 0)} mb converted to Mw by the Scordilis relations, "
        f"{scales.get('mw', 0)} already Mw, "
        f"{scales.get('unconverted', 0)} used as published"
    )


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
    # This is where the distance limit is reported, so the settings block below does
    # not repeat it. The key is present whenever a spatial cut was possible at all;
    # the other route refuses radius_km outright, because it has no positions.
    if "radius_km" in result:
        radius, removed = result["radius_km"], result.get("removed_by_radius", 0)
        farthest = result.get("farthest_km")
        reach = f"; farthest kept {farthest:.0f} km" if farthest else ""
        lines.append(
            f"distance limit       none, every event in the window is included{reach}"
            if radius is None
            else f"distance limit       {radius} km from the epicentre, "
            f"{removed} events removed{reach}"
        )
    lines += [
        f"events in window     {result['n_events']}",
        f"completeness Mc      {_show(result['mc'])}  (maximum curvature)",
        f"threshold used       {_show(result['threshold'])}"
        f"  ({_show(result['n_above'])} events in its completeness class, "
        f"at or above {_show(result.get('sample_floor'))})",
    ]
    methods = result.get("mc_methods") or {}
    if len(methods) > 1:
        lines.append(
            "Mc, other methods    "
            + ", ".join(
                f"{name.replace('_', ' ')} {_show(value)}"
                for name, value in methods.items()
                if name != "maximum_curvature"
            )
        )
        spread = [v for v in methods.values() if v is not None]
        if len(spread) > 1 and max(spread) - min(spread) > 1e-9:
            lines.append(
                f"                     the methods span M {min(spread)} to "
                f"{max(spread)}; b depends on which is used"
            )
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
        test = result.get("omori_fit_test")
        if test is not None and test.p_value == test.p_value:
            verdict = fit_verdict(test)
            se = test.p_value_se if test.p_value_se == test.p_value_se else 0.0
            # A separate name from the bootstrap `spread` above: the two once
            # shared one, and the reader had to prove they never overlapped.
            se_text = f" +/- {se:.3f}" if se else ""
            lines.append(
                f"decay fit test       KS {test.statistic:.4f}, "
                f"p {test.p_value:.3f}{se_text}  ({verdict})"
            )
            lines.append(f"                     p from a {test.method}")
        if result.get("omori_warning"):
            lines.append(f"CAUTION              {result['omori_warning']}")
    lines += [
        f"mainshock energy     {result['mainshock_energy_j']:.3e} J",
        f"Bath expectation     M {result['bath_mag']:.2f}",
    ]
    if result.get("magnitude_note"):
        lines.append(f"magnitudes           {result['magnitude_note']}")
    settings = _settings_lines(result)
    lines.append(f"settings             {settings[0]}")
    lines += [f"{_INDENT}{line}" for line in settings[1:]]
    return "\n".join(lines)


def fit_verdict(test) -> str:
    """Say what a decay-fit p-value does and does not settle.

    Two things stop a bare comparison against 0.05 from meaning anything.

    The p-value may be uncalibrated. The Omori curve is fitted to the very times
    being tested, so it hugs them, and the textbook Kolmogorov p-value is far too
    generous: on the reference catalogue it reads 0.508 where the simulated null
    reads about 0.03. That number is reported, because hiding it would be worse,
    but no verdict may be read from it.

    And a simulated p-value is itself an estimate. When it sits within two of its
    own standard errors of the threshold, the comparison is decided by how many
    replicates happened to run, not by the catalogue, and the honest answer is
    that more replicates are needed.
    """
    se = test.p_value_se if test.p_value_se == test.p_value_se else 0.0
    if "uncalibrated" in test.method:
        return "no verdict: this p-value is uncalibrated"
    if test.p_value + 2 * se < 0.05:
        return "NOT consistent with a single Omori decay"
    if test.p_value - 2 * se > 0.05:
        return "consistent with one Omori decay"
    return "borderline; raise the replicate count to decide"


def _settings_lines(result: dict[str, Any]) -> list[str]:
    """Every choice that moved a number, as a block a methods section can quote.

    A printed number that does not carry its settings cannot be traced back to the
    run that produced it. This once printed only the constants that differ from their
    published values, which is a small part of the answer: the window, the bin width,
    the completeness correction, the resample and replicate counts and the seed are
    just as capable of moving every number in the report, and when they were set in
    [analysis] rather than [constants] none of them appeared anywhere. Two runs that
    shared nothing but a catalogue could close with the same line. Each choice is now
    stated as the value that was used, whichever route it arrived by, and the list of
    changed constants is kept as well, because it says what the rest of the session
    will do.
    """
    used = result["settings_used"]
    given = used["threshold"]
    if given is not None:
        threshold = f"threshold M {_number(given)}, given in the settings file"
    elif result.get("threshold") is not None:
        threshold = (
            f"threshold M {_number(result['threshold'])}, the completeness magnitude "
            f"estimated from the catalogue"
        )
    else:
        threshold = "no threshold: the window was too sparse to estimate one"
    changed = [
        f"{name}={getattr(constants, name)}"
        for name in sorted(constants._PUBLISHED)
        if getattr(constants, name) != constants._PUBLISHED[name]
    ]
    # Only the work that was actually done. A sparse window produces no decay fit
    # and no bootstrap, and printing the counts regardless described replicates
    # that were never drawn, in the block whose whole purpose is to say what
    # produced the numbers above it.
    if result.get("omori") is None:
        effort = "no decay fit was made, so no resamples or replicates were drawn"
    else:
        effort = (
            f"{used['n_boot']} bootstrap resamples"
            if used["n_boot"]
            else "no bootstrap: the uncertainty on p and c was not estimated"
        )
        effort += (
            f"; {used['n_fit_simulations']} fit-test replicates"
            if used["n_fit_simulations"]
            else "; no fit-test replicates, so the p-value is the uncalibrated one"
        )
    reading = (
        "elapsed days read from the file"
        if used.get("already_windowed")
        else (
            "dates read day-first where the order is ambiguous"
            if used.get("dayfirst")
            else "dates read month-first where the order is ambiguous"
        )
    )
    return [
        f"window {_number(used['window_days'])} days; "
        f"magnitude bins {_number(used['dm'])} wide; "
        f"Mc correction {used['mc_correction']:+g}",
        threshold,
        reading,
        f"{effort}; seed {used['seed']}",
        "constants at their published defaults"
        if not changed
        else "constants changed: " + ", ".join(changed),
    ]


def _number(value) -> str:
    """A number as a reader would write it: 180 rather than 180.0."""
    return f"{value:g}"


def _show(value) -> str:
    return "not estimated" if value is None else f"{value}"


def _write_fmd(result: dict[str, Any], path: Path) -> None:
    edges, inc, cum = result["fmd"]
    pd.DataFrame({"magnitude": edges, "incremental": inc, "cumulative": cum}).to_csv(
        path, index=False
    )
