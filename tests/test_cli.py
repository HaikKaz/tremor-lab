"""The command-line tool: it must report the reference values and fail loudly."""

from pathlib import Path

import pandas as pd
import pytest

from tremor_lab import constants
from tremor_lab.cli import main

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"
EXAMPLE = ROOT / "examples" / "kahramanmaras.toml"

# TOML reads a backslash inside a double-quoted string as an escape character, so
# Windows paths are written with forward slashes.
FIXTURE = (DATA / "kahramanmaras_180d.csv").as_posix()
RAW_KOERI = (DATA / "raw_catalog_koeri.csv").as_posix()


@pytest.fixture(autouse=True)
def restore_constants():
    """A settings file may reassign constants; keep that inside the test that did it."""
    saved = {k: v for k, v in vars(constants).items() if k.isupper()}
    yield
    for name, value in saved.items():
        setattr(constants, name, value)


def write_settings(path, body):
    path.write_text(body)
    return path


def test_the_bundled_example_reports_the_reference_values(capsys):
    assert main(["run", str(EXAMPLE)]) == 0
    out = capsys.readouterr().out
    assert "events in window     3469" in out
    assert "completeness Mc      3.4" in out
    assert "b-value              0.844 +/- 0.019  on 1529 events" in out
    assert "Omori p              1.161" in out
    assert "Omori c              0.497" in out
    assert "Omori k              358.9" in out


def test_version_is_reported(capsys):
    with pytest.raises(SystemExit) as exit_status:
        main(["--version"])
    assert exit_status.value.code == 0
    assert "Tremor Lab 1.0.0" in capsys.readouterr().out


def test_a_constant_set_in_the_settings_file_changes_the_result(tmp_path, capsys):
    settings = write_settings(
        tmp_path / "no_correction.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n"
        "[analysis]\nn_boot = 0\n"
        "[constants]\nMC_CORRECTION = 0.0\n",
    )
    assert main(["run", str(settings)]) == 0
    # Without the +0.2 correction the completeness magnitude is the modal bin itself.
    assert "completeness Mc      3.2" in capsys.readouterr().out


def test_a_raw_catalogue_is_read_windowed_and_analysed(tmp_path, capsys):
    settings = write_settings(
        tmp_path / "raw.toml",
        f'[catalog]\npath = "{RAW_KOERI}"\n'
        '[catalog.columns]\ndate = "Tarih"\ntime = "Saat"\nlat = "Enlem"\n'
        'lon = "Boylam"\nmag = "Mag"\nmag_type = "Tip"\n'
        '[mainshock]\nt = "2023-02-06 01:17:32"\nlat = 37.0\nlon = 37.0\nmw = 7.8\n'
        "[analysis]\nwindow_days = 180\n",
    )
    assert main(["run", str(settings)]) == 0
    out = capsys.readouterr().out
    assert "events in window     4" in out
    # Four events cannot support an estimate, and the tool says so rather than fitting.
    assert "needs more than 50" in out


def test_the_frequency_magnitude_table_can_be_written(tmp_path, capsys):
    out_csv = tmp_path / "fmd.csv"
    assert main(["run", str(EXAMPLE), "--fmd-out", out_csv.as_posix()]) == 0
    table = pd.read_csv(out_csv)
    assert list(table.columns) == ["magnitude", "incremental", "cumulative"]
    assert table["incremental"].sum() == 3469
    assert f"written to {out_csv}" in capsys.readouterr().out


def test_an_unknown_constant_is_refused(tmp_path, capsys):
    settings = write_settings(
        tmp_path / "typo.toml",
        '[catalog]\npath = "x.csv"\n[mainshock]\nmw = 7.8\n'
        "[constants]\nDELTA_MB_TYPO = 1.2\n",
    )
    assert main(["run", str(settings)]) == 2
    assert "unknown constant 'DELTA_MB_TYPO'" in capsys.readouterr().err


def test_an_unknown_catalogue_key_is_refused(tmp_path, capsys):
    settings = write_settings(
        tmp_path / "stray.toml",
        f'[catalog]\npath = "{FIXTURE}"\nwindow_days = 180\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n",
    )
    # window_days belongs in [analysis]; silently ignoring it would misreport the run.
    assert main(["run", str(settings)]) == 2
    assert "unknown keys in [catalog]: ['window_days']" in capsys.readouterr().err


def test_a_missing_section_is_refused(tmp_path, capsys):
    settings = write_settings(tmp_path / "bare.toml", '[catalog]\npath = "x.csv"\n')
    assert main(["run", str(settings)]) == 2
    assert "missing [mainshock] section" in capsys.readouterr().err


def test_a_missing_settings_file_is_refused(tmp_path, capsys):
    assert main(["run", str(tmp_path / "absent.toml")]) == 2
    assert "tremor-lab:" in capsys.readouterr().err


def test_malformed_settings_are_refused(tmp_path, capsys):
    settings = write_settings(tmp_path / "broken.toml", "[catalog\npath = 'x.csv'\n")
    assert main(["run", str(settings)]) == 2
    assert "tremor-lab:" in capsys.readouterr().err


def test_a_shorter_window_in_the_settings_file_changes_the_result(tmp_path, capsys):
    # Without this the CLI could stop forwarding window_days and the suite would
    # not notice, because every other CLI test uses the default 180.
    settings = write_settings(
        tmp_path / "w30.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n"
        "[analysis]\nwindow_days = 30\nthreshold = 3.5\nn_boot = 0\n",
    )
    assert main(["run", str(settings)]) == 0
    out = capsys.readouterr().out
    assert "events in window     2501" in out
    assert "b-value              0.819 +/- 0.020  on 1207 events" in out


# ------------------------------------------------------- the provenance block
# The block the tool prints is what gets pasted into a methods section, so every
# choice that moved a number has to be in it. It used to list only the constants
# that differed from their published values, and a setting made in [analysis]
# reached the estimators without appearing anywhere in the report.


def test_a_choice_made_in_the_analysis_section_is_reported_as_the_value_used(
    tmp_path, capsys
):
    settings = write_settings(
        tmp_path / "elsewhere.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n"
        "[analysis]\nwindow_days = 30\ndm = 0.5\nmc_correction = 0.0\n"
        "n_boot = 0\nn_fit_simulations = 0\nseed = 7\n",
    )
    assert main(["run", str(settings)]) == 0
    out = capsys.readouterr().out
    # None of these went through [constants], so none of them used to be printed.
    assert "window 30 days" in out
    assert "magnitude bins 0.5 wide" in out
    assert "Mc correction +0" in out
    # Nothing was resampled and nothing was simulated, and the block says so in
    # those words rather than reporting counts of zero as though work had been
    # done and produced nothing.
    assert "no bootstrap: the uncertainty on p and c was not estimated" in out
    assert "no fit-test replicates, so the p-value is the uncalibrated one" in out
    assert "seed 7" in out
    # Nothing was set in [constants], so the correction is reported as the value
    # used without the constant ever having been touched.
    assert "MC_CORRECTION" not in out


def test_the_reference_run_reports_the_choices_that_produced_its_numbers(capsys):
    assert main(["run", str(EXAMPLE)]) == 0
    out = capsys.readouterr().out
    assert "window 180 days; magnitude bins 0.1 wide; Mc correction +0.2" in out
    assert "200 bootstrap resamples" in out
    # The replicate count is whatever is in force, not what the file happens to
    # name: the test session itself lowers this one to keep the suite quick.
    assert f"{constants.N_FIT_SIMULATIONS} fit-test replicates; seed 0" in out


def test_the_report_says_whether_the_threshold_was_given_or_estimated(tmp_path, capsys):
    # A b-value at M 3.5 means one thing when 3.5 was estimated from the catalogue
    # and another when it was typed into the settings file; the number is the same.
    assert main(["run", str(EXAMPLE)]) == 0
    assert "threshold M 3.5, given in the settings file" in capsys.readouterr().out
    settings = write_settings(
        tmp_path / "estimated.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n"
        "[analysis]\nn_boot = 0\nn_fit_simulations = 0\n",
    )
    assert main(["run", str(settings)]) == 0
    assert "estimated from the catalogue" in capsys.readouterr().out


def test_a_constant_changed_from_its_published_value_is_still_named(tmp_path, capsys):
    settings = write_settings(
        tmp_path / "changed.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n"
        "[analysis]\nn_boot = 0\nn_fit_simulations = 0\n"
        "[constants]\nMC_CORRECTION = 0.0\n",
    )
    assert main(["run", str(settings)]) == 0
    out = capsys.readouterr().out
    assert "constants changed: MC_CORRECTION=0.0" in out
    # The same value also has to appear as the correction that was actually applied,
    # so the two routes into the estimator read alike in the report.
    assert "Mc correction +0" in out


# ------------------------------------------------------------ what a file may set


def test_the_published_snapshot_cannot_be_set_from_the_settings_file(tmp_path, capsys):
    # _PUBLISHED is the record the report compares against. While the gate was
    # hasattr, a file could empty it and then change a constant, and the report
    # would state that the published defaults had been used.
    settings = write_settings(
        tmp_path / "hijack.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n"
        "[constants]\n_PUBLISHED = {}\nMC_CORRECTION = 0.0\n",
    )
    assert main(["run", str(settings)]) == 2
    assert "unknown constant '_PUBLISHED'" in capsys.readouterr().err
    assert len(constants._PUBLISHED) > 0


def test_a_constant_that_is_not_a_number_is_refused_by_name(tmp_path, capsys):
    # DELTA_MB = true ran to completion and published a Bath expectation of M 6.80,
    # because Python's True is arithmetically 1.
    settings = write_settings(
        tmp_path / "boolean.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n"
        "[constants]\nDELTA_MB = true\n",
    )
    assert main(["run", str(settings)]) == 2
    err = capsys.readouterr().err
    assert "DELTA_MB must be a number" in err
    assert "True" in err


def test_a_constant_given_as_text_names_the_constant_rather_than_failing_in_numpy(
    tmp_path, capsys
):
    settings = write_settings(
        tmp_path / "text.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n"
        '[constants]\nDM = "0.1"\n',
    )
    assert main(["run", str(settings)]) == 2
    err = capsys.readouterr().err
    assert "DM must be a number" in err
    assert "casting rule" not in err


def test_a_refused_settings_file_leaves_the_constants_alone(tmp_path, capsys):
    # The constants were applied before the file was validated, so a file the tool
    # refused still steered the next analysis run in the same session.
    settings = write_settings(
        tmp_path / "no_mainshock.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        "[constants]\nMC_CORRECTION = 1.0\nDELTA_MB = 2.0\n",
    )
    assert main(["run", str(settings)]) == 2
    assert "missing [mainshock] section" in capsys.readouterr().err
    assert constants.MC_CORRECTION == constants._PUBLISHED["MC_CORRECTION"]
    assert constants.DELTA_MB == constants._PUBLISHED["DELTA_MB"]


# --------------------------------------------------------- the magnitude note


def test_a_scale_column_is_refused_where_no_conversion_can_happen(tmp_path, capsys):
    # A catalogue of elapsed days is read straight through; to_mw is never called on
    # that route. Accepting mag_type meant a file of genuine Ms values was analysed
    # unconverted while the report named the scale column as though it had been used.
    settings = write_settings(
        tmp_path / "scale.toml",
        f'[catalog]\npath = "{FIXTURE}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\nmag_type = "Tip"\n'
        "[mainshock]\nmw = 7.8\n"
        "[analysis]\nn_boot = 0\n",
    )
    assert main(["run", str(settings)]) == 2
    assert "mag_type cannot be applied" in capsys.readouterr().err


def test_the_report_claims_no_conversion_on_a_file_that_is_already_windowed(capsys):
    assert main(["run", str(EXAMPLE)]) == 0
    out = capsys.readouterr().out
    assert "no scale conversion is performed" in out
    assert "Scordilis" not in out


# ------------------------------------------------ the frequency-magnitude table


def test_an_unwritable_table_target_is_refused_like_every_other_failure(
    tmp_path, capsys
):
    # A directory stands in for the ordinary cases: an output folder that does not
    # exist, or a CSV still open in Excel. This escaped as a raw traceback and exit 1
    # because the write sat outside the error handler.
    assert main(["run", str(EXAMPLE), "--fmd-out", tmp_path.as_posix()]) == 2
    err = capsys.readouterr().err
    assert err.startswith("tremor-lab: ")
    assert "could not be written" in err


def test_an_empty_window_says_the_table_was_not_written(tmp_path, capsys):
    catalog = tmp_path / "late.csv"
    catalog.write_text("dt_days,mw\n200,4.1\n250,3.9\n")
    out_csv = tmp_path / "fmd.csv"
    out_csv.write_text("magnitude,incremental,cumulative\n3.4,999,999\n")
    settings = write_settings(
        tmp_path / "empty.toml",
        f'[catalog]\npath = "{catalog.as_posix()}"\n'
        '[catalog.columns]\ndt_days = "dt_days"\nmag = "mw"\n'
        "[mainshock]\nmw = 7.8\n",
    )
    assert main(["run", str(settings), "--fmd-out", out_csv.as_posix()]) == 2
    err = capsys.readouterr().err
    assert "no frequency-magnitude table was written" in err
    # The stale file from an earlier run is still there, and the message says so,
    # because a plotting step pointed at that path would otherwise draw it as if it
    # belonged to this sequence.
    assert "999" in out_csv.read_text()


def test_the_report_says_how_much_of_the_file_was_usable(tmp_path, capsys):
    settings = write_settings(
        tmp_path / "raw.toml",
        f'[catalog]\npath = "{RAW_KOERI}"\n'
        '[catalog.columns]\ndate = "Tarih"\ntime = "Saat"\nlat = "Enlem"\n'
        'lon = "Boylam"\nmag = "Mag"\nmag_type = "Tip"\n'
        '[mainshock]\nt = "2023-02-06 01:17:32"\nlat = 37.0\nlon = 37.0\nmw = 7.8\n'
        "[analysis]\nwindow_days = 180\n",
    )
    assert main(["run", str(settings)]) == 0
    out = capsys.readouterr().out
    # three of the ten rows are deliberately unreadable in that fixture
    assert "10 read, 7 usable, 3 unreadable and left out" in out
    assert "converted to Mw by the Scordilis relations" in out


# --------------------------------------------------------------- fit verdict
# The rule that turns a p-value into a sentence had no test at all, which is
# how the browser page came to print "consistent with one Omori decay" from an
# uncalibrated p-value that the simulated null puts on the other side of 0.05.

from tremor_lab.cli import fit_verdict  # noqa: E402
from tremor_lab.omori import FitTest  # noqa: E402


def _test(p_value, se, method="parametric bootstrap, 200 replicates"):
    return FitTest(0.02, p_value, 1529, method, se)


def test_an_uncalibrated_p_value_gets_no_verdict_however_large_it_is():
    # 0.508 is the actual asymptotic value on the reference catalogue, where the
    # simulated null reads about 0.03. A verdict from it would be backwards.
    verdict = fit_verdict(_test(0.508, 0.0, "asymptotic, uncalibrated"))
    assert "uncalibrated" in verdict
    assert "consistent" not in verdict


def test_a_clear_rejection_is_stated_as_one():
    assert fit_verdict(_test(0.01, 0.007)) == "NOT consistent with a single Omori decay"


def test_a_clear_pass_is_stated_as_one():
    assert fit_verdict(_test(0.47, 0.035)) == "consistent with one Omori decay"


def test_a_p_value_within_monte_carlo_error_of_the_threshold_decides_nothing():
    # 0.06 looks like a pass and 0.04 like a failure, but at this many replicates
    # neither is separable from 0.05 and the difference is the seed, not the data.
    for p_value in (0.04, 0.05, 0.06):
        assert "borderline" in fit_verdict(_test(p_value, 0.02))


def test_a_missing_standard_error_falls_back_to_a_bare_comparison():
    nan = float("nan")
    assert fit_verdict(_test(0.60, nan)) == "consistent with one Omori decay"
    assert fit_verdict(_test(0.01, nan)) == "NOT consistent with a single Omori decay"


def test_the_browser_page_applies_the_same_rule_in_the_same_order():
    from pathlib import Path

    page = (Path(__file__).resolve().parents[1] / "web" / "tremor-lab.html").read_text(
        encoding="utf-8"
    )
    branch = page[
        page.index("var verdict;") : page.index('rows.push(["decay fit test"')
    ]
    order = [
        "uncalibrated",
        "pv + 2 * se < 0.05",
        "pv - 2 * se > 0.05",
        "borderline",
    ]
    positions = [branch.index(token) for token in order]
    assert positions == sorted(positions), "the page tests the branches out of order"


# ------------------------------------------- what a settings file may not outlive
# A run must not steer the next one. The first version of this guarantee covered
# only the refusals raised while reading the file, and every refusal raised while
# running it - a catalogue that is not there, a bad column mapping - happened after
# the constants had already been assigned.


def _settings_with(tmp_path, extra: str) -> Path:
    data = Path(__file__).parent / "data" / "kahramanmaras_180d.csv"
    text = f"""
[catalog]
path = "{data.as_posix()}"
[catalog.columns]
dt_days = "dt_days"
mag = "mw"
[mainshock]
t = "2023-02-06 01:17:32"
lat = 37.1578
lon = 36.8092
mw = 7.8
[analysis]
n_boot = 0
n_fit_simulations = 0
{extra}
"""
    path = tmp_path / "settings.toml"
    path.write_text(text, encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("extra", "expect_exit"),
    [
        ("[constants]\nMC_CORRECTION = 1.0", 0),  # accepted
        ("[constants]\nMC_CORRECTION = 1.0\nDELTA_MB = 2.0", 0),
    ],
)
def test_an_accepted_settings_file_puts_the_constants_back_when_it_is_done(
    tmp_path, capsys, extra, expect_exit
):
    from tremor_lab import constants

    before = (constants.MC_CORRECTION, constants.DELTA_MB)
    assert main(["run", str(_settings_with(tmp_path, extra))]) == expect_exit
    assert (constants.MC_CORRECTION, constants.DELTA_MB) == before


def test_a_file_refused_while_it_runs_puts_the_constants_back_too(tmp_path, capsys):
    from tremor_lab import constants

    before = (constants.MC_CORRECTION, constants.DELTA_MB)
    settings = tmp_path / "bad.toml"
    settings.write_text(
        """
[catalog]
path = "no_such_file.csv"
[catalog.columns]
dt_days = "dt_days"
mag = "mw"
[mainshock]
t = "2023-02-06 01:17:32"
lat = 37.1578
lon = 36.8092
mw = 7.8
[constants]
MC_CORRECTION = 1.0
DELTA_MB = 2.0
""",
        encoding="utf-8",
    )
    # Refused by _run, not by the reader: the catalogue path is only opened once
    # the file has been accepted, which is where the leak used to happen.
    assert main(["run", str(settings)]) == 2
    assert (constants.MC_CORRECTION, constants.DELTA_MB) == before


@pytest.mark.parametrize(
    "entry",
    ["mc_correction = true", "dm = true", "window_days = true", 'dm = "0.1"'],
)
def test_an_analysis_setting_that_is_not_a_number_is_refused_by_name(
    tmp_path, capsys, entry
):
    # These are the same quantities as the [constants] overrides by another route,
    # and they had no check at all: mc_correction = true ran to exit 0 and printed
    # a completeness magnitude of 4.2 instead of 3.4, because True counts as 1.
    assert main(["run", str(_settings_with(tmp_path, f"{entry}"))]) == 2
    message = capsys.readouterr().err
    assert "[analysis]" in message
    assert entry.split(" =")[0] in message
    assert "casting rule" not in message
