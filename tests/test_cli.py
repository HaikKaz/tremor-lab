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
