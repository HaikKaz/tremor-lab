"""The user guide must run, and its tables must say what the package does.

The README has the same guard in `test_documentation.py`, for the same reason: prose
about software drifts away from the software, and a reader who finds one wrong example
stops trusting the rest. Here the guide's Python is executed, and the tables a reader is
most likely to copy from, the published constants and the `[analysis]` keys, are checked
against the module and the command line that define them.
"""

import re
from pathlib import Path

import pytest

from tremor_lab import __version__, constants
from tremor_lab.cli import _ANALYSIS_NUMBERS

ROOT = Path(__file__).resolve().parents[1]
GUIDE = (ROOT / "docs" / "USER_GUIDE.md").read_text(encoding="utf-8")
README = (ROOT / "README.md").read_text(encoding="utf-8")

# Replicate counts are scaled down when the examples run, so the suite stays fast. The
# guide shows what a reader should type; this test checks that it executes.
_COUNTS = re.compile(r"\b(n_simulations|n_fit_simulations|n_boot)=\d+")


def _python_blocks() -> list[str]:
    return re.findall(r"^```python\n(.*?)^```$", GUIDE, re.S | re.M)


def test_every_python_example_in_the_guide_runs(monkeypatch):
    blocks = _python_blocks()
    assert len(blocks) >= 7, "the guide has lost its Python examples"
    # The examples name files by their path from the repository root, and one of them
    # reassigns a constant for the session. Neither may outlast this test.
    monkeypatch.chdir(ROOT)
    before = {name: getattr(constants, name) for name in constants._PUBLISHED}
    monkeypatch.setattr(constants, "N_BOOT", 20)
    monkeypatch.setattr(constants, "N_FIT_SIMULATIONS", 20)
    namespace: dict = {}
    try:
        for number, block in enumerate(blocks, start=1):
            code = _COUNTS.sub(lambda m: f"{m.group(1)}=20", block)
            source = f"USER_GUIDE.md, python block {number}"
            exec(compile(code, source, "exec"), namespace)
    finally:
        for name, value in before.items():
            setattr(constants, name, value)


def test_the_constants_table_lists_every_published_constant_with_its_default():
    rows = re.findall(r"^\| `([A-Z][A-Z0-9_]*)` \| ([0-9.e+-]+) \|", GUIDE, re.M)
    listed = {name: float(value) for name, value in rows}
    assert len(listed) == len(rows), "a constant is listed twice in the guide"
    assert set(listed) == set(constants._PUBLISHED), (
        "the guide's table and tremor_lab.constants name different constants"
    )
    for name, value in listed.items():
        assert value == pytest.approx(constants._PUBLISHED[name]), name


def test_the_analysis_table_names_exactly_the_keys_the_command_line_accepts():
    section = GUIDE.split("### The [analysis] section")[1].split("\n### ")[0]
    named = set(re.findall(r"^\| `([a-z_]+)` \|", section, re.M))
    assert named == _ANALYSIS_NUMBERS


def test_the_guide_states_the_version_the_command_line_prints():
    assert f"`Tremor Lab {__version__}`" in GUIDE
    assert f"Tremor Lab {__version__}\ncatalogue" in GUIDE


def test_the_readme_points_to_the_guide():
    assert "docs/USER_GUIDE.md" in README
