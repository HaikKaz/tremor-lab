"""The prose and the browser page must not drift away from the package.

Two kinds of drift have actually happened here, so both are now checked by a test
rather than by remembering.

*Numbers quoted in the README.* A constant was added and the README went on saying
there were twenty-three of them. Nothing broke, but a reader counting the list would
have found the documentation wrong, and a reader who trusts the documentation about
one thing trusts it about others.

*Constants duplicated in the browser page.* The page is a second implementation in a
different language; the same published values appear in its own `CONST` block. Nothing
makes the two agree except that someone checks, and once they did not: the page
calibrated the decay-fit test with 120 replicates where the package used 200.
"""

import re
from pathlib import Path

import pytest

from tremor_lab import constants

ROOT = Path(__file__).resolve().parents[1]
README = (ROOT / "README.md").read_text(encoding="utf-8")
PAGE = (ROOT / "web" / "tremor-lab.html").read_text(encoding="utf-8")

PUBLISHED = constants._PUBLISHED


def test_the_readme_states_the_true_number_of_published_constants():
    claimed = re.findall(r"All (\d+) live in\n`tremor_lab\.constants`", README)
    assert claimed, "the README no longer states how many constants there are"
    assert int(claimed[0]) == len(PUBLISHED)


def test_the_readme_arithmetic_about_the_three_exceptions_adds_up():
    # "Three of the N ... The other M can be set all three ways."
    total = int(re.search(r"Three of the (\d+) are the exception", README).group(1))
    others = int(
        re.search(r"The other (\d+) can be set all three ways", README).group(1)
    )
    assert total == len(PUBLISHED)
    assert others == total - 3


def test_the_three_constants_the_readme_calls_exceptions_really_have_no_keyword():
    # They judge whether a fitted decay is worth believing; they do not enter the
    # fit, so there is nothing to pass them to.
    import inspect

    from tremor_lab import fit_omori

    signature = inspect.signature(fit_omori).parameters
    for name in ("OMORI_C_FLOOR", "OMORI_P_MIN", "OMORI_P_MAX"):
        assert name in PUBLISHED
        assert name.lower() not in signature


# Every value the browser page hard-codes, and the published constant it must equal.
# A name here that the page stops using, or a value that moves on either side, fails.
PAGE_CONSTANTS = {
    "DELTA_MB": "DELTA_MB",
    "E_A": "ENERGY_A",
    "E_B": "ENERGY_B",
    "DM": "DM",
    "GRID_TOLERANCE": "GRID_TOLERANCE",
    "N_BOOT": "N_BOOT",
    "WINDOW_DAYS": "WINDOW_DAYS",
    "OMORI_C0": "OMORI_C0",
    "OMORI_P0": "OMORI_P0",
    "MC_CORR": "MC_CORRECTION",
    "MIN_FOR_MC": "MIN_EVENTS_FOR_MC",
    "MIN_FOR_B": "MIN_EVENTS_FOR_FIT",
    "N_FIT_SIMULATIONS": "N_FIT_SIMULATIONS",
    "OMORI_C_FLOOR": "OMORI_C_FLOOR",
    "OMORI_P_MIN": "OMORI_P_MIN",
    "OMORI_P_MAX": "OMORI_P_MAX",
}


def _page_const_block() -> dict[str, float]:
    block = re.search(r"var CONST = \{(.*?)\n\};", PAGE, re.S)
    assert block, "the page's CONST block has been renamed or removed"
    found = {}
    for name, value in re.findall(
        r"^\s*([A-Z][A-Z0-9_]*)\s*:\s*([0-9eE.+-]+)\s*,?\s*(?://.*)?$",
        block.group(1),
        re.M,
    ):
        found[name] = float(value)
    return found


@pytest.mark.parametrize(("in_page", "in_package"), sorted(PAGE_CONSTANTS.items()))
def test_the_browser_page_carries_the_same_value_as_the_package(in_page, in_package):
    page = _page_const_block()
    assert in_page in page, f"{in_page} has gone from the page's CONST block"
    assert page[in_page] == pytest.approx(PUBLISHED[in_package])


def test_no_constant_in_the_page_is_left_unchecked():
    unknown = set(_page_const_block()) - set(PAGE_CONSTANTS)
    assert not unknown, f"add these to PAGE_CONSTANTS so they cannot drift: {unknown}"


# Values the page writes as literals inside a formula rather than in CONST. They are
# just as capable of drifting, and a wrong Scordilis coefficient would be invisible.
PAGE_LITERALS = [
    (r"ms <= ([\d.]+)\)", "MS_MW_BRANCH"),
    (r"\? ([\d.]+) \* ms \+", "MS_MW_LOW_SLOPE"),
    (r"\* ms \+ ([\d.]+) :", "MS_MW_LOW_INTERCEPT"),
    (r": ([\d.]+) \* ms \+", "MS_MW_HIGH_SLOPE"),
    (r"\* ms \+ ([\d.]+); \}", "MS_MW_HIGH_INTERCEPT"),
    (r"return ([\d.]+) \* mb \+", "MB_MW_SLOPE"),
    (r"\* mb \+ ([\d.]+); \}", "MB_MW_INTERCEPT"),
    (r"var R = (\d+), r = Math\.PI", "EARTH_RADIUS_KM"),
    (r"var sigma = ([\d.]+) \* b \* b", "SHI_BOLT_K"),
]


@pytest.mark.parametrize(("pattern", "name"), PAGE_LITERALS)
def test_the_page_formulas_use_the_published_coefficients(pattern, name):
    match = re.search(pattern, PAGE)
    assert match, f"the page formula carrying {name} has been rewritten"
    assert float(match.group(1)) == pytest.approx(PUBLISHED[name])


def test_the_readme_names_exactly_the_settings_the_page_lets_a_reader_change():
    # The README once claimed every constant was editable in the browser. Seven are.
    claimed = re.search(
        r"(\w+) settings are editable on the page - (.+?) - and the", README, re.S
    )
    assert claimed, "the README no longer says which settings the page exposes"
    named = [
        part.strip() for part in re.split(r",| and ", claimed.group(2)) if part.strip()
    ]
    ids = [
        "windowDays",
        "threshold",
        "dm",
        "mcCorr",
        "deltaMb",
        "nBoot",
        "nFitSim",
        "seed",
    ]
    assert claimed.group(1) == "Eight"
    assert len(named) == len(ids)
    for element_id in ids:
        assert f'id="{element_id}"' in PAGE, f"{element_id} is no longer on the page"
