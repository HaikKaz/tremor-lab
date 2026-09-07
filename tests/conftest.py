"""Shared test setup.

The decay fit test simulates its own null distribution, and the published default
is 600 replicates - about three seconds per call on the reference sequence. Most
tests here call `analyze_case` for something else entirely and would pay that
without ever looking at the p-value, so the default is lowered for the session.

Any test that is actually about the calibration passes its own count, which
overrides this, and `test_the_published_default_is_paid_in_full` runs one call at
the real 600 so the shipped default is exercised rather than merely declared.
"""

import pytest

from tremor_lab import constants

CHEAP_SIMULATIONS = 24


@pytest.fixture(autouse=True)
def _cheap_fit_calibration(request):
    """Lower the replicate count unless a test asks for the published one."""
    if request.node.get_closest_marker("full_calibration"):
        yield
        return
    published = constants.N_FIT_SIMULATIONS
    constants.N_FIT_SIMULATIONS = CHEAP_SIMULATIONS
    try:
        yield
    finally:
        constants.N_FIT_SIMULATIONS = published


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "full_calibration: run with the published replicate count, not the cheap one",
    )
