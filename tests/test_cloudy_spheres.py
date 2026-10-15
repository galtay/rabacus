"""
Compare rabacus Stromgren spheres against stored Cloudy 13.05 results.

The models and reference profiles live in ``cloudy/`` (see
``cloudy/spheres.py`` and ``cloudy/make_reference.sh``).  The two codes
are not expected to agree exactly: Cloudy uses more detailed atomic
physics and different rate fits, and only converges the ionization
fractions of each element to sum to unity within ~1e-4.  The tolerances
below are roughly 1.5 times the differences seen with Cloudy 13.05, so
these tests catch regressions in rabacus rather than measure accuracy.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "cloudy"))
import spheres  # noqa: E402


# Neutral fractions inside this radius are below Cloudy's convergence
# floor (xHI ~ 1e-7 to 1e-5), so profiles are only compared outside it.
R_MIN_KPC = 0.25

# (front radius relative tolerance, max absolute xHI difference)
TOLERANCES = {
    "iliev_test_1": (0.01, 0.05),
    "iliev_test_2": (0.01, 0.05),
    "raicevic_case_A": (0.01, 0.05),
    "raicevic_case_B": (0.01, 0.05),
    # rabacus' outward-only approximation and Cloudy's "diffuse fields
    # outward" treat recombination photons differently
    "raicevic_outward": (0.05, 0.25),
}
T_RTOL = 0.06  # iliev_test_2 temperature, outside R_MIN_KPC


@pytest.fixture(scope="module")
def sources():
    return spheres.make_sources()


def solve(name, sources):
    ref = spheres.read_reference(name)
    sphere = spheres.solve(name, ref, sources)
    r_cl = ref["r"].rescale("kpc").magnitude
    r_ra = sphere.r_c.rescale("kpc").magnitude
    return ref, sphere, r_cl, r_ra


@pytest.mark.parametrize("name", list(spheres.MODELS))
def test_ionization_front(name, sources):
    ref, sphere, r_cl, r_ra = solve(name, sources)
    rf_cl = spheres.front_radius(r_cl, ref["xHI"])
    rf_ra = spheres.front_radius(r_ra, sphere.xH1)
    rtol, _ = TOLERANCES[name]
    assert rf_ra == pytest.approx(rf_cl, rel=rtol)


@pytest.mark.parametrize("name", list(spheres.MODELS))
def test_neutral_fraction_profile(name, sources):
    ref, sphere, r_cl, r_ra = solve(name, sources)
    x_ra = np.interp(r_cl, r_ra, np.asarray(sphere.xH1))
    mask = r_cl > R_MIN_KPC
    _, atol = TOLERANCES[name]
    np.testing.assert_allclose(x_ra[mask], ref["xHI"][mask], rtol=0, atol=atol)


def test_iliev_test_2_temperature(sources):
    ref, sphere, r_cl, r_ra = solve("iliev_test_2", sources)
    T_ra = np.interp(r_cl, r_ra, sphere.T.magnitude)
    mask = r_cl > R_MIN_KPC
    np.testing.assert_allclose(
        T_ra[mask], ref["T"].magnitude[mask], rtol=T_RTOL)
