"""
Compare rabacus photo-collisional equilibrium against stored Cloudy 13.05
results for primordial gas exposed to the HM12 background at z = 3.

The models and reference results live in ``cloudy/cooling/`` (see
``models.py`` and ``make_reference.py``).  Only temperatures where
photoionization dominates are compared.  At higher temperatures the
codes differ by tens of percent because of their collisional rates and
line cooling (Cloudy uses multi-level atoms), and Cloudy includes
heating terms that rabacus does not, so these tests only check
quantities that depend mostly on photoionization and recombination.
The tolerances are roughly 1.5 times the differences seen with Cloudy
13.05.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "cloudy", "cooling"))
import models  # noqa: E402


T_MAX = 10**4.3  # photoionization dominated below this temperature

# relative tolerances for ionization fractions and heating
TOLERANCES = {
    "primordial_pce": dict(
        H1=0.04, H2=0.01, He1=0.10, He2=0.03, He3=0.13, heat=0.045),
    # the attenuated model agrees less well in xHI (~12%); not investigated
    "primordial_pce_18": dict(
        H1=0.18, H2=0.03, He1=0.13, He2=0.03, He3=0.13, heat=0.07),
}


@pytest.fixture(scope="module", params=list(TOLERANCES))
def model(request):
    name = request.param
    T, x_cl, h_cl, _ = models.read_reference(name)
    mask = T < T_MAX
    x_ra, h_ra, _ = models.solve_rabacus(
        T[mask], models.photo_rates(*models.MODELS[name]))
    x_cl = {k: v[mask] for k, v in x_cl.items()}
    return name, x_ra, h_ra, x_cl, h_cl[mask]


@pytest.mark.parametrize("species", ["H1", "H2", "He1", "He2", "He3"])
def test_ionization_fractions(model, species):
    name, x_ra, _, x_cl, _ = model
    rtol = TOLERANCES[name][species]
    np.testing.assert_allclose(x_ra[species], x_cl[species], rtol=rtol)


def test_heating(model):
    name, _, h_ra, _, h_cl = model
    np.testing.assert_allclose(h_ra, h_cl, rtol=TOLERANCES[name]["heat"])
