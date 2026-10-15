"""
Single zone primordial gas models compared between rabacus and Cloudy.

Each entry in ``MODELS`` matches a Cloudy grid written by
write_cooling_table.py (fixed temperature, nH = 1e-2 cm^-3, HM12 at
z = 3).  Reference Cloudy results for the photoionized models are kept
in ``reference/<name>.txt`` and are produced by ``make_reference.py``.
"""

import os

import numpy as np
import rabacus as ra

from read_cloudy import (read_grid_file, read_heating_file,
                         read_overview_file)


HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_DIR = os.path.join(HERE, "reference")

z = 3.0
nH_val = 1.0e-2
y_He = 10**(-1.0701)  # He/H by number in Cloudy's GASS10 abundances

# Cloudy is given the spectrum at 100 energies (write_cooling_table.py) and
# interpolates it on a fine mesh.  rabacus integrates rates with the
# trapezoid rule on its own energy samples, which overestimates the H I
# photoionization rate by ~40% with 100 samples over this range, so the
# rates are computed from a finely sampled spectrum instead.
NNU_RATES = 5000

# Cloudy model name -> (spectrum scale factor, NH1, NHe1, NHe2 [cm^-2])
# as set up in write_cooling_table.py
MODELS = {
    "primordial_ce": (1.0e-20, 0.0, 0.0, 0.0),
    "primordial_pce": (1.0, 0.0, 0.0, 0.0),
    "primordial_pce_18": (1.0, 1.0e18, 1.0e17, 1.0e17),
}

SPECIES = ["H1", "H2", "He1", "He2", "He3"]
REFERENCE_COLUMNS = ["T"] + ["x" + k for k in SPECIES] + ["heat", "cool"]


def read_cloudy_model(fbase):
    """Return T [K], ionization fractions, heating and cooling
    [erg/s/cm^3] from the Cloudy output files of grid `fbase`."""
    NT = read_grid_file(fbase)
    T, xH1, xH2, xHe1, xHe2, xHe3 = read_overview_file(fbase, NT)
    _, htot, ctot, _, _, _ = read_heating_file(fbase, NT)
    x = dict(H1=xH1, H2=xH2, He1=xHe1, He2=xHe2, He3=xHe3)
    return T, x, htot, ctot


def write_reference(fbase, ref_fname):
    T, x, heat, cool = read_cloudy_model(fbase)
    dat = np.column_stack([T] + [x[k] for k in SPECIES] + [heat, cool])
    np.savetxt(ref_fname, dat, fmt="%.5e", delimiter="\t",
               header="\t".join(REFERENCE_COLUMNS))


def read_reference(name, ref_dir=REFERENCE_DIR):
    """Return (T, x, heat, cool) as stored by write_reference."""
    dat = np.loadtxt(os.path.join(ref_dir, name + ".txt"))
    cols = dict(zip(REFERENCE_COLUMNS, dat.T))
    x = {k: cols["x" + k] for k in SPECIES}
    return cols["T"], x, cols["heat"], cols["cool"]


def photo_rates(fac, NH1, NHe1, NHe2, Nnu=NNU_RATES):
    """HM12 photoionization and heating rates behind the given columns."""
    uvb = ra.BackgroundSource(1.0e-2, 1.0e6, "hm12", z=z, Nnu=Nnu)
    uvb.scale_spectrum(fac)
    if NH1 == NHe1 == NHe2 == 0.0:
        t = uvb.thin
        return dict(H1i=t.H1i, He1i=t.He1i, He2i=t.He2i,
                    H1h=t.H1h, He1h=t.He1h, He2h=t.He2h)
    cm2 = ra.u.cm**2
    tau = (float(NH1 * uvb.th.sigma_H1.rescale(cm2).magnitude),
           float(NHe1 * uvb.th.sigma_He1.rescale(cm2).magnitude),
           float(NHe2 * uvb.th.sigma_He2.rescale(cm2).magnitude))
    return dict(H1i=uvb.shld_H1i(*tau), He1i=uvb.shld_He1i(*tau),
                He2i=uvb.shld_He2i(*tau), H1h=uvb.shld_H1h(*tau),
                He1h=uvb.shld_He1h(*tau), He2h=uvb.shld_He2h(*tau))


def solve_rabacus(T, rates):
    """Return (x, heat, cool) from rabacus for gas at temperatures T [K]
    in photo-collisional equilibrium with case A recombination rates."""
    T = T * ra.u.K
    nH = np.ones(T.size) * nH_val / ra.u.cm**3
    nHe = nH * y_He
    ones = np.ones(T.size)
    kchem = ra.ChemistryRates(
        T, 1.0, 1.0, 1.0, H1i=ones * rates["H1i"],
        He1i=ones * rates["He1i"], He2i=ones * rates["He2i"])
    kcool = ra.CoolingRates(
        T, 1.0, 1.0, 1.0, H1h=ones * rates["H1h"],
        He1h=ones * rates["He1h"], He2h=ones * rates["He2h"])
    x = ra.Solve_PCE(nH, nHe, kchem)
    heat = kcool.return_heating(nH, nHe, x).rescale("erg/s/cm**3")
    cool = kcool.return_cooling(nH, nHe, x, z).rescale("erg/s/cm**3")
    xs = {k: np.asarray(getattr(x, k)) for k in SPECIES}
    return xs, heat.magnitude, cool.magnitude
