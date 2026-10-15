"""
Overlay rabacus and Cloudy ionization fractions, heating and cooling for
primordial (H + He) gas in a single zone at fixed temperature.

Run write_cooling_table.py first to produce the Cloudy grids, then run
this script from this directory.  It reproduces each Cloudy model with
rabacus (same HM12 spectrum at z = 3, same attenuating columns, nH =
1e-2 cm^-3, case A rates), saves compare_cooling.png and prints the
median relative difference of each quantity.
"""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pylab as plt
import rabacus as ra

from read_cloudy import (read_grid_file, read_heating_file,
                         read_overview_file)


z = 3.0
nH_val = 1.0e-2
y_He = 10**(-1.0701)  # He/H by number in Cloudy's GASS10 abundances

# Cloudy model name -> (spectrum scale factor, NH1, NHe1, NHe2 [cm^-2])
# as set up in write_cooling_table.py
MODELS = {
    "primordial_ce": (1.0e-20, 0.0, 0.0, 0.0),
    "primordial_pce": (1.0, 0.0, 0.0, 0.0),
    "primordial_pce_18": (1.0, 1.0e18, 1.0e17, 1.0e17),
}

SPECIES = [("H1", "HI", "red", "-"), ("H2", "HII", "red", "--"),
           ("He1", "HeI", "blue", "-"), ("He2", "HeII", "blue", "--"),
           ("He3", "HeIII", "blue", ":")]


def read_cloudy_model(fbase):
    NT = read_grid_file(fbase)
    T, xH1, xH2, xHe1, xHe2, xHe3 = read_overview_file(fbase, NT)
    _, htot, ctot, _, _, _ = read_heating_file(fbase, NT)
    x = dict(H1=xH1, H2=xH2, He1=xHe1, He2=xHe2, He3=xHe3)
    return T, x, htot, ctot


def photo_rates(fac, NH1, NHe1, NHe2):
    """HM12 photoionization and heating rates behind the given columns."""
    uvb = ra.BackgroundSource(1.0e-2, 1.0e6, "hm12", z=z, Nnu=100)
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
    xs = {k: np.asarray(getattr(x, k)) for k, _, _, _ in SPECIES}
    return xs, heat.magnitude, cool.magnitude


def median_rel(a, b):
    return np.median(np.abs(a - b) / np.abs(b))


fig, axes = plt.subplots(len(MODELS), 3, figsize=(13, 3.6 * len(MODELS)),
                         squeeze=False)

for row, (fbase, (fac, NH1, NHe1, NHe2)) in zip(axes, MODELS.items()):
    T, x_cl, h_cl, c_cl = read_cloudy_model(fbase)
    x_ra, h_ra, c_ra = solve_rabacus(T, photo_rates(fac, NH1, NHe1, NHe2))

    print(f"{fbase}: median |rabacus - cloudy| / cloudy")
    for k, label, _, _ in SPECIES:
        print(f"  x{label:6s} {median_rel(x_ra[k], x_cl[k]):.3f}")
    if fac == 1.0:
        print(f"  heating {median_rel(h_ra, h_cl):.3f}")
    print(f"  cooling {median_rel(c_ra, c_cl):.3f}")

    ax = row[0]
    for k, label, color, ls in SPECIES:
        ax.loglog(T, x_ra[k], color=color, ls=ls, lw=2.5, label=label)
        ax.loglog(T, x_cl[k], color="black", ls=ls, lw=1.0)
    ax.set_ylim(1.0e-6, 2.0)
    ax.set_ylabel("x")
    ax.set_title(f"{fbase}: ionization fractions", fontsize=10)
    ax.legend(fontsize=7, ncol=2, loc="lower left")

    ax = row[1]
    if fac == 1.0:
        ax.loglog(T, h_ra / nH_val**2, color="red", lw=2.5, label="rabacus")
        ax.loglog(T, h_cl / nH_val**2, color="black", lw=1.0, label="cloudy")
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "no radiation field", ha="center",
                transform=ax.transAxes)
    ax.set_ylabel(r"$\mathcal{H} / n_{\rm H}^2$")
    ax.set_title(f"{fbase}: heating", fontsize=10)

    ax = row[2]
    ax.loglog(T, c_ra / nH_val**2, color="blue", lw=2.5, label="rabacus")
    ax.loglog(T, c_cl / nH_val**2, color="black", lw=1.0, label="cloudy")
    ax.set_ylabel(r"$\Lambda / n_{\rm H}^2$")
    ax.set_title(f"{fbase}: cooling", fontsize=10)
    ax.legend(fontsize=7)

for ax in axes[-1]:
    ax.set_xlabel("T [K]")

fig.tight_layout()
fig.savefig("compare_cooling.png")
