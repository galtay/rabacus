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

from models import (MODELS, nH_val, photo_rates, read_cloudy_model,
                    solve_rabacus)


# key, label, color, line style
SPECIES = [("H1", "HI", "red", "-"), ("H2", "HII", "red", "--"),
           ("He1", "HeI", "blue", "-"), ("He2", "HeII", "blue", "--"),
           ("He3", "HeIII", "blue", ":")]


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
