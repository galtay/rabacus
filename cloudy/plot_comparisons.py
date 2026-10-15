"""
Plot rabacus Stromgren spheres against the Cloudy reference profiles.

Reads ``reference/*.txt`` (see ``make_reference.sh``), solves the
matching rabacus spheres (see ``spheres.py``), saves plots to ``png/``
and prints the radius at which the neutral fraction crosses 0.5 in
both codes.  Run from any directory.
"""

import os

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pylab as plt

import spheres


png_dir = os.path.join(spheres.HERE, "png")
os.makedirs(png_dir, exist_ok=True)

sources = spheres.make_sources()
refs = {}
sols = {}
for name in spheres.MODELS:
    refs[name] = spheres.read_reference(name)
    sols[name] = spheres.solve(name, refs[name], sources)

    r_ra = sols[name].r_c.rescale("kpc").magnitude
    r_cl = refs[name]["r"].magnitude
    rf_ra = spheres.front_radius(r_ra, sols[name].xH1)
    rf_cl = spheres.front_radius(r_cl, refs[name]["xHI"])
    print(
        f"{name}: xHI=0.5 at r = {rf_ra:.4f} kpc (rabacus), "
        f"{rf_cl:.4f} kpc (cloudy), rel diff = {(rf_ra - rf_cl) / rf_cl:+.2%}"
    )


def plot_xH1(ax, name, ls="-", label=""):
    s = sols[name]
    ref = refs[name]
    ax.plot(s.r_c.rescale("kpc"), np.log10(s.xH1), color="red", ls=ls,
            label=f"rabacus {label}")
    ax.plot(ref["r"], np.log10(ref["xHI"]), color="blue", ls=ls,
            label=f"cloudy {label}")


def finish(fig, ax, title, fname, ylabel=r"$\log_{10} x_{\rm HI}$"):
    ax.set_xlabel("r [kpc]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)
    fig.savefig(os.path.join(png_dir, fname))


# Iliev 06 tests
#=================================================================

for name in ["iliev_test_1", "iliev_test_2"]:
    fig, ax = plt.subplots()
    plot_xH1(ax, name)
    finish(fig, ax, name, f"{name}_xH1.png")

name = "iliev_test_2"
fig, ax = plt.subplots()
ax.plot(sols[name].r_c.rescale("kpc"), np.log10(sols[name].T.magnitude),
        color="red", label="rabacus")
ax.plot(refs[name]["r"], np.log10(refs[name]["T"].magnitude),
        color="blue", label="cloudy")
finish(fig, ax, name, f"{name}_T.png", ylabel=r"$\log_{10} T$ [K]")


# Raicevic 14 tests
#=================================================================

fig, ax = plt.subplots()
for name, ls in [("raicevic_case_A", "-"), ("raicevic_case_B", ":"),
                 ("raicevic_outward", "--")]:
    plot_xH1(ax, name, ls=ls, label=name.replace("raicevic_", ""))
ax.set_xlim(1.6, 6.8)
ax.set_ylim(-2.5, 0.1)
finish(fig, ax, "raicevic 14", "raicevic_xH1.png")
