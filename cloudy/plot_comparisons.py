"""
Compare rabacus Stromgren spheres against Cloudy.

Run Cloudy on the *.in files in this directory first (see ``run``), then
run this script from this directory.  For each test it solves a rabacus
sphere on the Cloudy zone grid, saves a plot of the neutral fraction
(and temperature for iliev_test_2) to ``png/``, and prints the radius at
which the neutral fraction crosses 0.5 in both codes.
"""

import os
import time

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pylab as plt
import rabacus as ra


Rsphere = 6.6 * ra.u.kpc
R_inner = 1.0 * ra.u.pc  # Cloudy "radius 1 parsec"
nH_val = 1.0e-3
nHe_null_val = 1.0e-15

Ln = 5.0e48 / ra.u.s  # photon luminosity, Q(H) = 10**48.69897
Nnu = 128
png_dir = "png"


# sources
#=================================================================

q_mono = 1.1
src_mono = ra.rad_src.PointSource(q_mono, q_mono, "monochromatic")
src_mono.normalize_Ln(Ln)

T_eff = 1.0e5 * ra.u.K
src_thrm = ra.rad_src.PointSource(1.0, 10.0, "thermal", T_eff=T_eff, Nnu=Nnu)
src_thrm.normalize_Ln(Ln)


# helpers
#=================================================================

def read_cloudy(fbase):
    """Return the last iteration of a Cloudy overview (.ovr) file.

    Column 0 is depth [cm], 1 is log T, 6 is log xHI.
    """
    with open(fbase + ".ovr") as f:
        lines = f.readlines()

    # each iteration starts with a commented header line
    iterations = []
    for line in lines:
        if line.startswith("#"):
            iterations.append([])
        elif line.strip() and not line.startswith("#"):
            iterations[-1].append([float(v) for v in line.split("\t")])
    return np.array(iterations[-1])


def cloudy_radius(dat):
    return dat[:, 0] * ra.u.cm + R_inner.rescale("cm")


def edges_from_cloudy(dat):
    """Layer edges halfway between Cloudy zone centers."""
    r_c = cloudy_radius(dat).magnitude
    edges = np.zeros(r_c.size + 1)
    edges[1:-1] = 0.5 * (r_c[:-1] + r_c[1:])
    edges[-1] = Rsphere.rescale("cm").magnitude
    return edges * ra.u.cm


def sphere_inputs(Nl):
    T = np.ones(Nl) * 1.0e4 * ra.u.K
    nH = np.ones(Nl) * nH_val / ra.u.cm**3
    nHe_null = np.ones(Nl) * nHe_null_val / ra.u.cm**3
    return T, nH, nHe_null


def front_radius(r, xH1):
    """Radius [kpc] where the neutral fraction first exceeds 0.5."""
    r = np.asarray(r)
    xH1 = np.asarray(xH1)
    i = np.argmax(xH1 > 0.5)
    if i == 0:
        return np.nan
    return np.interp(0.5, xH1[i - 1:i + 1], r[i - 1:i + 1])


def solve(key, Edges, src, **kwargs):
    Nl = Edges.size - 1
    T, nH, nHe_null = sphere_inputs(Nl)
    t1 = time.time()
    s = ra.f2py.SphereStromgren(Edges, T, nH, nHe_null, src, **kwargs)
    print(f"{key}: solved {Nl} layers in {time.time() - t1:.2f} s")
    return s


def report(key, s, dat):
    r_ra = s.r_c.rescale("kpc").magnitude
    r_cl = cloudy_radius(dat).rescale("kpc").magnitude
    rf_ra = front_radius(r_ra, s.xH1)
    rf_cl = front_radius(r_cl, 10**dat[:, 6])
    print(
        f"{key}: xHI=0.5 at r = {rf_ra:.4f} kpc (rabacus), "
        f"{rf_cl:.4f} kpc (cloudy), rel diff = {(rf_ra - rf_cl) / rf_cl:+.2%}"
    )


def plot_xH1(ax, s, dat, ls="-", label=None):
    ax.plot(s.r_c.rescale("kpc"), np.log10(s.xH1), color="red", ls=ls,
            label=None if label is None else "rabacus " + label)
    ax.plot(cloudy_radius(dat).rescale("kpc"), dat[:, 6], color="blue",
            ls=ls, label=None if label is None else "cloudy " + label)


os.makedirs(png_dir, exist_ok=True)


#=================================================================
# Iliev 06 tests
#=================================================================

# test 1: case B, monochromatic, fixed T
#-----------------------------------------------------------------
key = "iliev_test_1"
dat = read_cloudy(key)
s = solve(key, edges_from_cloudy(dat), src_mono, fixed_fcA=0.0)
report(key, s, dat)

fig, ax = plt.subplots()
plot_xH1(ax, s, dat, label="")
ax.set_xlabel("r [kpc]")
ax.set_ylabel(r"$\log_{10} x_{\rm HI}$")
ax.set_title(key)
ax.legend(loc="best")
fig.savefig(f"{png_dir}/{key}_xH1.png")


# test 2: case B, thermal spectrum, equilibrium T
#-----------------------------------------------------------------
key = "iliev_test_2"
dat = read_cloudy(key)
s = solve(key, edges_from_cloudy(dat), src_thrm,
          fixed_fcA=0.0, find_Teq=True, z=0.0)
report(key, s, dat)

fig, ax = plt.subplots()
plot_xH1(ax, s, dat, label="")
ax.set_xlabel("r [kpc]")
ax.set_ylabel(r"$\log_{10} x_{\rm HI}$")
ax.set_title(key)
ax.legend(loc="best")
fig.savefig(f"{png_dir}/{key}_xH1.png")

fig, ax = plt.subplots()
ax.plot(s.r_c.rescale("kpc"), np.log10(s.T.magnitude), color="red",
        label="rabacus")
ax.plot(cloudy_radius(dat).rescale("kpc"), dat[:, 1], color="blue",
        label="cloudy")
ax.set_xlabel("r [kpc]")
ax.set_ylabel(r"$\log_{10} T$ [K]")
ax.set_title(key)
ax.legend(loc="best")
fig.savefig(f"{png_dir}/{key}_T.png")


#=================================================================
# Raicevic 14 tests
#=================================================================

spheres = {}
dats = {}

# case A, monochromatic, fixed T
key = "raicevic_case_A"
dats[key] = read_cloudy(key)
spheres[key] = solve(key, edges_from_cloudy(dats[key]), src_mono,
                     fixed_fcA=1.0)

# case B, monochromatic, fixed T
key = "raicevic_case_B"
dats[key] = read_cloudy(key)
spheres[key] = solve(key, edges_from_cloudy(dats[key]), src_mono,
                     fixed_fcA=0.0)

# outward-only diffuse field, monochromatic, fixed T
key = "raicevic_outward"
dats[key] = read_cloudy(key)
Nl = dats[key][:, 0].size
Edges = np.linspace(0.0 * ra.u.kpc, Rsphere.rescale("kpc"), Nl + 1)
spheres[key] = solve(key, Edges, src_mono, rec_meth="outward")

fig, ax = plt.subplots()
for key, ls in [("raicevic_case_A", "-"), ("raicevic_case_B", ":"),
                ("raicevic_outward", "--")]:
    report(key, spheres[key], dats[key])
    plot_xH1(ax, spheres[key], dats[key], ls=ls,
             label=key.replace("raicevic_", ""))
ax.set_xlim(1.6, 6.8)
ax.set_ylim(-2.5, 0.1)
ax.set_xlabel("r [kpc]")
ax.set_ylabel(r"$\log_{10} x_{\rm HI}$")
ax.set_title("raicevic 14")
ax.legend(loc="best", fontsize=8)
fig.savefig(f"{png_dir}/raicevic_xH1.png")
