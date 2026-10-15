"""
Stromgren sphere models compared between rabacus and Cloudy.

Each entry in ``MODELS`` matches a Cloudy input file ``<name>.in`` in
this directory.  The reference Cloudy results (last iteration of the
overview file, temperature and H/He ionization columns only) are kept in
``reference/<name>.txt`` and are produced by ``make_reference.sh``.
"""

import os

import numpy as np
import rabacus as ra


HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_DIR = os.path.join(HERE, "reference")

# columns kept from the Cloudy overview file
REFERENCE_COLUMNS = ["depth", "Te", "HI", "HII", "HeI", "HeII", "HeIII"]

R_SPHERE = 6.6 * ra.u.kpc
R_INNER = 1.0 * ra.u.pc  # Cloudy "radius 1 parsec"
NH = 1.0e-3 / ra.u.cm**3
NHE_NULL = 1.0e-15 / ra.u.cm**3
T_FIX = 1.0e4 * ra.u.K
LN = 5.0e48 / ra.u.s  # photon luminosity, Q(H) = 10**48.69897


def make_sources():
    """Return the monochromatic (1.1 Ry) and 1e5 K thermal point sources."""
    q_mono = 1.1
    mono = ra.rad_src.PointSource(q_mono, q_mono, "monochromatic")
    mono.normalize_Ln(LN)

    thrm = ra.rad_src.PointSource(
        1.0, 10.0, "thermal", T_eff=1.0e5 * ra.u.K, Nnu=128)
    thrm.normalize_Ln(LN)
    return {"mono": mono, "thermal": thrm}


# name -> (source, grid, SphereStromgren keyword arguments)
#   grid "cloudy": layer edges halfway between Cloudy zone centers
#   grid "uniform": as many evenly spaced layers as Cloudy zones
MODELS = {
    "iliev_test_1": ("mono", "cloudy", dict(fixed_fcA=0.0)),
    "iliev_test_2": ("thermal", "cloudy",
                     dict(fixed_fcA=0.0, find_Teq=True, z=0.0)),
    "raicevic_case_A": ("mono", "cloudy", dict(fixed_fcA=1.0)),
    "raicevic_case_B": ("mono", "cloudy", dict(fixed_fcA=0.0)),
    "raicevic_outward": ("mono", "uniform", dict(rec_meth="outward")),
}


def read_ovr(fname):
    """Return (header, data) for the last iteration of a Cloudy .ovr file."""
    header = None
    iterations = []
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                if header is None:
                    header = line[1:].rstrip("\n").split("\t")
                iterations.append([])
            elif line.strip():
                iterations[-1].append([float(v) for v in line.split("\t")])
    return header, np.array(iterations[-1])


def write_reference(ovr_fname, ref_fname):
    """Write the reference columns of the last iteration of a .ovr file."""
    header, dat = read_ovr(ovr_fname)
    icols = [header.index(c) for c in REFERENCE_COLUMNS]
    np.savetxt(
        ref_fname, dat[:, icols],
        fmt="%.5e", delimiter="\t", header="\t".join(REFERENCE_COLUMNS))


def read_reference(name, ref_dir=REFERENCE_DIR):
    """Return a dict of Cloudy reference profiles with units attached.

    Temperatures and ionization fractions are returned linearly (Cloudy
    writes log10 values).
    """
    dat = np.loadtxt(os.path.join(ref_dir, name + ".txt"))
    cols = dict(zip(REFERENCE_COLUMNS, dat.T))
    ref = {"r": (cols["depth"] * ra.u.cm + R_INNER).rescale("kpc"),
           "T": 10**cols["Te"] * ra.u.K}
    for k in REFERENCE_COLUMNS[2:]:
        ref["x" + k] = 10**cols[k]
    return ref


def layer_edges(ref, grid):
    """Layer edges for a rabacus sphere matched to a Cloudy profile."""
    r_c = ref["r"].rescale("kpc").magnitude
    if grid == "cloudy":
        edges = np.zeros(r_c.size + 1)
        edges[1:-1] = 0.5 * (r_c[:-1] + r_c[1:])
        edges[-1] = R_SPHERE.rescale("kpc").magnitude
        return edges * ra.u.kpc
    if grid == "uniform":
        return np.linspace(0.0 * ra.u.kpc, R_SPHERE.rescale("kpc"),
                           r_c.size + 1)
    raise ValueError(grid)


def solve(name, ref, sources=None):
    """Solve the rabacus sphere for model `name` on the grid of `ref`."""
    if sources is None:
        sources = make_sources()
    src, grid, kwargs = MODELS[name]
    edges = layer_edges(ref, grid)
    Nl = edges.size - 1
    T = np.ones(Nl) * T_FIX
    nH = np.ones(Nl) * NH
    nHe = np.ones(Nl) * NHE_NULL
    return ra.f2py.SphereStromgren(edges, T, nH, nHe, sources[src], **kwargs)


def front_radius(r, xH1):
    """Radius where the neutral fraction first exceeds 0.5."""
    r = np.asarray(r)
    xH1 = np.asarray(xH1)
    i = np.argmax(xH1 > 0.5)
    if i == 0:
        return np.nan
    return np.interp(0.5, xH1[i - 1:i + 1], r[i - 1:i + 1])
