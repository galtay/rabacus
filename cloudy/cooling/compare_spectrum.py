"""
Compare the HM12 spectrum Cloudy actually used with the rabacus spectrum.

Run write_cooling_table.py first, then run this script from this
directory.  It reads the incident continuum of the unattenuated
primordial_pce model, plots it against the rabacus spectrum that was
written to spectrum.ini, and prints the H I, He I and He II
photoionization rates obtained by integrating Cloudy's continuum and by
rabacus with increasing numbers of energy samples.
"""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pylab as plt
import rabacus as ra


z = 3.0
q_min = 1.0e-2
q_max = 1.0e6
NNU_CLOUDY = 100  # samples written to spectrum.ini by write_cooling_table.py
SPECIES = ["H1", "He1", "He2"]


def read_incident(fname):
    """Return energy [Ry] and incident 4 pi nu J_nu [erg/s/cm^2] for the
    first grid point of a Cloudy "save continuum" file."""
    rows = []
    with open(fname) as f:
        f.readline()
        for line in f:
            if line.startswith("#"):
                break
            fields = line.split("\t")
            rows.append((float(fields[0]), float(fields[1])))
    dat = np.array(rows)
    return dat[:, 0], dat[:, 1]


def cloudy_rates(E_Ry, nu4piJ, px):
    """Photoionization rates [1/s] from Cloudy's incident continuum."""
    h = ra.pc.h.rescale("erg*s").magnitude
    nu = (E_Ry * ra.u.Ry_inf / ra.pc.h).rescale("Hz").magnitude
    rates = {}
    for k in SPECIES:
        Eth = getattr(px, "Eth_" + k).rescale("Ry_inf").magnitude
        m = E_Ry >= Eth
        sigma = getattr(px, "sigma_" + k)(E_Ry[m] * ra.u.Ry_inf)
        dG_dlnnu = nu4piJ[m] / (h * nu[m]) * sigma.rescale("cm**2").magnitude
        rates[k] = np.trapezoid(dG_dlnnu, np.log(nu[m]))
    return rates


def rabacus_spectrum(Nnu):
    uvb = ra.BackgroundSource(q_min, q_max, "hm12", z=z, Nnu=Nnu)
    E_Ry = uvb.E.rescale("Ry_inf").magnitude
    nu = (uvb.E / ra.pc.h).rescale("Hz").magnitude
    Inu = uvb.Inu.rescale("erg/(s*cm**2*Hz*sr)").magnitude
    return uvb, E_Ry, 4.0 * np.pi * nu * Inu


E_cl, nu4piJ_cl = read_incident("primordial_pce.con")
uvb, E_ra, nu4piJ_ra = rabacus_spectrum(NNU_CLOUDY)

# spectrum shape at the rabacus sample energies
ratio = np.exp(np.interp(np.log(E_ra), np.log(E_cl),
                         np.log(np.maximum(nu4piJ_cl, 1.0e-300))))
ratio = ratio / nu4piJ_ra
# Cloudy's energy mesh ends near 1e6 Ry, so skip the last samples
ionizing = (E_ra >= 1.0) & (E_ra <= 1.0e5)
print(f"cloudy / rabacus spectrum from 1 to 1e5 Ry: "
      f"min {ratio[ionizing].min():.4f}, max {ratio[ionizing].max():.4f}")

# photoionization rates
G_cl = cloudy_rates(E_cl, nu4piJ_cl, uvb.PX)
print("photoionization rates, rabacus / cloudy continuum:")
for Nnu in [NNU_CLOUDY, 1000, 5000]:
    u = ra.BackgroundSource(q_min, q_max, "hm12", z=z, Nnu=Nnu)
    txt = "  ".join(
        f"{k}i {getattr(u.thin, k + 'i').rescale('1/s').magnitude / G_cl[k]:.3f}"
        for k in SPECIES)
    print(f"  Nnu = {Nnu:5d}: {txt}")

fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True,
                         gridspec_kw=dict(height_ratios=[3, 1]))
ax = axes[0]
ax.loglog(E_cl, nu4piJ_cl, color="black", lw=1.0, label="cloudy incident")
ax.loglog(E_ra, nu4piJ_ra, "o", color="red", ms=3,
          label=f"rabacus (Nnu={NNU_CLOUDY})")
ax.set_xlim(q_min, q_max)
ax.set_ylim(1.0e-8, 1.0e-2)
ax.set_ylabel(r"$4 \pi \nu J_\nu$ [erg s$^{-1}$ cm$^{-2}$]")
ax.legend(loc="best")
ax = axes[1]
ax.semilogx(E_ra, ratio, "o", color="red", ms=3)
ax.axhline(1.0, color="black", lw=0.5)
ax.set_ylim(0.9, 1.1)
ax.set_xlabel("E [Ry]")
ax.set_ylabel("cloudy / rabacus")
fig.tight_layout()
fig.savefig("compare_spectrum.png")
