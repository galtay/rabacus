"""
Store the Cloudy results for the photoionized primordial models used by
tests/test_cloudy_cooling.py.  Run write_cooling_table.py first, then
run this script from this directory.
"""

import os

import models


os.makedirs(models.REFERENCE_DIR, exist_ok=True)
for name in ["primordial_pce", "primordial_pce_18"]:
    fname = os.path.join(models.REFERENCE_DIR, name + ".txt")
    models.write_reference(name, fname)
    print("wrote", fname)
