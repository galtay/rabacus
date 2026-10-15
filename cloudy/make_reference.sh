#!/bin/sh
# Run Cloudy on the Stromgren sphere inputs and store the reference
# profiles used by plot_comparisons.py and tests/test_cloudy_spheres.py.
#
#   ./make_reference.sh
#
# Uses $CLOUDY_EXE, or the build made by install_cloudy.sh.
set -eu

here=$(cd "$(dirname "$0")" && pwd)
cd "$here"

models="iliev_test_1 iliev_test_2 raicevic_case_A raicevic_case_B raicevic_outward"

./run $models
mkdir -p reference
for name in $models; do
    uv run python -c "import spheres; spheres.write_reference('$name.ovr', 'reference/$name.txt')"
done
ls -l reference
