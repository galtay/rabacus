#!/bin/sh
# Download and build Cloudy 13.05 (https://nublado.org), the version used
# for the comparisons in this directory.
#
#   ./install_cloudy.sh [install_dir]
#
# install_dir defaults to cloudy/build.  The data directory location is
# compiled into the executable, so do not move the build afterwards.
# Requires curl, make and a C++ compiler.
set -eu

version=c13.05
sha256=0369fab1bc231a44f1e81e15e345e40b8dcd6912d96b8f15d96d96053e3feae1
url=https://data.nublado.org/cloudy_releases/c13/$version.tar.gz

here=$(cd "$(dirname "$0")" && pwd)
dest=${1:-$here/build}
mkdir -p "$dest"
dest=$(cd "$dest" && pwd)
cd "$dest"

if [ ! -f $version.tar.gz ]; then
    echo "downloading $url"
    curl -fL -o $version.tar.gz.part "$url"
    mv $version.tar.gz.part $version.tar.gz
fi

if command -v sha256sum > /dev/null; then
    echo "$sha256  $version.tar.gz" | sha256sum -c -
else
    echo "$sha256  $version.tar.gz" | shasum -a 256 -c -
fi

if [ ! -d $version ]; then
    tar xzf $version.tar.gz
fi

ncpu=$(getconf _NPROCESSORS_ONLN 2> /dev/null || echo 2)
make -C $version/source -j"$ncpu"

echo
echo "Cloudy built.  To use it with ./run, set"
echo "  export CLOUDY_EXE=$dest/$version/source/cloudy.exe"
