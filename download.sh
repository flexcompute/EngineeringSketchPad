#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error
set -x # echo commands

version=${1:-129}
tar=ESP${version}-linux-x86_64.tgz

# A given ESP version lives in exactly ONE of PreBuilts/ or archive/: MIT moves
# a release from PreBuilts/ to archive/ once it ages out, so a 404 from the
# other directory is EXPECTED and must not abort the build. Reaching
# acdl.mit.edu can also fail transiently (connection timeouts, 5xx) -- that is
# what we retry. curl's --retry covers those transient errors but never a 404,
# and --connect-timeout keeps a hung endpoint from stalling the job for minutes.
prebuilts_url=https://acdl.mit.edu/ESP/PreBuilts/${tar}
archive_url=https://acdl.mit.edu/ESP/archive/${tar}

mit_curl() {
  curl --fail --location --connect-timeout 30 \
       --retry 5 --retry-delay 5 --retry-connrefused "$@"
}

if mit_curl --silent --head --max-time 120 "$prebuilts_url" > /dev/null; then
  url=$prebuilts_url
  echo "Found ESP${version} in PreBuilts"
elif mit_curl --silent --head --max-time 120 "$archive_url" > /dev/null; then
  url=$archive_url
  echo "Not in PreBuilts (expected); found ESP${version} in archive"
else
  echo "ERROR: ESP${version} unreachable in PreBuilts or archive after retries" >&2
  exit 1
fi

mit_curl --silent --show-error --max-time 1800 -o "${tar}" "$url"

rm -rf ESP ESP${version}
tar xzf ${tar}
mv ESP${version} ESP

# Version-specific source patches
if [ "$version" = "128" ]; then
  echo "Patching ESP128: evaluate.c (periodic parameter support for EG_splinePCDeriv)"
  cp replace128evaluate.c ESP/EngSketchPad/src/EGADS/util/evaluate.c
fi
if [ "$version" = "129" ]; then
  echo "Patching ESP129: EGADS try/catch for OCC calls (FXC-6881), null surface checks, OuterShell guard (FXC-8846)"
  patch -p0 < egadsTopo129.patch
  patch -p0 < egadsIO129.patch
  patch -p0 < egadsTopoOuterShell129.patch
fi
