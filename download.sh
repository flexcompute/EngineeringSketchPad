#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error
set -x # echo commands

version=${1:-129}
tar=ESP${version}-linux-x86_64.tgz

# MIT rotates older releases from PreBuilts/ to archive/
prebuilts_url=https://acdl.mit.edu/ESP/PreBuilts/${tar}
archive_url=https://acdl.mit.edu/ESP/archive/${tar}

if curl --head --silent --fail "$prebuilts_url" > /dev/null 2>&1; then
  curl -O "$prebuilts_url" -o ${tar}
elif curl --head --silent --fail "$archive_url" > /dev/null 2>&1; then
  echo "Not found in PreBuilts, downloading from archive"
  curl -O "$archive_url" -o ${tar}
else
  echo "ERROR: ESP${version} not found in PreBuilts or archive" >&2
  exit 1
fi

rm -rf ESP ESP${version}
tar xzf ${tar}

# Version-specific source patches
if [ "$version" = "128" ]; then
  echo "Version is 128: patching evaluate.c (periodic parameter support for EG_splinePCDeriv)"
  cp replace128evaluate.c ESP${version}/EngSketchPad/src/EGADS/util/evaluate.c
fi

mv ESP${version} ESP
