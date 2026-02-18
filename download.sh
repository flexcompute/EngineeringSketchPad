#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error
set -x # echo commands

version=128
tar=ESP${version}-linux-x86_64.tgz

curl -O https://acdl.mit.edu/ESP/PreBuilts/${tar} -o ${tar}

rm -rf ESP ESP${version}
tar xzf ${tar}

# Replace evaluate.c with a custom version for ESP128 to fix
# add periodic parameter support to EG_splinePCDeriv
if [ "$version" = "128" ]; then
  echo "Version is 128: copying replace128evaluate.c to ESP${version}/EngSketchPad/src/EGADS/util/evaluate.c"
  cp replace128evaluate.c ESP${version}/EngSketchPad/src/EGADS/util/evaluate.c
else
  echo "Version is ${version} (not 128): skipping replace128evaluate.c copy"
fi

mv ESP${version} ESP
