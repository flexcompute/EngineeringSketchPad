#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error
set -x # echo commands

version=128
tar=ESP${version}-linux-x86_64.tgz

curl -O https://acdl.mit.edu/ESP/PreBuilts/${tar} -o ${tar}

rm -rf ESP ESP${version}
tar xzf ${tar}
mv ESP${version} ESP
