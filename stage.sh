#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error
set -x # echo commands

esp=ESP

rm -rf include lib

mkdir -p lib
cp ${esp}/OpenCASCADE-*/lib/lib* lib

export ESP_ROOT=${PWD}/ESP/EngSketchPad
export ESP_ARCH=LINUX64
export CASREV=7.8
export CASROOT=${PWD}/ESP/OpenCASCADE-${CASREV}.1

export RPATH=-Wl,-rpath='$$ORIGIN:$$ORIGIN/../lib',--disable-new-dtags

( cd ${esp}/EngSketchPad/src/EGADS/src && make clean && \
  make RPATH=\'${RPATH}\' )

cp ${esp}/EngSketchPad/lib/* lib
cp -r ${esp}/EngSketchPad/include include
