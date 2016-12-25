#!/bin/bash
set -e

ln -s /io openmm

source openmm/devtools/packaging/scripts/source/prepare.sh
source openmm/devtools/packaging/scripts/source/build.sh
source openmm/devtools/packaging/scripts/source/package.sh

# Copy packages out
cp -r packaging/* /io
