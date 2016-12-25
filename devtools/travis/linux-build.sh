#!/bin/bash
set -e

ln -s /io openmm

source openmm/devtools/packaging/scripts/linux/prepare.sh
source openmm/devtools/packaging/scripts/linux/build.sh
source openmm/devtools/packaging/scripts/linux/package.sh

# Copy packages out
cp -r packaging/* /io
