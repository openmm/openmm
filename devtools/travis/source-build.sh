#!/bin/bash
set -e

source devtools/packaging/scripts/source/prepare.sh
source devtools/packaging/scripts/source/build.sh
source devtools/packaging/scripts/source/package.sh

# Copy packages out
cp -r packaging/* /io
