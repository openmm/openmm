#!/bin/bash
set -e

source devtools/packaging/scripts/linux/prepare.sh
source devtools/packaging/scripts/linux/build.sh
source devtools/packaging/scripts/linux/package.sh

# Copy packages out
cp -r packaging/* /io
