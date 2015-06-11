#!/bin/bash

# Test Mac OS X distribution package.
# Package must have already been downloaded.

# Set relative workspace path.
export WORKSPACE=`pwd`

# Add conda binaries to path.
PATH=$WORKSPACE/miniconda/bin:$PATH

# Get Python executable.
PYTHON=`which python`

# Set install directory.
INSTALL=`pwd`/install
if [ -e $INSTALL ]; then
    rm -rf $INSTALL
fi

# Unzip package
unzip OpenMM-*.zip

# Install
bash -x install.sh $INSTALL $PYTHON

# Test.
python -m simtk.testInstallation
