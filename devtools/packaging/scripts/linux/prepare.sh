#!/bin/bash

# Prepare for build by ensuring necessary prerequisites are locally installed.

# Set relative workspace path.
export WORKSPACE=`pwd`

# Install miniconda
export VERSION="latest"
export PLATFORM="Linux"
export ARCH="x86_64"
export MINICONDA="Miniconda2-$VERSION-$PLATFORM-$ARCH.sh"
if [ -f miniconda ];
then
   echo "miniconda already exists"
else
   echo "Downloading miniconda..."
   rm -rf Miniconda-* miniconda ~/.condarc
   wget --quiet http://repo.continuum.io/miniconda/${MINICONDA}
   bash ${MINICONDA} -b -p miniconda
   PIP_ARGS="-U"
fi

# Add to path.
export PATH=$WORKSPACE/miniconda/bin:$PATH

# Ensure configuration is up to date.
conda config --add channels omnia
conda install --yes --quiet swig fftw3f pip doxygen sphinx sphinxcontrib-bibtex sphinxcontrib-lunrsearch sphinxcontrib-autodoc_doxygen lxml cmake

