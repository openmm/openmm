#!/bin/tcsh

# Prepare for build by ensuring necessary prerequisites are locally installed.

# Set relative workspace path.
export WORKSPACE=`pwd`

# Install miniconda
export VERSION="Latest"
export PLATFORM="MacOSX"
export ARCH="x86_64"
export MINICONDA="Miniconda-$VERSION-$PLATFORM-$ARCH.sh"
if [ -f miniconda ];
then
   echo "miniconda already exists"
else
   echo "Downloading miniconda..."
   rm -rf Miniconda-*
   wget --quiet http://repo.continuum.io/miniconda/${MINICONDA}
   bash ${MINICONDA} -b -p miniconda
   PIP_ARGS="-U"
fi

# Add to path.
export PATH=$WORKSPACE/miniconda/bin:$PATH

# Ensure configuration is up to date.
conda config --add channels http://conda.binstar.org/omnia
conda install --yes --quiet swig fftw3f pip
pip install sphinxcontrib-bibtex sphinxcontrib-lunrsearch sphinxcontrib-autodoc_doxygen
