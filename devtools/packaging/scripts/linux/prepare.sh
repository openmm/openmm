#!/bin/tcsh

# Prepare for build by ensuring necessary prerequisites are locally installed.

# Install miniconda
export VERSION="3.7.0"
export PLATFORM="Linux"
export ARCH="x86_64"
export MINICONDA="Miniconda-$VERSION-$PLATFORM-$ARCH.sh"
if [ -f $MINICONDA ];
then
   echo "File $MINICONDA exists, not downloading."
   export PATH=$WORKSPACE/miniconda/bin:$PATH
else
   echo "Downloading miniconda..."
   wget http://repo.continuum.io/miniconda/${MINICONDA}
   bash ${MINICONDA} -b -p miniconda
   PIP_ARGS="-U"
   export PATH=$WORKSPACE/miniconda/bin:$PATH
   conda config --add channels http://conda.binstar.org/omnia
   conda install --yes --quiet swig fftw3f pip
   pip install sphinxcontrib-bibtex
fi
