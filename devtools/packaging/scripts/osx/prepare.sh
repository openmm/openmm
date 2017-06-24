#!/bin/tcsh

# Prepare for build by ensuring necessary prerequisites are locally installed.

# Set relative workspace path.
export WORKSPACE=`pwd`

# Install miniconda
export VERSION="latest"
export PLATFORM="MacOSX"
export ARCH="x86_64"
export MINICONDA="Miniconda3-$VERSION-$PLATFORM-$ARCH.sh"
if [ -f $WORKSPACE/miniconda ];
then
   echo "miniconda already exists"
else
   echo "Downloading miniconda..."
   rm -rf $WORKSPACE/Miniconda3-*
   wget https://repo.continuum.io/miniconda/${MINICONDA}
   bash ${MINICONDA} -b -p $WORKSPACE/miniconda
   PIP_ARGS="-U"
fi

# Add to path.
export PATH=$WORKSPACE/miniconda/bin:$PATH

# Ensure configuration is up to date.
conda config --add channels http://conda.binstar.org/omnia
conda install --yes --quiet swig fftw3f pip doxygen sphinx sphinxcontrib-bibtex sphinxcontrib-lunrsearch sphinxcontrib-autodoc_doxygen lxml cmake
pip install sphinxcontrib-bibtex sphinxcontrib-lunrsearch sphinxcontrib-autodoc_doxygen
