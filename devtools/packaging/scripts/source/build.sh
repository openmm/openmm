#!/bin/bash

# Build script for Linux distribution, for use in automated packaging.
# Note that this must be run from outside the checked-out openmm/ directory.

# Set relative workspace path.
export WORKSPACE=`pwd`

# Add conda binaries to path.
PATH=$WORKSPACE/miniconda/bin:$PATH

INSTALL=`pwd`/install
if [ -e $INSTALL ]; then
    rm -rf $INSTALL
fi

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$INSTALL"

# setting the rpath so that libOpenMMPME.so finds the right libfftw3
#CMAKE_FLAGS+=" -DCMAKE_INSTALL_RPATH=.."
CMAKE_FLAGS+=" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
CMAKE_FLAGS+=" -DOPENMM_BUILD_AMOEBA_CUDA_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_CPU_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_CUDA_COMPILER_PLUGIN=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_CUDA_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_DRUDE_CUDA_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_DRUDE_OPENCL_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_OPENCL_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_RPMD_CUDA_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_RPMD_OPENCL_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_GENERATE_API_DOCS=ON"

# Set location for FFTW3
#PREFIX="$WORKSPACE/miniconda"
#CMAKE_FLAGS+=" -DFFTW_INCLUDES=$PREFIX/include"
#CMAKE_FLAGS+=" -DFFTW_LIBRARY=$PREFIX/lib/libfftw3f.so"
#CMAKE_FLAGS+=" -DFFTW_THREADS_LIBRARY=$PREFIX/lib/libfftw3f_threads.so"

# Build in subdirectory.
if [ -e build ]; then
    rm -rf build
fi
mkdir build
cd build
cmake ../openmm $CMAKE_FLAGS
make -j4 all install
make -j4 PythonInstall C++ApiDocs PythonApiDocs sphinxpdf

# Install.
make install
