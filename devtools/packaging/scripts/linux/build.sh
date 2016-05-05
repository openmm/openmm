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
CMAKE_FLAGS+=" -DCUDA_CUDART_LIBRARY=/usr/local/cuda-7.5/lib64/libcudart.so"
CMAKE_FLAGS+=" -DCUDA_NVCC_EXECUTABLE=/usr/local/cuda-7.5/bin/nvcc"
CMAKE_FLAGS+=" -DCUDA_SDK_ROOT_DIR=/usr/local/cuda-7.5/"
CMAKE_FLAGS+=" -DCUDA_TOOLKIT_INCLUDE=/usr/local/cuda-7.5/include"
CMAKE_FLAGS+=" -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-7.5/"
CMAKE_FLAGS+=" -DOPENCL_INCLUDE_DIR=/opt/AMDAPPSDK-3.0/include/"
CMAKE_FLAGS+=" -DOPENCL_LIBRARY=/opt/AMDAPPSDK-3.0/lib/x86_64/libOpenCL.so"
CMAKE_FLAGS+=" -DOPENMM_GENERATE_API_DOCS=ON"

# Set location for FFTW3
PREFIX="$WORKSPACE/miniconda"
CMAKE_FLAGS+=" -DFFTW_INCLUDES=$PREFIX/include"
CMAKE_FLAGS+=" -DFFTW_LIBRARY=$PREFIX/lib/libfftw3f.so"
CMAKE_FLAGS+=" -DFFTW_THREADS_LIBRARY=$PREFIX/lib/libfftw3f_threads.so"

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
