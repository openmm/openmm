#!/bin/bash

# Build script for Linux distribution, for use in automated packaging.
# Note that this must be run from outside the checked-out openmm/ directory.

# Add conda binaries to path.
PATH=${HOME}/miniconda/bin:${PATH}

INSTALL=`pwd`/install
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$INSTALL"

# setting the rpath so that libOpenMMPME.so finds the right libfftw3
#CMAKE_FLAGS+=" -DCMAKE_INSTALL_RPATH=.."
CMAKE_FLAGS+=" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
CMAKE_FLAGS+=" -DCUDA_CUDART_LIBRARY=/usr/local/cuda-6.5/lib64/libcudart.so"
CMAKE_FLAGS+=" -DCUDA_NVCC_EXECUTABLE=/usr/local/cuda-6.5/bin/nvcc"
CMAKE_FLAGS+=" -DCUDA_SDK_ROOT_DIR=/usr/local/cuda-6.5/"
CMAKE_FLAGS+=" -DCUDA_TOOLKIT_INCLUDE=/usr/local/cuda-6.5/include"
CMAKE_FLAGS+=" -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-6.5/"
CMAKE_FLAGS+=" -DOPENCL_INCLUDE_DIR=/usr/local/cuda-6.5/include"
CMAKE_FLAGS+=" -DOPENCL_LIBRARY=/usr/local/cuda-6.5/lib64/libOpenCL.so"

# Set location for FFTW3
PREFIX="${HOME}/miniconda"
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
make -j4 all DoxygenApiDocs sphinxpdf
make install

# Install Python wrappers.
OPENMM_INCLUDE_PATH=$INSTALL/include
OPENMM_LIB_PATH=$INSTALL/lib
cd python
python setup.py install --prefix=$INSTALL
cd ..

# Copy all tests to bin directory so they will be distributed with install package.
#cp `find . -name "Test*" -type f -maxdepth 1` $PREFIX/bin
