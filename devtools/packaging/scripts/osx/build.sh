#!/bin/bash

# Build script for Mac OS X distribution, for use in automated packaging.
# Note that this must be run from outside the checked-out openmm/ directory.

INSTALL=`pwd`/install
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$INSTALL"

# setting the rpath so that libOpenMMPME.so finds the right libfftw3
#CMAKE_FLAGS+=" -DCMAKE_INSTALL_RPATH=.."
CMAKE_FLAGS+=" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
CMAKE_FLAGS+=" -DCMAKE_OSX_DEPLOYMENT_TARGET=10.9"
CMAKE_FLAGS+=" -DCMAKE_OSX_SYSROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk"
CMAKE_FLAGS+=" -DOPENMM_BUILD_OPENCL_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_DRUDE_OPENCL_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_RPMD_OPENCL_LIB=OFF"
CMAKE_FLAGS+=" -DOPENMM_BUILD_OPENCL_TESTS=FALSE"
CMAKE_FLAGS+=" -DOPENMM_BUILD_OPENCL_DOUBLE_PRECISION_TESTS=FALSE"

# Build in subdirectory.
# Set location for FFTW3
PREFIX="${HOME}/miniconda"
CMAKE_FLAGS+=" -DFFTW_INCLUDES=$PREFIX/include"
CMAKE_FLAGS+=" -DFFTW_LIBRARY=$PREFIX/lib/libfftw3f.so"
CMAKE_FLAGS+=" -DFFTW_THREADS_LIBRARY=$PREFIX/lib/libfftw3f_threads.so"

mkdir build
cd build
cmake ../openmm $CMAKE_FLAGS
make -j4
make install

# Install Python wrappers.
export OPENMM_INCLUDE_PATH=$INSTALL/include
export OPENMM_LIB_PATH=$INSTALL/lib
cd python
$PYTHON setup.py install --prefix=$INSTALL
cd ..

# Copy all tests to bin directory so they will be distributed with install package.
#cp `find . -name "Test*" -type f -maxdepth 1` $PREFIX/bin
