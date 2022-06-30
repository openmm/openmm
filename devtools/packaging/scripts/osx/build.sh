#!/bin/bash

# Build script for Mac OS X distribution, for use in automated packaging.
# Note that this must be run from outside the checked-out openmm/ directory.

# Set relative workspace path.
export WORKSPACE=`pwd`

# Add conda binaries to path.
PATH=$WORKSPACE/miniconda/bin:$PATH

# Set install directory.
INSTALL=`pwd`/install
if [ -e $INSTALL ]; then
    rm -rf $INSTALL
fi

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$INSTALL"

CMAKE_FLAGS+=" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
CMAKE_FLAGS+=" -DCMAKE_OSX_DEPLOYMENT_TARGET=10.9"
CMAKE_FLAGS+=" -DCMAKE_OSX_SYSROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk"
CMAKE_FLAGS+=" -DOPENMM_GENERATE_API_DOCS=ON"

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

# Return to directory
cd $WORKSPACE
