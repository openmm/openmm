#!/bin/bash

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$PREFIX"

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    # setting the rpath so that libOpenMMPME.so finds the right libfftw3
    CMAKE_FLAGS+=" -DCMAKE_INSTALL_RPATH=.."
    CMAKE_FLAGS+=" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
    CMAKE_FLAGS+=" -DOPENCL_LIBRARY=/opt/AMDAPP/lib/x86_64/libOpenCL.so" # TEST

elif [[ "$OSTYPE" == "darwin"* ]]; then
    export MACOSX_DEPLOYMENT_TARGET="10.7"
    CMAKE_FLAGS+=" -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
fi

# Set location for FFTW3 on both linux and mac
CMAKE_FLAGS+=" -DFFTW_INCLUDES=$PREFIX/include"
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    CMAKE_FLAGS+=" -DFFTW_LIBRARY=$PREFIX/lib/libfftw3f.so"
    CMAKE_FLAGS+=" -DFFTW_THREADS_LIBRARY=$PREFIX/lib/libfftw3f_threads.so"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    CMAKE_FLAGS+=" -DFFTW_LIBRARY=$PREFIX/lib/libfftw3f.dylib"
    CMAKE_FLAGS+=" -DFFTW_THREADS_LIBRARY=$PREFIX/lib/libfftw3f_threads.dylib"
fi

# Copy source to current directory.
cp -r $RECIPE_DIR/../.. .

# Build in subdirectory.
mkdir build
cd build
cmake .. $CMAKE_FLAGS
make -j4
make install

# Run C tests.
# Exclude OpenCL tests because @peastman suspects mesa on travis implementation is broken.
# @jchodera and @pgrinaway suspect travis is working, but AMD OpenCL tests are actually failing due to a bug.
#ctest -j2 -V -E "[A-Za-z]+OpenCL[A-Za-z]+"

# Install Python wrappers.
export OPENMM_INCLUDE_PATH=$PREFIX/include
export OPENMM_LIB_PATH=$PREFIX/lib
cd python
$PYTHON setup.py install
cd ..

# Remove one random file
#rm $PREFIX/bin/TestReferenceHarmonicBondForce

# Copy all tests to bin directory so they will be distributed with install package.
cp `find . -name "Test*" -type f -maxdepth 1` $PREFIX/bin

