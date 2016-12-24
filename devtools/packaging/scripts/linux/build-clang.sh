#!/bin/bash

# Build script for Linux distribution, for use in automated packaging.
# Note that this must be run from outside the checked-out openmm/ directory.

# Fix hbb issues.
# If statements needed because multiple Python versions are built in same docker image.
if [ ! -e /opt/rh/devtoolset-2/root/usr/lib/gcc/x86_64-redhat-linux ]; then
   ln -s /opt/rh/devtoolset-2/root/usr/lib/gcc/x86_64-CentOS-linux/ /opt/rh/devtoolset-2/root/usr/lib/gcc/x86_64-redhat-linux
fi
if [ ! -e /opt/rh/devtoolset-2/root/usr/include/c++/4.8.2/x86_64-redhat-linux ]; then
   ln -s /opt/rh/devtoolset-2/root/usr/include/c++/4.8.2/x86_64-CentOS-linux/ /opt/rh/devtoolset-2/root/usr/include/c++/4.8.2/x86_64-redhat-linux
fi

# Clang paths
export CLANG_PREFIX="/opt/clang"
export PATH=$PATH:$CLANG_PREFIX/bin

# enable devtoolset-2
source /opt/rh/devtoolset-2/enable

export MINIMAL_CFLAGS="-g -O3"
export CFLAGS="$MINIMAL_CFLAGS"
export CXXFLAGS="$MINIMAL_CFLAGS"
export LDFLAGS="$LDPATHFLAGS"

# Set relative workspace path.
export WORKSPACE=`pwd`

# Add conda binaries to path.
PATH=$WORKSPACE/miniconda/bin:$PATH

INSTALL=`pwd`/install
if [ -e $INSTALL ]; then
    rm -rf $INSTALL
fi

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$INSTALL"

# Don't build tests
CMAKE_FLAGS+=" -DBUILD_TESTING=OFF"

# Use clang
CMAKE_FLAGS+=" -DCMAKE_C_COMPILER=$CLANG_PREFIX/bin/clang -DCMAKE_CXX_COMPILER=$CLANG_PREFIX/bin/clang++"

# Ensure we build a release
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Release"

# setting the rpath so that libOpenMMPME.so finds the right libfftw3
#CMAKE_FLAGS+=" -DCMAKE_INSTALL_RPATH=.."
# Use NVIDIA CUDA 8.0
CMAKE_FLAGS+=" -DCUDA_CUDART_LIBRARY=/usr/local/cuda-8.0/lib64/libcudart.so"
CMAKE_FLAGS+=" -DCUDA_NVCC_EXECUTABLE=/usr/local/cuda-8.0/bin/nvcc"
CMAKE_FLAGS+=" -DCUDA_SDK_ROOT_DIR=/usr/local/cuda-8.0/"
CMAKE_FLAGS+=" -DCUDA_TOOLKIT_INCLUDE=/usr/local/cuda-8.0/include"
CMAKE_FLAGS+=" -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-8.0/"
# Use AMD APP SDK 3.0
CMAKE_FLAGS+=" -DOPENCL_INCLUDE_DIR=/opt/AMDAPPSDK-3.0/include/"
CMAKE_FLAGS+=" -DOPENCL_LIBRARY=/opt/AMDAPPSDK-3.0/lib/x86_64/libOpenCL.so"
# Generate API docs
CMAKE_FLAGS+=" -DOPENMM_GENERATE_API_DOCS=ON"
# Set location for FFTW3
PREFIX="$WORKSPACE/miniconda"
CMAKE_FLAGS+=" -DFFTW_INCLUDES=$PREFIX/include"
CMAKE_FLAGS+=" -DFFTW_LIBRARY=$PREFIX/lib/libfftw3f.so"
CMAKE_FLAGS+=" -DFFTW_THREADS_LIBRARY=$PREFIX/lib/libfftw3f_threads.so"
# Necessary to find GL headers
CMAKE_FLAGS+=" -DCMAKE_CXX_FLAGS_RELEASE=-I/usr/include/nvidia/"

# Build in subdirectory.
if [ -e build ]; then
    rm -rf build
fi
mkdir build
cd build
cmake ../openmm $CMAKE_FLAGS
make -j4 all install

export CFLAGS="$MINIMAL_CFLAGS"
export CXXFLAGS="$MINIMAL_CFLAGS"
export LDFLAGS="$LDPATHFLAGS"

make -j4 PythonInstall C++ApiDocs PythonApiDocs sphinxpdf

# Install.
make install

cd ..
