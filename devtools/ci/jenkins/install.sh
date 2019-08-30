#!/bin/bash -ex
# This script is executed via the line:
#   source devtools/ci/jenkins/install.sh
# in a bash shell with the -lex options turned on

echo "Using the following SWIG (`which swig`) version:"
swig -version

echo "Using cmake (`which cmake`) version":
cmake --version

echo "Using g++ (`which g++`) version:"
g++ --version

if [ ! -z "$OPENMM_CUDA_COMPILER" ]; then
    echo "Using nvcc ($OPENMM_CUDA_COMPILER) version:"
    $OPENMM_CUDA_COMPILER --version
    CUDA_ARGS="-DCUDA_TOOLKIT_ROOT_DIR=${CUDA_HOME} -DOPENMM_BUILD_CUDA_LIB=true"
fi

cmake -DCMAKE_INSTALL_PREFIX="`pwd`/install" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc \
      -DSWIG_EXECUTABLE=`which swig` $CUDA_ARGS $EXTRA_CMAKE_ARGS .
make -j6 install
