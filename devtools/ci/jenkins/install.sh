# This script is executed via the line:
#   source devtools/ci/jenkins/install.sh
# in a bash shell with the -lex options turned on

echo "Using the following SWIG (`which swig`) version:"
swig -version

echo "Using cmake (`which cmake`) version":
cmake --version

echo "Using g++ (`which g++`) version:"
g++ --version

if [ -x `which nvcc 2>/dev/null || true` ]; then
    echo "Using nvcc (`which nvcc`) version:"
    nvcc --version
    export OPENMM_CUDA_COMPILER=`which nvcc`
fi

test -z "$INSTALL_DIRECTORY" && INSTALL_DIRECTORY="`pwd`/install"

# Build OpenMM
cmake -DCMAKE_INSTALL_PREFIX="${INSTALL_DIRECTORY}" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc \
      -DSWIG_EXECUTABLE=`which swig` -DCUDA_TOOLKIT_ROOT_DIR=${CUDA_HOME} $* .
make -j6 install
