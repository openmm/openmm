# This script is executed via the line:
#   source devtools/ci/jenkins/install.sh
# in a bash shell with the -lex options turned on

echo "I am `whoami`"

conda install -y cmake numpy scipy pytest swig

echo "Using the following SWIG (`which swig`) version:"
swig -version

echo "Using cmake (`which cmake`) version":
cmake --version

echo "Using g++ (`which g++`) version:"
g++ --version

echo "Using nvcc (`which nvcc`) version:"
nvcc --version

export OPENMM_CUDA_COMPILER=`which nvcc`

# Constants
INSTALL_DIRECTORY="`pwd`/install"

# Build OpenMM
cmake -DCMAKE_INSTALL_PREFIX="${INSTALL_DIRECTORY}" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc .
make -j6 install
make PythonInstall

# Now run the tests
python -m simtk.testInstallation
cd python/tests && py.test -v && cd ../..
python devtools/run-ctest.py --job-duration=120 --timeout 300 --parallel=2
