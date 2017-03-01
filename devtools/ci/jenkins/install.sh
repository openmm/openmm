# This script is executed via the line:
#   source devtools/ci/jenkins/install.sh
# in a bash shell with the -lex options turned on

echo "Using the following SWIG (`which swig`) version:"
swig -version

echo "Using cmake (`which cmake`) version":
cmake --version

echo "Using clang (`which clang`) version:"
clang --version

module load cuda conda/jenkins

# Constants
CONDAENV=openmm-test-3.5
INSTALL_DIRECTORY="${WORKSPACE}/openmm-install"
SRC_DIRECTORY="${WORKSPACE}/openmm-src" # set in the Jenkins configuration
export OPENMM_CUDA_COMPILER=`which nvcc`

# Create a conda environment, but clean up after one first. If it doesn't exist, don't complain.
# But since we are invoking this shell with -e (exit on all errors), we need || true to prevent this
# command from crashing the whole shell
conda remove -yn ${CONDAENV} --all --quiet || true
conda create -yn ${CONDAENV} python=3.5 --no-default-packages --quiet
conda install -yn ${CONDAENV} numpy scipy pytest --quiet
source activate ${CONDAENV} # enter our new environment

# Build OpenMM
cd "${SRC_DIRECTORY}"
cmake -DCMAKE_INSTALL_PREFIX="${INSTALL_DIRECTORY}" -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang .
make -j4 install
make PythonInstall

# Now run the tests
python -m simtk.testInstallation
cd python/tests && py.test -v && cd ../..
python devtools/run-ctest.py --job-duration=120 --timeout 300

# Now remove the conda environment
source deactivate
conda remove -yn ${CONDAENV} --all --quiet