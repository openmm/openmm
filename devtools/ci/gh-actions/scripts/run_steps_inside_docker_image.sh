source /opt/conda/etc/profile.d/conda.sh
set -eo pipefail


WORKSPACE="$HOME/workspace"

echo "Prepare build environment..."
# Remove gromacs from dependencies
sed -E "s/.*gromacs.*//" ${WORKSPACE}/devtools/ci/gh-actions/conda-envs/build-ubuntu-latest.yml > conda-env.yml

conda create -y -n build python=${PYTHON_VER} compilers
conda env update -n build -f conda-env.yml
conda activate build || true

echo "Configure with CMake..."
if [[ -d /usr/local/cuda ]]; then
    export CUDA_PATH="/usr/local/cuda"
    export CUDA_LIB_PATH="${CUDA_PATH}/lib64/stubs"
    export LD_LIBRARY_PATH="${CUDA_PATH}/lib64/stubs:${LD_LIBRARY_PATH:-}"
    export PATH="${CUDA_PATH}/bin:${PATH}"
fi

rm -rf build || true
mkdir -p build
cd build
cmake ${WORKSPACE} \
    -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} \
    -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} \
    -DOPENMM_BUILD_CUDA_TESTS=OFF \
    -DOPENMM_BUILD_OPENCL_TESTS=OFF

# Build
echo "Build with make..."
make -j2 install PythonInstall

# Core tests
echo "Run core tests..."
python ${WORKSPACE}/devtools/run-ctest.py --parallel 2 --timeout 600 --job-duration 300
test -f ${CONDA_PREFIX}/lib/libOpenMM.so
test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMCPU.so
test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMPME.so
if [[ ! -z ${CUDA_VER} ]]; then
    test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMCUDA.so
    test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMOpenCL.so
fi

# Python tests
echo "Run Python tests..."
python -m simtk.testInstallation
python -c "import simtk.openmm as mm; print('---Loaded---', *mm.pluginLoadedLibNames, '---Failed---', *mm.Platform.getPluginLoadFailures(), sep='\n')"
cd python/tests
python -m pytest -v -n 2 -k "not gromacs"

echo "We are done!"
touch "${WORKSPACE}/docker_steps_run_successfully"
