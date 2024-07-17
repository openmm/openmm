source /opt/conda/etc/profile.d/conda.sh
set -eo pipefail

WORKSPACE="$HOME/workspace"

# This endgroup closes the open tag in CI.yml
echo "::endgroup::"

echo "::group::Prepare build environment..."
extra_conda_packages=""
if [[ ${COMPILERS} == devtoolset* ]]; then
    sudo yum install -y centos-release-scl
    sudo yum install -y ${COMPILERS}
    source /opt/rh/${COMPILERS}/enable
else
    extra_conda_packages="${COMPILERS}"
fi

# Patch environment file
sed -E -e "s/.*gromacs.*//" \
       -e "s/^- python$/- python ${PYTHON_VER}.*/" \
       ${WORKSPACE}/devtools/ci/gh-actions/conda-envs/build-ubuntu-latest.yml > conda-env.yml
for package in $extra_conda_packages; do
    if [[ -n ${package// } ]]; then
        echo "- ${package}" >> conda-env.yml
    fi
done
conda env create -n build -f conda-env.yml
conda activate build || true
echo "::endgroup::"

echo "::group::Prepare ccache..."
export CCACHE_BASEDIR=${WORKSPACE}
export CCACHE_DIR=${WORKSPACE}/.ccache
export CCACHE_COMPRESS=true
export CCACHE_COMPRESSLEVEL=6
export CCACHE_MAXSIZE=400M
ccache -p
ccache -z
echo "::endgroup::"

echo "::group::Configure with CMake..."
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
    -DCMAKE_C_COMPILER_LAUNCHER=ccache \
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
    -DOPENMM_BUILD_CUDA_TESTS=OFF \
    -DOPENMM_BUILD_OPENCL_TESTS=OFF
echo "::endgroup::"

# Build
echo "::group::Build with make..."
make -j2 install PythonInstall
echo "::endgroup::"

echo "::group::Check ccache performance..."
ccache -s
echo "::endgroup::"

# Core tests

echo "::group::Run core tests..."
python ${WORKSPACE}/devtools/run-ctest.py --parallel 2 --timeout 1500 --job-duration 360 --attempts 3
test -f ${CONDA_PREFIX}/lib/libOpenMM.so
test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMCPU.so
test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMPME.so
if [[ ! -z ${CUDA_VER} ]]; then
    test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMCUDA.so
    test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMOpenCL.so
fi
echo "::endgroup::"

# Python tests
echo "::group::Run Python tests..."
python -m openmm.testInstallation
python -c "import openmm as mm; print('---Loaded---', *mm.pluginLoadedLibNames, '---Failed---', *mm.Platform.getPluginLoadFailures(), sep='\n')"
cd python/tests
# Gromacs is not available on condaforge for PPC/ARM
# Membrane an MTS Langevin Integrator tests timeout (>6h!), possibly due to the emulation slowdown
python -m pytest -v -k "not gromacs and not membrane and not MTSLangevinIntegrator" -n 2
echo "::endgroup::"

echo "We are done!"
touch "${WORKSPACE}/docker_steps_run_successfully"
