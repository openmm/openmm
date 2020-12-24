set -euxo pipefail
WORKSPACE="$HOME/workspace"

# Remove gromacs from dependencies
sed -E "s/.*gromacs.*//" ${WORKSPACE}/devtools/ci/gh-actions/conda-envs/build-ubuntu-latest.yml > conda-env.yml

conda env create -n build -f conda-env.yml
conda activate build

CMAKE_FLAGS=""
if [[ ! -z ${CUDA_VER} ]]; then
    export CUDA_PATH="/usr/local/cuda-${CUDA_VER}"
    CMAKE_FLAGS+=" -DOPENCL_LIBRARY=${CUDA_PATH}/lib64/libOpenCL.so"
    CMAKE_FLAGS+=" -DCUDA_CUDART_LIBRARY=${CUDA_PATH}/lib64/libcudart.so"
    CMAKE_FLAGS+=" -DCUDA_NVCC_EXECUTABLE=${CUDA_PATH}/bin/nvcc"
    CMAKE_FLAGS+=" -DCUDA_SDK_ROOT_DIR=${CUDA_PATH}/"
    CMAKE_FLAGS+=" -DCUDA_TOOLKIT_INCLUDE=${CUDA_PATH}/include"
    CMAKE_FLAGS+=" -DCUDA_TOOLKIT_ROOT_DIR=${CUDA_PATH}/"
fi

mkdir build
cd build
cmake ${WORKSPACE} \
    -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} \
    -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} \
    -DOPENMM_BUILD_CUDA_TESTS=OFF \
    -DOPENMM_BUILD_OPENCL_TESTS=OFF \
    ${CMAKE_FLAGS}

# Build
make -j2 install PythonInstall

# Core tests
python ${WORKSPACE}/devtools/run-ctest.py --parallel 2 --timeout 600 --job-duration 300
test -f ${CONDA_PREFIX}/lib/libOpenMM.so
test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMCPU.so
test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMPME.so
if [[ ! -z ${CUDA_VER} ]]; then
    test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMCUDA.so
    test -f ${CONDA_PREFIX}/lib/plugins/libOpenMMOpenCL.so
fi

# Python tests
python -m simtk.testInstallation
python -c "import simtk.openmm as mm; print('---Loaded---', *mm.pluginLoadedLibNames, '---Failed---', *mm.Platform.getPluginLoadFailures(), sep='\n')"
cd python/tests
python -m pytest -v -n 2 -k "not gromacs"

touch "${WORKSPACE}/docker_steps_run_successfully"
