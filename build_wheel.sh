#!/bin/bash
set -e
mkdir build
cd build
cmake3 .. \
    -DCMAKE_INSTALL_PREFIX=$(pwd)/python/openmm \
    -DCMAKE_PREFIX_PATH=$(pwd)/python/openmm \
    -DOPENMM_BUILD_CUDA_TESTS=OFF \
    -DOPENMM_BUILD_OPENCL_TESTS=OFF \
    -DOPENMM_BUILD_CPU_LIB=ON \
    -DOPENMM_BUILD_PYTHON_WRAPPERS=ON \
    -DOPENMM_BUILD_REFERENCE_TESTS=OFF \
    -DOPENMM_BUILD_SERIALIZATION_TESTS=OFF \
    -DOPENMM_BUILD_C_AND_FORTRAN_WRAPPERS=OFF \
    -DOPENMM_BUILD_EXAMPLES=OFF \
    -DOPENMM_BUILD_UNVERSIONED=ON \
    -DOPENCL_LIBRARY=/usr/local/cuda-11.3/lib64/libOpenCL.so \
    -DCUDA_CUDART_LIBRARY=/usr/local/cuda-11.3/lib64/libcudart.so \
    -DCUDA_NVCC_EXECUTABLE=/usr/local/cuda-11.3/bin/nvcc \
    -DCUDA_SDK_ROOT_DIR=/usr/local/cuda-11.3/ \
    -DCUDA_TOOLKIT_INCLUDE=/usr/local/cuda-11.3/include \
    -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11.3/ \

export OPENMM_MAKE_WHEEL=1

make -j24
make install

cd python
export OPENMM_INCLUDE_PATH=$(pwd)/openmm/include
export OPENMM_LIB_PATH=$(pwd)/openmm/lib
python setup.py bdist_wheel