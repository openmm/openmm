cmake_minimum_required(VERSION 3.18)
project(OpenMMCudaTensorBridge LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import torch; print(torch.utils.cmake_prefix_path)"
    OUTPUT_VARIABLE TORCH_CMAKE_PREFIX_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
list(APPEND CMAKE_PREFIX_PATH "${TORCH_CMAKE_PREFIX_PATH}")
find_package(Torch REQUIRED)

find_package(CUDAToolkit REQUIRED)

set(OPENMM_LM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH "Path to the OpenMM-LM source root")
set(OPENMM_INCLUDE_DIR "${OPENMM_LM_ROOT}/openmmapi/include" CACHE PATH "OpenMM API include directory")
set(OPENMM_CUDA_INCLUDE_DIR "${OPENMM_LM_ROOT}/platforms/cuda/include" CACHE PATH "OpenMM CUDA include directory")
set(OPENMM_COMMON_INCLUDE_DIR "${OPENMM_LM_ROOT}/platforms/common/include" CACHE PATH "OpenMM common include directory")
set(OPENMM_LIBRARY OpenMM CACHE STRING "OpenMM shared library target or linker name")
set(OPENMM_CUDA_LIBRARY OpenMMCUDA CACHE STRING "OpenMM CUDA shared library target or linker name")

pybind11_add_module(openmm_cuda_bridge MODULE bridge.cpp bridge.cu)

target_include_directories(
    openmm_cuda_bridge
    PRIVATE
        "${OPENMM_INCLUDE_DIR}"
        "${OPENMM_INCLUDE_DIR}/openmm"
        "${OPENMM_INCLUDE_DIR}/openmm/internal"
        "${OPENMM_CUDA_INCLUDE_DIR}"
        "${OPENMM_COMMON_INCLUDE_DIR}"
)

target_link_libraries(
    openmm_cuda_bridge
    PRIVATE
        "${TORCH_LIBRARIES}"
        CUDA::cudart
        CUDA::cuda_driver
        "${OPENMM_LIBRARY}"
        "${OPENMM_CUDA_LIBRARY}"
)

target_compile_definitions(openmm_cuda_bridge PRIVATE TORCH_API_INCLUDE_EXTENSION_H)
target_compile_options(openmm_cuda_bridge PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)

set_target_properties(
    openmm_cuda_bridge
    PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        BUILD_RPATH_USE_ORIGIN ON
        INSTALL_RPATH "$ORIGIN"
)
