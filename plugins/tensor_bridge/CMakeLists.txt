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
    COMMAND "${Python_EXECUTABLE}" -c "import torch, os; print(os.path.join(os.path.dirname(torch.__file__), 'lib'))"
    OUTPUT_VARIABLE TORCH_PY_LIBDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
list(PREPEND CMAKE_PREFIX_PATH "${TORCH_PY_LIBDIR}")

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import torch; print(torch.utils.cmake_prefix_path)"
    OUTPUT_VARIABLE TORCH_CMAKE_PREFIX_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(DEFINED ENV{CUDA_HOME} AND NOT "$ENV{CUDA_HOME}" STREQUAL "")
    list(PREPEND CMAKE_PREFIX_PATH "$ENV{CUDA_HOME}")
endif()
list(APPEND CMAKE_PREFIX_PATH "${TORCH_CMAKE_PREFIX_PATH}")
find_package(Torch REQUIRED)

find_library(TORCH_PYTHON_LIBRARY NAMES torch_python PATHS "${TORCH_PY_LIBDIR}" REQUIRED NO_DEFAULT_PATH)

find_package(CUDAToolkit REQUIRED)

set(OPENMM_LM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH "Path to the OpenMM source root")
set(OPENMM_BUILD_DIR "" CACHE PATH "OpenMM CMake build directory (must contain libOpenMM.so and libOpenMMCUDA.so)")
set(OPENMM_INCLUDE_DIR "${OPENMM_LM_ROOT}/openmmapi/include" CACHE PATH "OpenMM API include directory")
set(OPENMM_CUDA_INCLUDE_DIR "${OPENMM_LM_ROOT}/platforms/cuda/include" CACHE PATH "OpenMM CUDA include directory")
set(OPENMM_COMMON_INCLUDE_DIR "${OPENMM_LM_ROOT}/platforms/common/include" CACHE PATH "OpenMM common include directory")

if(NOT OPENMM_BUILD_DIR)
    message(FATAL_ERROR "OPENMM_BUILD_DIR must be set to the OpenMM build tree (same dir as libOpenMM.so).")
endif()
if(NOT EXISTS "${OPENMM_BUILD_DIR}/libOpenMM.so")
    message(FATAL_ERROR "libOpenMM.so not found under OPENMM_BUILD_DIR='${OPENMM_BUILD_DIR}'")
endif()
if(NOT EXISTS "${OPENMM_BUILD_DIR}/libOpenMMCUDA.so")
    message(FATAL_ERROR "libOpenMMCUDA.so not found under OPENMM_BUILD_DIR='${OPENMM_BUILD_DIR}'")
endif()

add_library(OpenMM SHARED IMPORTED)
set_target_properties(OpenMM PROPERTIES IMPORTED_LOCATION "${OPENMM_BUILD_DIR}/libOpenMM.so")
add_library(OpenMMCUDA SHARED IMPORTED)
set_target_properties(OpenMMCUDA PROPERTIES IMPORTED_LOCATION "${OPENMM_BUILD_DIR}/libOpenMMCUDA.so")

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
        "${TORCH_PYTHON_LIBRARY}"
        CUDA::cudart
        CUDA::cuda_driver
        "${CUDAToolkit_LIBRARY_DIR}/libcudadevrt.a"
        OpenMM
        OpenMMCUDA
)

target_compile_definitions(openmm_cuda_bridge PRIVATE TORCH_API_INCLUDE_EXTENSION_H)
target_compile_options(openmm_cuda_bridge PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)

set_target_properties(
    openmm_cuda_bridge
    PROPERTIES
        CUDA_SEPARABLE_COMPILATION OFF
        BUILD_RPATH_USE_ORIGIN ON
        BUILD_RPATH "${OPENMM_BUILD_DIR}"
        INSTALL_RPATH "${OPENMM_BUILD_DIR};$ORIGIN"
)
