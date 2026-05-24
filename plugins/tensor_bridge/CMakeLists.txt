cmake_minimum_required(VERSION 3.18)
project(OpenMMCudaTensorBridge LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(GENAI_TPS_STABLE_TORCH_ABI "Build openmm_cuda_bridge with LibTorch Stable ABI" ON)

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
if(DEFINED CUDAToolkit_ROOT)
    list(PREPEND CMAKE_PREFIX_PATH "${CUDAToolkit_ROOT}")
endif()
list(APPEND CMAKE_PREFIX_PATH "${TORCH_CMAKE_PREFIX_PATH}")
find_package(Torch REQUIRED)

find_library(TORCH_PYTHON_LIBRARY NAMES torch_python PATHS "${TORCH_PY_LIBDIR}" REQUIRED NO_DEFAULT_PATH)

# conda-forge pytorch ships a minimal Stable ABI shim (no stack/mul_Scalar/cat). Fall back to
# torch/extension.h when required symbols are absent so pixi builds succeed out of the box.
if(GENAI_TPS_STABLE_TORCH_ABI)
    set(_genai_tps_shim_header "")
    foreach(_inc_dir IN ITEMS ${TORCH_INCLUDE_DIRS} "$ENV{CONDA_PREFIX}/include")
        if(EXISTS "${_inc_dir}/torch/csrc/inductor/aoti_torch/generated/c_shim_aten.h")
            set(_genai_tps_shim_header "${_inc_dir}/torch/csrc/inductor/aoti_torch/generated/c_shim_aten.h")
            break()
        endif()
    endforeach()
    if(_genai_tps_shim_header)
        file(READ "${_genai_tps_shim_header}" _genai_tps_shim_contents)
        if(NOT _genai_tps_shim_contents MATCHES "aoti_torch_aten_stack")
            message(STATUS "openmm_cuda_bridge: incomplete LibTorch Stable ABI shim; using legacy extension ABI")
            set(GENAI_TPS_STABLE_TORCH_ABI OFF CACHE BOOL "Build openmm_cuda_bridge with LibTorch Stable ABI" FORCE)
        endif()
    else()
        message(STATUS "openmm_cuda_bridge: Stable ABI shim header missing; using legacy extension ABI")
        set(GENAI_TPS_STABLE_TORCH_ABI OFF CACHE BOOL "Build openmm_cuda_bridge with LibTorch Stable ABI" FORCE)
    endif()
endif()

find_package(CUDAToolkit REQUIRED)

set(OPENMM_LM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH "Path to the OpenMM source root")
set(OPENMM_PREFIX "" CACHE PATH "Installed OpenMM prefix (pixi/conda); used when OPENMM_BUILD_DIR unset")
set(OPENMM_BUILD_DIR "" CACHE PATH "OpenMM CMake build directory (must contain libOpenMM.so and libOpenMMCUDA.so)")
set(OPENMM_INCLUDE_DIR "${OPENMM_LM_ROOT}/openmmapi/include" CACHE PATH "OpenMM API include directory")
set(OPENMM_CUDA_INCLUDE_DIR "${OPENMM_LM_ROOT}/platforms/cuda/include" CACHE PATH "OpenMM CUDA include directory")
set(OPENMM_COMMON_INCLUDE_DIR "${OPENMM_LM_ROOT}/platforms/common/include" CACHE PATH "OpenMM common include directory")

if(NOT OPENMM_BUILD_DIR)
    if(OPENMM_PREFIX)
        set(_openmm_lib_dir "${OPENMM_PREFIX}/lib")
    elseif(DEFINED ENV{CONDA_PREFIX})
        set(_openmm_lib_dir "$ENV{CONDA_PREFIX}/lib")
    else()
        set(_openmm_lib_dir "${CMAKE_INSTALL_PREFIX}/lib")
    endif()
    if(EXISTS "${_openmm_lib_dir}/libOpenMM.so")
        set(OPENMM_BUILD_DIR "${_openmm_lib_dir}")
        message(STATUS "openmm_cuda_bridge: resolved OPENMM_BUILD_DIR=${OPENMM_BUILD_DIR}")
    endif()
endif()

if(NOT OPENMM_BUILD_DIR)
    message(FATAL_ERROR "OPENMM_BUILD_DIR must be set to a directory containing libOpenMM.so (or set OPENMM_PREFIX).")
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

set(BRIDGE_SOURCES bridge.cpp bridge.cu)
if(GENAI_TPS_STABLE_TORCH_ABI)
    list(APPEND BRIDGE_SOURCES bridge_stable_ops.cpp)
    message(STATUS "openmm_cuda_bridge: LibTorch Stable ABI enabled (min torch 2.10)")
endif()

pybind11_add_module(openmm_cuda_bridge MODULE ${BRIDGE_SOURCES})

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

if(GENAI_TPS_STABLE_TORCH_ABI)
    target_compile_definitions(openmm_cuda_bridge PRIVATE GENAI_TPS_STABLE_TORCH_ABI)
    set(_STABLE_TORCH_FLAGS "-DTORCH_STABLE_ONLY -DTORCH_TARGET_VERSION=0x020a000000000000")
    set_source_files_properties(bridge_stable_ops.cpp PROPERTIES COMPILE_FLAGS "${_STABLE_TORCH_FLAGS}")
else()
    target_compile_definitions(openmm_cuda_bridge PRIVATE GENAI_TPS_LEGACY_TORCH_ABI)
endif()

target_compile_definitions(openmm_cuda_bridge PRIVATE TORCH_API_INCLUDE_EXTENSION_H)
target_compile_options(openmm_cuda_bridge PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)

set(_openmm_rpath "${OPENMM_BUILD_DIR}")
set_target_properties(
    openmm_cuda_bridge
    PROPERTIES
        CUDA_SEPARABLE_COMPILATION OFF
        BUILD_RPATH_USE_ORIGIN ON
        BUILD_RPATH "${_openmm_rpath}"
        INSTALL_RPATH "${_openmm_rpath};$ORIGIN"
)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import site; print(site.getsitepackages()[0])"
    OUTPUT_VARIABLE GENAI_TPS_PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

install(
    TARGETS openmm_cuda_bridge
    LIBRARY DESTINATION "${GENAI_TPS_PYTHON_SITE_PACKAGES}"
    RUNTIME DESTINATION "${GENAI_TPS_PYTHON_SITE_PACKAGES}"
)
