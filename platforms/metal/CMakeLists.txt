# --------------------------------------------------------------------------
# OpenMM
#
# This is part of the OpenMM molecular simulation toolkit.
# See https://openmm.org/development.
#
# Metal Platform code:
# Portions copyright (c) 2026 Chun-Chi Hung.
# Authors: Chun-Chi Hung
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------------

# Metal Platform: Common Compute with an OpenCL-semantic GPU baseline.

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    message(FATAL_ERROR "The Metal Platform is only available on macOS")
endif()
if(NOT CMAKE_OSX_ARCHITECTURES STREQUAL "arm64")
    message(FATAL_ERROR "The Metal Platform currently requires an arm64-only build")
endif()
if(NOT OPENMM_BUILD_SHARED_LIB)
    message(FATAL_ERROR "The Metal Platform currently requires OPENMM_BUILD_SHARED_LIB=ON")
endif()

# Keep the minimum target local to Metal and preserve a newer user target.
if(CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS "13.0")
    set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0")
endif()
enable_language(OBJCXX)

option(OPENMM_METAL_RECORD_AND_COMMIT
    "Commit each Metal compute/blit operation instead of batching (no additional CPU wait)" ON)
option(OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
    "Use scoped GPU floating accumulators for minimization beyond the Q32.32 force range" ON)
option(OPENMM_METAL_NATIVE_FLOAT_ATOMICS
    "Use native Metal float atomic addition instead of the OpenCL-style CAS loop" OFF)
option(OPENMM_METAL_FAST_MATH
    "Allow broad Metal fast math except where compensated or robust arithmetic is required" OFF)
option(OPENMM_METAL_FAST_MINIMIZE_SHUFFLE
    "Use the Common CUDA/HIP shuffle reduction in the minimizer" OFF)
option(OPENMM_METAL_FAST_SHORT_LIST_SORT
    "Use the CUDA short-list sorting strategy for small arrays" OFF)
option(OPENMM_METAL_FAST_CONSTANT_POTENTIAL_REDUCTION
    "Use the Common CUDA/HIP shuffle reduction in the base ConstantPotential kernels" OFF)
option(OPENMM_METAL_FAST_CONSTANT_POTENTIAL_CG_REDUCTION
    "Use the Common CUDA/HIP shuffle reductions in the ConstantPotential CG solver" OFF)
option(OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_REDUCTION
    "Use the Common CUDA/HIP shuffle reduction in the ConstantPotential matrix solver" OFF)
option(OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_BROADCAST
    "Use the Common CUDA/HIP lane broadcasts in the ConstantPotential matrix solver" OFF)
option(OPENMM_METAL_FAST_CUSTOM_NONBONDED_GROUPS_SHUFFLE
    "Use the Common CUDA xor-shuffle maximum in CustomNonbonded interaction groups" OFF)
option(OPENMM_METAL_FAST_CUSTOM_MANY_PARTICLE_BALLOT
    "Use the Common CUDA neighbor-block ballot in CustomManyParticle" OFF)
option(OPENMM_METAL_FAST_LCPO_BALLOT
    "Use the Common CUDA neighbor-block ballot in LCPO" OFF)
option(OPENMM_METAL_FAST_PME_FLOAT_SPREAD "Use the Common floating-grid spreading path for PME" OFF)
option(OPENMM_METAL_FAST_LJPME_FLOAT_SPREAD "Use floating grids for both reciprocal terms of LJPME" OFF)
option(OPENMM_METAL_FAST_CONSTANT_POTENTIAL_FLOAT_SPREAD "Use floating-grid spreading in ConstantPotential PME" OFF)

# Each optional algorithm can be compared independently with the OpenCL path.
set(METAL_OPTIMIZATION_SWITCHES
    BLOCK_BOUNDS FP16_BOUNDS NEIGHBOR_BALLOT SPARSE_PAIRS NONBONDED_SHUFFLE CUSTOM_GB_VALUE_SHUFFLE CUSTOM_GB_ENERGY_SHUFFLE GBSA_BORN_SHUFFLE GBSA_FORCE_SHUFFLE DPD_PARTICLE_SHUFFLE DPD_TILE_BROADCAST CUSTOM_HBOND_SHUFFLE CENTROID_REDUCTION RG_REDUCTION RMSD_REDUCTION ORIENTATION_REDUCTION LCPO_NEIGHBOR_SCAN CUSTOM_MANY_PARTICLE_NEIGHBOR_SCAN SORT_BUCKET_SCAN)
foreach(feature IN LISTS METAL_OPTIMIZATION_SWITCHES)
    option(OPENMM_METAL_FAST_${feature} "Enable the independent Metal ${feature} algorithm" OFF)
endforeach()

find_library(METAL_FRAMEWORK Metal)
find_library(FOUNDATION_FRAMEWORK Foundation)
find_library(QUARTZCORE_FRAMEWORK QuartzCore)
if(NOT METAL_FRAMEWORK OR NOT FOUNDATION_FRAMEWORK)
    message(FATAL_ERROR "The Metal Platform requires the Metal and Foundation frameworks")
endif()

# VkFFT's Metal backend uses metal-cpp. Keep C++17 private to that translation
# unit; public interfaces and the Objective-C++ runtime remain C++11.
find_path(METAL_CPP_INCLUDE_DIR Metal/Metal.hpp DOC "Path to Apple's metal-cpp headers")
if(NOT METAL_CPP_INCLUDE_DIR OR NOT QUARTZCORE_FRAMEWORK)
    message(FATAL_ERROR "Metal VkFFT requires metal-cpp headers (METAL_CPP_INCLUDE_DIR) and QuartzCore")
endif()

# Keep this dependency check's C++17 requirement private to the probe and VkFFT.
function(check_metal_cpp_math_options)
    include(CheckCXXSourceCompiles)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
    set(CMAKE_REQUIRED_INCLUDES "${METAL_CPP_INCLUDE_DIR}")
    # A reconfigure may point at a different external header release.
    unset(OPENMM_METAL_CPP_HAS_MATH_OPTIONS CACHE)
    check_cxx_source_compiles([=[
        #include <Metal/Metal.hpp>
        void check(MTL::CompileOptions* options) {
            options->setMathMode(MTL::MathModeSafe);
            options->setMathMode(MTL::MathModeFast);
            options->setMathFloatingPointFunctions(MTL::MathFloatingPointFunctionsPrecise);
            options->setMathFloatingPointFunctions(MTL::MathFloatingPointFunctionsFast);
        }
    ]=] OPENMM_METAL_CPP_HAS_MATH_OPTIONS)
    if(NOT OPENMM_METAL_CPP_HAS_MATH_OPTIONS)
        message(FATAL_ERROR "Metal requires metal-cpp headers with mathMode and mathFloatingPointFunctions (macOS 15 release or newer)")
    endif()
endfunction()
check_metal_cpp_math_options()

function(encode_metal_sources class source_dir extension)
    if(ARGC GREATER 3)
        set(source_file "${source_dir}/${ARGV3}.${extension}")
        set(kernels "${source_file}")
    else()
        set(source_file "")
        file(GLOB kernels "${source_dir}/*.${extension}")
    endif()
    add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/src/${class}.cpp" "${CMAKE_CURRENT_BINARY_DIR}/src/${class}.h"
        COMMAND ${CMAKE_COMMAND} "-DSOURCE_DIR=${source_dir}" "-DEXTENSION=${extension}"
            "-DCLASS=${class}" "-DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/src" "-DSOURCE_FILE=${source_file}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/EncodeSources.cmake"
        DEPENDS ${kernels} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/EncodeSources.cmake")
endfunction()
encode_metal_sources(MetalKernelSources "${CMAKE_CURRENT_SOURCE_DIR}/src/kernels" metal)
encode_metal_sources(MetalOpenCLKernelSources "${CMAKE_CURRENT_SOURCE_DIR}/../opencl/src/kernels" cl)
encode_metal_sources(MetalCudaKernelSources "${CMAKE_CURRENT_SOURCE_DIR}/../cuda/src/kernels" cu nonbonded)
file(GLOB COMMON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/../common/src/*.cpp")
list(REMOVE_ITEM COMMON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/../common/src/CommonParallelKernels.cpp")

add_library(OpenMMMetalVkFFT OBJECT src/MetalVkFFT.cpp)
target_compile_definitions(OpenMMMetalVkFFT PRIVATE
    OPENMM_METAL_FAST_MATH=$<BOOL:${OPENMM_METAL_FAST_MATH}>)
target_include_directories(OpenMMMetalVkFFT SYSTEM PRIVATE
    "${METAL_CPP_INCLUDE_DIR}" "${CMAKE_SOURCE_DIR}/libraries/vkfft/include")
set_target_properties(OpenMMMetalVkFFT PROPERTIES CXX_STANDARD 17 POSITION_INDEPENDENT_CODE ON)

set(OPENMMMETAL_LIBRARY_NAME OpenMMMetal)
set_source_files_properties("${CMAKE_CURRENT_BINARY_DIR}/../common/src/CommonKernelSources.cpp" PROPERTIES GENERATED TRUE)
add_library(${OPENMMMETAL_LIBRARY_NAME} SHARED
    src/MetalContext.mm
    src/MetalArray.mm
    src/MetalKernel.mm
    src/MetalProgram.mm
    src/MetalQueue.mm
    src/MetalEvent.mm
    src/MetalSourceAdapter.cpp
    src/MetalPairwiseOptimizations.cpp
    src/MetalReductionOptimizations.cpp
    src/MetalPlatform.mm
    src/MetalKernelFactory.cpp
    src/MetalKernels.cpp
    src/MetalNonbondedUtilities.cpp
    src/MetalIntegrationUtilities.cpp
    src/MetalSort.cpp
    src/MetalFFT3D.mm
    $<TARGET_OBJECTS:OpenMMMetalVkFFT>
    ${COMMON_SOURCES}
    "${CMAKE_CURRENT_BINARY_DIR}/src/MetalKernelSources.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/src/MetalOpenCLKernelSources.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/src/MetalCudaKernelSources.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/../common/src/CommonKernelSources.cpp")
add_dependencies(${OPENMMMETAL_LIBRARY_NAME} CommonKernels)
target_include_directories(${OPENMMMETAL_LIBRARY_NAME}
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
           "${CMAKE_CURRENT_SOURCE_DIR}/../common/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_BINARY_DIR}/src"
            "${CMAKE_CURRENT_BINARY_DIR}/../common/src")
target_link_libraries(${OPENMMMETAL_LIBRARY_NAME}
    PUBLIC ${OPENMM_LIBRARY_NAME}
    PRIVATE ${METAL_FRAMEWORK} ${FOUNDATION_FRAMEWORK} ${QUARTZCORE_FRAMEWORK})
target_compile_options(${OPENMMMETAL_LIBRARY_NAME} PRIVATE
    "$<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-arc>")
target_compile_definitions(${OPENMMMETAL_LIBRARY_NAME} PRIVATE
    OPENMM_COMMON_BUILDING_SHARED_LIBRARY
    OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS=$<BOOL:${OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS}>
    OPENMM_METAL_NATIVE_FLOAT_ATOMICS=$<BOOL:${OPENMM_METAL_NATIVE_FLOAT_ATOMICS}>
    OPENMM_METAL_FAST_PME_FLOAT_SPREAD=$<BOOL:${OPENMM_METAL_FAST_PME_FLOAT_SPREAD}>
    OPENMM_METAL_FAST_LJPME_FLOAT_SPREAD=$<BOOL:${OPENMM_METAL_FAST_LJPME_FLOAT_SPREAD}>
    OPENMM_METAL_FAST_CONSTANT_POTENTIAL_FLOAT_SPREAD=$<BOOL:${OPENMM_METAL_FAST_CONSTANT_POTENTIAL_FLOAT_SPREAD}>

    OPENMM_METAL_FAST_MATH=$<BOOL:${OPENMM_METAL_FAST_MATH}>
    OPENMM_METAL_FAST_MINIMIZE_SHUFFLE=$<BOOL:${OPENMM_METAL_FAST_MINIMIZE_SHUFFLE}>
    OPENMM_METAL_FAST_SHORT_LIST_SORT=$<BOOL:${OPENMM_METAL_FAST_SHORT_LIST_SORT}>
    OPENMM_METAL_FAST_CONSTANT_POTENTIAL_REDUCTION=$<BOOL:${OPENMM_METAL_FAST_CONSTANT_POTENTIAL_REDUCTION}>
    OPENMM_METAL_FAST_CONSTANT_POTENTIAL_CG_REDUCTION=$<BOOL:${OPENMM_METAL_FAST_CONSTANT_POTENTIAL_CG_REDUCTION}>
    OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_REDUCTION=$<BOOL:${OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_REDUCTION}>
    OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_BROADCAST=$<BOOL:${OPENMM_METAL_FAST_CONSTANT_POTENTIAL_MATRIX_BROADCAST}>
    OPENMM_METAL_FAST_CUSTOM_NONBONDED_GROUPS_SHUFFLE=$<BOOL:${OPENMM_METAL_FAST_CUSTOM_NONBONDED_GROUPS_SHUFFLE}>
    OPENMM_METAL_FAST_CUSTOM_MANY_PARTICLE_BALLOT=$<BOOL:${OPENMM_METAL_FAST_CUSTOM_MANY_PARTICLE_BALLOT}>
    OPENMM_METAL_FAST_LCPO_BALLOT=$<BOOL:${OPENMM_METAL_FAST_LCPO_BALLOT}>
    OPENMM_METAL_RECORD_AND_COMMIT=$<BOOL:${OPENMM_METAL_RECORD_AND_COMMIT}>)
foreach(feature IN LISTS METAL_OPTIMIZATION_SWITCHES)
    target_compile_definitions(${OPENMMMETAL_LIBRARY_NAME} PRIVATE
        OPENMM_METAL_FAST_${feature}=$<BOOL:${OPENMM_METAL_FAST_${feature}}>)
endforeach()
set_target_properties(${OPENMMMETAL_LIBRARY_NAME} PROPERTIES
    CXX_STANDARD 11
    OBJCXX_STANDARD 11
    COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_LINK_FLAGS}")

# The core Platform loads through OpenMM's normal platform-library mechanism;
# this does not enable optional force plugins such as AMOEBA or Drude.
install(TARGETS ${OPENMMMETAL_LIBRARY_NAME} LIBRARY DESTINATION lib/plugins)
file(GLOB METAL_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
install(FILES ${METAL_HEADERS} DESTINATION include/openmm/metal)

if(BUILD_TESTING)
    add_executable(TestMetalPlatform tests/TestMetalPlatform.cpp)
    target_link_libraries(TestMetalPlatform PRIVATE ${OPENMM_LIBRARY_NAME})
    target_compile_definitions(TestMetalPlatform PRIVATE METAL_LIBRARY_PATH="$<TARGET_FILE:OpenMMMetal>")
    add_dependencies(TestMetalPlatform OpenMMMetal)
    add_test(NAME TestMetalPlatform COMMAND TestMetalPlatform)
    set_tests_properties(TestMetalPlatform PROPERTIES SKIP_RETURN_CODE 77 LABELS "Metal")
    foreach(test IN ITEMS CommonSources SourceAdapter Sort FFT3D Primitives PairwiseOptimizations ReductionOptimizations NonbondedOptimizations)
        add_executable(TestMetal${test} tests/TestMetal${test}.cpp)
        target_include_directories(TestMetal${test} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/../common/src"
            "${CMAKE_CURRENT_BINARY_DIR}/src")
        target_link_libraries(TestMetal${test} PRIVATE ${OPENMMMETAL_LIBRARY_NAME})
        add_test(NAME TestMetal${test} COMMAND TestMetal${test})
        set_tests_properties(TestMetal${test} PROPERTIES SKIP_RETURN_CODE 77)
    endforeach()
    target_include_directories(TestMetalNonbondedOptimizations PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    foreach(feature IN ITEMS BLOCK_BOUNDS FP16_BOUNDS NEIGHBOR_BALLOT SPARSE_PAIRS NONBONDED_SHUFFLE)
        target_compile_definitions(TestMetalNonbondedOptimizations PRIVATE
            OPENMM_METAL_FAST_${feature}=$<BOOL:${OPENMM_METAL_FAST_${feature}}>)
    endforeach()
    foreach(option IN ITEMS CUSTOM_NONBONDED_GROUPS_SHUFFLE CUSTOM_MANY_PARTICLE_BALLOT LCPO_BALLOT)
        target_compile_definitions(TestMetalSourceAdapter PRIVATE
            OPENMM_METAL_FAST_${option}=$<BOOL:${OPENMM_METAL_FAST_${option}}>)
    endforeach()
    # Reuse the existing single-precision core tests without copying algorithms.
    target_compile_definitions(TestMetalSourceAdapter PRIVATE
        OPENMM_METAL_NATIVE_FLOAT_ATOMICS=$<BOOL:${OPENMM_METAL_NATIVE_FLOAT_ATOMICS}>
        OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS=$<BOOL:${OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS}>)
    set(METAL_CORE_TESTS HarmonicBondForce HarmonicAngleForce PeriodicTorsionForce
        RBTorsionForce CMAPTorsionForce CustomBondForce CustomAngleForce
        CustomTorsionForce CustomExternalForce CustomCompoundBondForce
        CustomCentroidBondForce NonbondedForce Ewald DispersionPME
        GBSAOBCForce CustomGBForce CustomNonbondedForce CustomHbondForce
        CustomManyParticleForce GayBerneForce LCPOForce ConstantPotentialForce
        CustomCVForce ATMForce CustomCPPForce PythonForce OrientationRestraintForce RGForce
        RMSDForce VerletIntegrator LangevinIntegrator LangevinMiddleIntegrator
        BrownianIntegrator VariableVerletIntegrator VariableLangevinIntegrator
        CustomIntegrator CompoundIntegrator DPDIntegrator QTBIntegrator
        NoseHooverIntegrator AndersenThermostat CMMotionRemover VirtualSites
        Settle Checkpoints MultipleForces LocalEnergyMinimizer
        MonteCarloBarostat MonteCarloAnisotropicBarostat MonteCarloFlexibleBarostat)
    foreach(test IN LISTS METAL_CORE_TESTS)
        add_executable(TestMetal${test} tests/TestMetalCommon.cpp)
        target_link_libraries(TestMetal${test} PRIVATE ${OPENMMMETAL_LIBRARY_NAME})
        target_compile_definitions(TestMetal${test} PRIVATE METAL_COMMON_TEST_HEADER="Test${test}.h")
        add_test(NAME TestMetal${test} COMMAND TestMetal${test})
        set_tests_properties(TestMetal${test} PROPERTIES LABELS "Metal;Core;Single" TIMEOUT 600)
    endforeach()
    target_compile_definitions(TestMetalNonbondedForce PRIVATE METAL_TEST_NONBONDED)
    target_compile_definitions(TestMetalCheckpoints PRIVATE METAL_TEST_CHECKPOINTS)
    target_compile_definitions(TestMetalLocalEnergyMinimizer PRIVATE METAL_TEST_MINIMIZER)
    target_compile_definitions(TestMetalLocalEnergyMinimizer PRIVATE
        OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS=$<BOOL:${OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS}>)
    target_compile_definitions(TestMetalConstantPotentialForce PRIVATE METAL_TEST_CONSTANT_POTENTIAL)
    target_compile_definitions(TestMetalMonteCarloBarostat PRIVATE METAL_TEST_MONTE_CARLO_BAROSTAT)
    target_compile_definitions(TestMetalMonteCarloAnisotropicBarostat PRIVATE METAL_TEST_MONTE_CARLO_ANISOTROPIC_BAROSTAT)
    # Validate the private GPU-only minimization representation with unchanged
    # force assertions, not only the extreme-overlap minimizer regression.
    set(METAL_FLOAT_TESTS HarmonicBondForce HarmonicAngleForce PeriodicTorsionForce
        RBTorsionForce CMAPTorsionForce CustomBondForce CustomAngleForce
        CustomTorsionForce CustomExternalForce CustomCompoundBondForce
        NonbondedForce Ewald DispersionPME GBSAOBCForce CustomGBForce
        CustomNonbondedForce CustomHbondForce CustomManyParticleForce GayBerneForce
        LCPOForce ConstantPotentialForce CustomCVForce ATMForce CustomCPPForce
        PythonForce CustomCentroidBondForce OrientationRestraintForce RGForce RMSDForce
        MultipleForces VirtualSites LocalEnergyMinimizer)
    if(NOT OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS)
        # A permanently floating test Context is not the fixed-point minimizer baseline.
        list(REMOVE_ITEM METAL_FLOAT_TESTS LocalEnergyMinimizer)
    endif()
    foreach(test IN LISTS METAL_FLOAT_TESTS)
        add_executable(TestMetalFloat${test} tests/TestMetalCommon.cpp)
        target_link_libraries(TestMetalFloat${test} PRIVATE ${OPENMMMETAL_LIBRARY_NAME})
        target_compile_definitions(TestMetalFloat${test} PRIVATE
            METAL_TEST_FLOAT_ACCUMULATORS METAL_COMMON_TEST_HEADER="Test${test}.h")
        add_test(NAME TestMetalFloat${test} COMMAND TestMetalFloat${test})
        set_tests_properties(TestMetalFloat${test} PROPERTIES LABELS "Metal;FloatAccumulators;Single" TIMEOUT 600)
    endforeach()
    target_compile_definitions(TestMetalFloatNonbondedForce PRIVATE METAL_TEST_NONBONDED)
    target_compile_definitions(TestMetalFloatConstantPotentialForce PRIVATE METAL_TEST_CONSTANT_POTENTIAL)
    if(TARGET TestMetalFloatLocalEnergyMinimizer)
        target_compile_definitions(TestMetalFloatLocalEnergyMinimizer PRIVATE METAL_TEST_MINIMIZER
            OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS=$<BOOL:${OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS}>)
    endif()
    add_executable(TestMetalComputeContext tests/TestMetalComputeContext.cpp)
    target_link_libraries(TestMetalComputeContext PRIVATE ${OPENMMMETAL_LIBRARY_NAME})
    target_compile_definitions(TestMetalComputeContext PRIVATE
        METAL_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/kernels")
    set_target_properties(TestMetalComputeContext PROPERTIES
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
        LINK_FLAGS "${EXTRA_LINK_FLAGS}")
    add_test(NAME TestMetalComputeContext COMMAND TestMetalComputeContext)
    set_tests_properties(TestMetalComputeContext PROPERTIES SKIP_RETURN_CODE 77)

    # Exercise both policies independently of the selected platform build option.
    foreach(mode IN ITEMS Immediate Batched)
        if(mode STREQUAL "Immediate")
            set(record_and_commit 1)
        else()
            set(record_and_commit 0)
        endif()
        add_executable(TestMetalQueue${mode} tests/TestMetalQueue.mm src/MetalQueue.mm)
        target_include_directories(TestMetalQueue${mode} PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/../common/include")
        target_link_libraries(TestMetalQueue${mode} PRIVATE
            ${OPENMM_LIBRARY_NAME} ${METAL_FRAMEWORK} ${FOUNDATION_FRAMEWORK})
        target_compile_options(TestMetalQueue${mode} PRIVATE -fobjc-arc)
        target_compile_definitions(TestMetalQueue${mode} PRIVATE
            OPENMM_METAL_RECORD_AND_COMMIT=${record_and_commit})
        set_target_properties(TestMetalQueue${mode} PROPERTIES
            OBJCXX_STANDARD 11
            COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
            LINK_FLAGS "${EXTRA_LINK_FLAGS}")
        add_test(NAME TestMetalQueue${mode} COMMAND TestMetalQueue${mode})
        set_tests_properties(TestMetalQueue${mode} PROPERTIES SKIP_RETURN_CODE 77)
    endforeach()
endif()
