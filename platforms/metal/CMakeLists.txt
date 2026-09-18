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

# Metal Platform: minimal implementation of the Common compute interfaces.
# This target does not yet register an OpenMM simulation Platform.

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

find_library(METAL_FRAMEWORK Metal)
find_library(FOUNDATION_FRAMEWORK Foundation)
if(NOT METAL_FRAMEWORK OR NOT FOUNDATION_FRAMEWORK)
    message(FATAL_ERROR "The Metal Platform requires the Metal and Foundation frameworks")
endif()

set(OPENMMMETAL_LIBRARY_NAME OpenMMMetal)
add_library(${OPENMMMETAL_LIBRARY_NAME} SHARED
    src/MetalContext.mm
    src/MetalArray.mm
    src/MetalKernel.mm
    src/MetalProgram.mm
    src/MetalQueue.mm
    src/MetalEvent.mm
    ../common/src/ComputeArray.cpp
    ../common/src/ComputeContext.cpp
    ../common/src/ComputeForceInfo.cpp)
target_include_directories(${OPENMMMETAL_LIBRARY_NAME}
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
           "${CMAKE_CURRENT_SOURCE_DIR}/../common/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(${OPENMMMETAL_LIBRARY_NAME}
    PUBLIC ${OPENMM_LIBRARY_NAME}
    PRIVATE ${METAL_FRAMEWORK} ${FOUNDATION_FRAMEWORK})
target_compile_options(${OPENMMMETAL_LIBRARY_NAME} PRIVATE
    "$<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-arc>")
target_compile_definitions(${OPENMMMETAL_LIBRARY_NAME} PRIVATE
    OPENMM_COMMON_BUILDING_SHARED_LIBRARY)
set_target_properties(${OPENMMMETAL_LIBRARY_NAME} PROPERTIES
    CXX_STANDARD 11
    OBJCXX_STANDARD 11
    COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_LINK_FLAGS}")

if(BUILD_TESTING)
    add_executable(TestMetalComputeContext tests/TestMetalComputeContext.cpp)
    target_link_libraries(TestMetalComputeContext PRIVATE ${OPENMMMETAL_LIBRARY_NAME})
    target_compile_definitions(TestMetalComputeContext PRIVATE
        METAL_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/kernels")
    set_target_properties(TestMetalComputeContext PROPERTIES
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
        LINK_FLAGS "${EXTRA_LINK_FLAGS}")
    add_test(NAME TestMetalComputeContext COMMAND TestMetalComputeContext)
    set_tests_properties(TestMetalComputeContext PROPERTIES SKIP_RETURN_CODE 77)
endif()
