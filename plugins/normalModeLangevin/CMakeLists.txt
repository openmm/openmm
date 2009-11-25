#/* -------------------------------------------------------------------------- *
# *                                   OpenMM                                   *
# * -------------------------------------------------------------------------- *
# * This is part of the OpenMM molecular simulation toolkit originating from   *
# * Simbios, the NIH National Center for Physics-Based Simulation of           *
# * Biological Structures at Stanford, funded under the NIH Roadmap for        *
# * Medical Research, grant U54 GM072970. See https://simtk.org.               *
# *                                                                            *
# * Portions copyright (c) 2009 Stanford University and the Authors.           *
# * Authors: Chris Sweet, Christopher Bruns                                    *
# * Contributors:                                                              *
# *                                                                            *
# * Permission is hereby granted, free of charge, to any person obtaining a    *
# * copy of this software and associated documentation files (the "Software"), *
# * to deal in the Software without restriction, including without limitation  *
# * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
# * and/or sell copies of the Software, and to permit persons to whom the      *
# * Software is furnished to do so, subject to the following conditions:       *
# *                                                                            *
# * The above copyright notice and this permission notice shall be included in *
# * all copies or substantial portions of the Software.                        *
# *                                                                            *
# * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
# * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
# * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
# * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
# * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
# * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
# * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
# * -------------------------------------------------------------------------- */

# Create Normal Mode Langevin (NML) plugin for OpenMM

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
file(GLOB nml_headers "include/*.h")
file(GLOB nml_sources "src/*.cpp")

# Reference platform implementation
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/platforms/reference/include)
include_directories(${CMAKE_SOURCE_DIR}/platforms/reference/src)
file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/platforms/reference/include/*.h")
file(GLOB sources "${CMAKE_CURRENT_SOURCE_DIR}/platforms/reference/src/*.cpp")
set(nml_headers ${nml_headers} ${headers})
set(nml_sources ${nml_sources} ${sources})

add_library(NormalModeLangevin SHARED ${nml_sources} ${nml_headers})
set_target_properties(NormalModeLangevin PROPERTIES COMPILE_FLAGS "-DOPENMM_BUILDING_SHARED_LIBRARY")
# On Unix or Cygwin we have to add the debug suffix manually
if(UNIX)
    set_target_properties(NormalModeLangevin PROPERTIES DEBUG_OUTPUT_NAME NormalModeLangevin_d)
endif(UNIX)
target_link_libraries( NormalModeLangevin OpenMM )

if(BUILD_TESTING)
    add_subdirectory(test)
endif(BUILD_TESTING)

install(TARGETS NormalModeLangevin RUNTIME DESTINATION lib/plugins)

