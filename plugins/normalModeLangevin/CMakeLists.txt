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

set(NML_PLUGIN NormalModeLangevin)
# On Unix or Cygwin we have to add the debug suffix manually
if(UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    set(NML_PLUGIN ${NML_PLUGIN}_d)
endif(UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
add_library(${NML_PLUGIN} SHARED ${nml_sources} ${nml_headers})
set_target_properties(${NML_PLUGIN} PROPERTIES COMPILE_FLAGS "-DOPENMM_BUILDING_SHARED_LIBRARY")
target_link_libraries( ${NML_PLUGIN} ${SHARED_TARGET} )

if(BUILD_TESTING)
    # Copy plugin to test_plugin_dir
    set(test_plugin_dir "${CMAKE_BINARY_DIR}/test_plugin_dir")
    file(MAKE_DIRECTORY "${test_plugin_dir}")
    # On Windows we need to copy the correct Release/Debug plugin for testing
    if(MSVC)
        add_custom_command(TARGET ${NML_PLUGIN} POST_BUILD
            DEPENDS $(TargetPath)
            COMMAND "${CMAKE_COMMAND}"
            ARGS -E copy \"$\(TargetPath\)\" \"${test_plugin_dir}\"
            COMMENT "Copying normal mode langevin plugin for testing (WIN32)")
    else(MSVC)
        get_target_property(old_loc ${NML_PLUGIN} LOCATION)
        add_custom_command(TARGET ${NML_PLUGIN} POST_BUILD
            COMMAND "${CMAKE_COMMAND}"
            ARGS -E copy ${old_loc} ${test_plugin_dir}
            COMMENT "Copying normal mode langevin plugin for testing")
    endif(MSVC)
    add_subdirectory(test)
endif(BUILD_TESTING)

if(WIN32)
    # install DLL but not LIB
    install(TARGETS ${NML_PLUGIN} RUNTIME DESTINATION lib/plugins)
else(WIN32)
    install(TARGETS ${NML_PLUGIN} LIBRARY DESTINATION lib/plugins)
endif(WIN32)
