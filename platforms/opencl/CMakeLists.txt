#---------------------------------------------------
# OpenMM OpenCL Platform
#
# Creates OpenMM library, base name=OpenMMOpenCL.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   OpenMMOpenCL[_d].dll
#   OpenMMOpenCL[_d].lib
#   OpenMMOpenCL_static[_d].lib
# Unix:
#   libOpenMMOpenCL[_d].so
#   libOpenMMOpenCL_static[_d].a
#----------------------------------------------------

IF (APPLE)
   SET (CMAKE_CXX_FLAGS "-arch i386")
   SET (CMAKE_C_FLAGS "-arch i386")
ENDIF (APPLE)

# Only run tests if this machine has a OpenCL-capable GPU
# So run a little test program at configuration time to sniff for GPUs
# But first, we need to set DYLD_LIBRARY_PATH here.
get_filename_component(OPENCL_LIB_DIR "${OPENCL_LIBRARY}" PATH)
file(TO_NATIVE_PATH "${OPENCL_LIB_DIR}" OPENCL_LIB_NATIVE_DIR)
if(APPLE)
        set(ENV{DYLD_LIBRARY_PATH} "$ENV{DYLD_LIBRARY_PATH}:${OPENCL_LIB_NATIVE_DIR}")
elseif(UNIX)
        set(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:${OPENCL_LIB_NATIVE_DIR}")
elseif(MSVC)
        set(ENV{PATH} "$ENV{PATH};${OPENCL_LIB_NATIVE_DIR}")
endif(APPLE)
try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
    ${CMAKE_BINARY_DIR} 
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/has_opencl_gpu.c
    CMAKE_FLAGS 
        -DINCLUDE_DIRECTORIES:STRING=${OPENCL_INCLUDE_DIR}
        -DLINK_LIBRARIES:STRING=${OPENCL_LIBRARY}
    COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT_VAR
    RUN_OUTPUT_VARIABLE RUN_OUTPUT_VAR
)
# message("${RUN_OUTPUT_VAR}") # Display number of GPUs found
# COMPILE_RESULT_VAR is TRUE when compile succeeds
# RUN_RESULT_VAR is zero when a GPU is found
set(OPENCL_HAVE_GPU_MAYBE TRUE)
if(NOT COMPILE_RESULT_VAR)
    set(OPENCL_HAVE_GPU_MAYBE FALSE)
endif(NOT COMPILE_RESULT_VAR)
if(RUN_RESULT_VAR)
    set(OPENCL_HAVE_GPU_MAYBE FALSE)
endif(RUN_RESULT_VAR)
if(OPENCL_HAVE_GPU_MAYBE)
    set(OPENCL_HAVE_GPU TRUE CACHE BOOL "Whether OpenCL-capable GPU is present")
    # Remain silent in case where GPU is found.
    # message("GPU check succeeded")
    # message("${COMPILE_OUTPUT_VAR}")
    # message("${RUN_OUTPUT_VAR}")
else(OPENCL_HAVE_GPU_MAYBE)
    set(OPENCL_HAVE_GPU FALSE CACHE BOOL "Whether OpenCL-capable GPU is present")
    message("OpenCL GPU check failed")
    message("${COMPILE_OUTPUT_VAR}")
    message("${RUN_OUTPUT_VAR}")
endif(OPENCL_HAVE_GPU_MAYBE)
if(OPENCL_HAVE_GPU)
    SUBDIRS (tests)
endif(OPENCL_HAVE_GPU)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)


# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMOPENCL_LIBRARY_NAME OpenMMOpenCL)

SET(SHARED_TARGET ${OPENMMOPENCL_LIBRARY_NAME})
SET(STATIC_TARGET ${OPENMMOPENCL_LIBRARY_NAME}_static)


# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.
IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_TARGET ${SHARED_TARGET}_d)
    SET(STATIC_TARGET ${STATIC_TARGET}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)


# These are all the places to search for header files which are
# to be part of the API.
SET(API_INCLUDE_DIRS) # start empty
FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include
                         ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/internal)
ENDFOREACH(subdir)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
SET(API_REL_INCLUDE_FILES)   # start these out empty
SET(API_ABS_INCLUDE_FILES)

FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES} ${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    FILE(GLOB_RECURSE src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.c)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Set variables needed for encoding kernel sources into a C++ class

SET(CL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(CL_KERNELS_CPP ${CMAKE_CURRENT_BINARY_DIR}/src/OpenCLKernelSources.cpp)
SET(CL_KERNELS_H ${CMAKE_CURRENT_BINARY_DIR}/src/OpenCLKernelSources.h)
SET(SOURCE_FILES ${SOURCE_FILES} ${CL_KERNELS_CPP} ${CL_KERNELS_H})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/src)

SUBDIRS (sharedTarget)
