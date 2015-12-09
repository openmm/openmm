#---------------------------------------------------
# OpenMM OpenCL Platform
#
# Creates OpenMMOpenCL library.
#
# Windows:
#   OpenMMOpenCL.dll
#   OpenMMOpenCL.lib
#   OpenMMOpenCL_static.lib
# Unix:
#   libOpenMMOpenCL.so
#   libOpenMMOpenCL_static.a
#----------------------------------------------------

set(OPENMM_BUILD_OPENCL_TESTS TRUE CACHE BOOL "Whether to build OpenCL test cases")
if(BUILD_TESTING AND OPENMM_BUILD_OPENCL_TESTS)
    SUBDIRS (tests)
endif(BUILD_TESTING AND OPENMM_BUILD_OPENCL_TESTS)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)


# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMOPENCL_LIBRARY_NAME OpenMMOpenCL)

SET(SHARED_TARGET ${OPENMMOPENCL_LIBRARY_NAME})
SET(STATIC_TARGET ${OPENMMOPENCL_LIBRARY_NAME}_static)


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
    IF(MSVC)
         FILE(GLOB_RECURSE kernel_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/kernels/*.cl)
         SET(SOURCE_FILES ${SOURCE_FILES} ${kernel_files})
    ENDIF(MSVC)
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Set variables needed for encoding kernel sources into a C++ class

SET(CL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(CL_SOURCE_CLASS OpenCLKernelSources)
SET(CL_KERNELS_CPP ${CMAKE_CURRENT_BINARY_DIR}/src/${CL_SOURCE_CLASS}.cpp)
SET(CL_KERNELS_H ${CMAKE_CURRENT_BINARY_DIR}/src/${CL_SOURCE_CLASS}.h)
SET(SOURCE_FILES ${SOURCE_FILES} ${CL_KERNELS_CPP} ${CL_KERNELS_H})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/src)

# Install headers

FILE(GLOB CORE_HEADERS include/*.h)
INSTALL_FILES(/include/openmm/opencl FILES ${CORE_HEADERS})

SUBDIRS (sharedTarget)
IF(OPENMM_BUILD_STATIC_LIB)
    SUBDIRS (staticTarget)
ENDIF(OPENMM_BUILD_STATIC_LIB)
