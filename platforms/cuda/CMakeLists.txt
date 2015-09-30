#---------------------------------------------------
# OpenMM CUDA Platform
#
# Creates OpenMMCUDA library.
#
# Windows:
#   OpenMMCUDA.dll
#   OpenMMCUDA.lib
#   OpenMMCUDA_static.lib
# Unix:
#   libOpenMMCUDA.so
#   libOpenMMCUDA_static.a
#----------------------------------------------------

set(OPENMM_BUILD_CUDA_TESTS TRUE CACHE BOOL "Whether to build CUDA test cases")
if(BUILD_TESTING AND OPENMM_BUILD_CUDA_TESTS)
    SUBDIRS (tests)
endif(BUILD_TESTING AND OPENMM_BUILD_CUDA_TESTS)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)


# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMCUDA_LIBRARY_NAME OpenMMCUDA)

SET(SHARED_TARGET ${OPENMMCUDA_LIBRARY_NAME})
SET(STATIC_TARGET ${OPENMMCUDA_LIBRARY_NAME}_static)


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
         FILE(GLOB_RECURSE kernel_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/kernels/*.cu)
         SET(SOURCE_FILES ${SOURCE_FILES} ${kernel_files})
    ENDIF(MSVC)
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Set variables needed for encoding kernel sources into a C++ class

SET(CUDA_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(CUDA_SOURCE_CLASS CudaKernelSources)
SET(CUDA_KERNELS_CPP ${CMAKE_CURRENT_BINARY_DIR}/src/${CUDA_SOURCE_CLASS}.cpp)
SET(CUDA_KERNELS_H ${CMAKE_CURRENT_BINARY_DIR}/src/${CUDA_SOURCE_CLASS}.h)
SET(SOURCE_FILES ${SOURCE_FILES} ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/src)

# Install headers

FILE(GLOB CORE_HEADERS include/*.h)
INSTALL_FILES(/include/openmm/cuda FILES ${CORE_HEADERS})

SUBDIRS (sharedTarget)
IF(OPENMM_BUILD_STATIC_LIB)
    SUBDIRS (staticTarget)
ENDIF(OPENMM_BUILD_STATIC_LIB)
