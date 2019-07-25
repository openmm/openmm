#---------------------------------------------------
# OpenMM CUDA Drude Integrator
#
# Creates OpenMMDrudeCUDA library.
#
# Windows:
#   OpenMMDrudeCUDA.dll
#   OpenMMDrudeCUDA.lib
# Unix:
#   libOpenMMDrudeCUDA.so
#----------------------------------------------------

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)


# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMDRUDECUDA_LIBRARY_NAME OpenMMDrudeCUDA)

SET(SHARED_TARGET ${OPENMMDRUDECUDA_LIBRARY_NAME})

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
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/cuda/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/cuda/src)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_BINARY_DIR}/platforms/cuda/src)

# Set variables needed for encoding kernel sources into a C++ class

SET(CUDA_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(CUDA_SOURCE_CLASS CudaDrudeKernelSources)
SET(CUDA_KERNELS_CPP ${CMAKE_CURRENT_BINARY_DIR}/src/${CUDA_SOURCE_CLASS}.cpp)
SET(CUDA_KERNELS_H ${CMAKE_CURRENT_BINARY_DIR}/src/${CUDA_SOURCE_CLASS}.h)
SET(SOURCE_FILES ${SOURCE_FILES} ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/src)

# Create the library

INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})

FILE(GLOB CUDA_KERNELS ${CUDA_SOURCE_DIR}/kernels/*.cu)
ADD_CUSTOM_COMMAND(OUTPUT ${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H}
    COMMAND ${CMAKE_COMMAND}
    ARGS -D CUDA_SOURCE_DIR=${CUDA_SOURCE_DIR} -D CUDA_KERNELS_CPP=${CUDA_KERNELS_CPP} -D CUDA_KERNELS_H=${CUDA_KERNELS_H} -D CUDA_SOURCE_CLASS=${CUDA_SOURCE_CLASS} -P ${CMAKE_SOURCE_DIR}/platforms/cuda/EncodeCUDAFiles.cmake
    DEPENDS ${CUDA_KERNELS}
)
SET_SOURCE_FILES_PROPERTIES(${CUDA_KERNELS_CPP} ${CUDA_KERNELS_H} PROPERTIES GENERATED TRUE)
ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${OPENMM_LIBRARY_NAME} ${PTHREADS_LIB})
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${OPENMM_LIBRARY_NAME}CUDA)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${SHARED_DRUDE_TARGET})
SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DOPENMM_BUILDING_SHARED_LIBRARY")
IF (APPLE)
    SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_COMPILE_FLAGS} -F/Library/Frameworks -framework CUDA")
ELSE (APPLE)
    SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}")
ENDIF (APPLE)

INSTALL(TARGETS ${SHARED_TARGET} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/plugins)
# Ensure that links to the main CUDA library will be resolved.
IF (APPLE)
    SET(CUDA_LIBRARY libOpenMMCUDA.dylib)
    INSTALL(CODE "EXECUTE_PROCESS(COMMAND install_name_tool -change ${CUDA_LIBRARY} @loader_path/${CUDA_LIBRARY} ${CMAKE_INSTALL_PREFIX}/lib/plugins/lib${SHARED_TARGET}.dylib)")
ENDIF (APPLE)
if(OPENMM_BUILD_CUDA_TESTS)
    SUBDIRS (tests)
endif(OPENMM_BUILD_CUDA_TESTS)
