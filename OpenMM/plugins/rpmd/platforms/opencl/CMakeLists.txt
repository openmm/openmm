#---------------------------------------------------
# OpenMM OpenCL RPMD Integrator
#
# Creates OpenMMRPMDOpenCL library.
#
# Windows:
#   OpenMMRPMDOpenCL.dll
#   OpenMMRPMDOpenCL.lib
# Unix:
#   libOpenMMRPMDOpenCL.so
#----------------------------------------------------

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)


# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMRPMDOPENCL_LIBRARY_NAME OpenMMRPMDOpenCL)

SET(SHARED_TARGET ${OPENMMRPMDOPENCL_LIBRARY_NAME})

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
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/opencl/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/opencl/src)

# Set variables needed for encoding kernel sources into a C++ class

SET(CL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
SET(CL_SOURCE_CLASS OpenCLRpmdKernelSources)
SET(CL_KERNELS_CPP ${CMAKE_CURRENT_BINARY_DIR}/src/${CL_SOURCE_CLASS}.cpp)
SET(CL_KERNELS_H ${CMAKE_CURRENT_BINARY_DIR}/src/${CL_SOURCE_CLASS}.h)
SET(SOURCE_FILES ${SOURCE_FILES} ${CL_KERNELS_CPP} ${CL_KERNELS_H})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR}/src)

# Create the library

INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})

FILE(GLOB OPENCL_KERNELS ${CL_SOURCE_DIR}/kernels/*.cl)
ADD_CUSTOM_COMMAND(OUTPUT ${CL_KERNELS_CPP} ${CL_KERNELS_H}
    COMMAND ${CMAKE_COMMAND}
    ARGS -D CL_SOURCE_DIR=${CL_SOURCE_DIR} -D CL_KERNELS_CPP=${CL_KERNELS_CPP} -D CL_KERNELS_H=${CL_KERNELS_H} -D CL_SOURCE_CLASS=${CL_SOURCE_CLASS} -P ${CMAKE_SOURCE_DIR}/platforms/opencl/EncodeCLFiles.cmake
    DEPENDS ${OPENCL_KERNELS}
)
SET_SOURCE_FILES_PROPERTIES(${CL_KERNELS_CPP} ${CL_KERNELS_H} PROPERTIES GENERATED TRUE)
ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${OPENMM_LIBRARY_NAME}  ${OPENCL_LIBRARIES} ${PTHREADS_LIB})
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${OPENMM_LIBRARY_NAME}OpenCL)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${SHARED_RPMD_TARGET})
SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DOPENMM_BUILDING_SHARED_LIBRARY")

INSTALL(TARGETS ${SHARED_TARGET} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/plugins)
# Ensure that links to the main OpenCL library will be resolved.
IF (APPLE)
    SET(OPENCL_LIBRARY libOpenMMOpenCL.dylib)
    INSTALL(CODE "EXECUTE_PROCESS(COMMAND install_name_tool -change ${OPENCL_LIBRARY} @loader_path/${OPENCL_LIBRARY} ${CMAKE_INSTALL_PREFIX}/lib/plugins/lib${SHARED_TARGET}.dylib)")
ENDIF (APPLE)
if(OPENMM_BUILD_OPENCL_TESTS)
    SUBDIRS (tests)
endif(OPENMM_BUILD_OPENCL_TESTS)
