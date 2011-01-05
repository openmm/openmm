#---------------------------------------------------
# OpenMMFreeEnergy CUDA Platform
#
# Creates OpenMM library, base name=OpenMMFreeEnergyCuda.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   OpenMMFreeEnergyCuda[_d].dll
#   OpenMMFreeEnergyCuda[_d].lib
#   OpenMMFreeEnergyCuda_static[_d].lib
# Unix:
#   libOpenMMFreeEnergyCuda[_d].so
#   libOpenMMFreeEnergyCuda_static[_d].a
#----------------------------------------------------

# Only run tests when cuda-enabled GPU is found
SET(DO_TESTS FALSE)
IF(NOT DEFINED CUDA_HAVE_GPU)
    SET(DO_TESTS TRUE)
ENDIF(NOT DEFINED CUDA_HAVE_GPU)
IF(CUDA_HAVE_GPU)
    SET(DO_TESTS TRUE)
ENDIF(CUDA_HAVE_GPU)
IF(DO_TESTS)
    SUBDIRS (tests)
ENDIF(DO_TESTS)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_FREE_ENERGY_SOURCE_SUBDIRS .)

# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMCUDA_LIBRARY_NAME OpenMMFreeEnergyCuda)

SET(SHARED_TARGET ${OPENMMCUDA_LIBRARY_NAME})
SET(STATIC_TARGET ${OPENMMCUDA_LIBRARY_NAME}_static)

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
SET(API_FREE_ENERGY_INCLUDE_DIRS) # start empty
FOREACH(subdir ${OPENMM_FREE_ENERGY_SOURCE_SUBDIRS})
    # append
    SET(API_FREE_ENERGY_INCLUDE_DIRS ${API_FREE_ENERGY_INCLUDE_DIRS}
                                     ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/internal)
ENDFOREACH(subdir)

# We'll need both *relative* path names, starting with their API_FREE_ENERGY_INCLUDE_DIRS,
# and absolute pathnames.
SET(API_FREE_ENERGY_REL_INCLUDE_FILES)   # start these out empty
SET(API_FREE_ENERGY_ABS_INCLUDE_FILES)

FOREACH(dir ${API_FREE_ENERGY_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_FREE_ENERGY_ABS_INCLUDE_FILES ${API_FREE_ENERGY_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_FREE_ENERGY_REL_INCLUDE_FILES ${API_FREE_ENERGY_REL_INCLUDE_FILES} ${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_FREE_ENERGY_SOURCE_SUBDIRS})
    FILE(GLOB_RECURSE src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.c)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)
INCLUDE_DIRECTORIES(${OPENMM_DIR}/platforms/cuda/include)
INCLUDE_DIRECTORIES(${OPENMM_DIR}/platforms/cuda/src)
INCLUDE_DIRECTORIES(${OPENMM_DIR}/platforms/cuda/src/kernels)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)

# SET(FINDCUDA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cuda-cmake)

# SUBDIRS (sharedTarget staticTarget)
SUBDIRS (sharedTarget)

