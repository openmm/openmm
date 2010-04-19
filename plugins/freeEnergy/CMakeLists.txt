#---------------------------------------------------
# OpenMM FreeEnergy Plugin
#
# Creates OpenMM FreeEnergy plugin library, base name=OpenMMFreeEnergy.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   OpenMMFreeEnergy[_d].dll
#   OpenMMFreeEnergy[_d].lib
#   OpenMMFreeEnergy_static[_d].lib
# Unix:
#   libOpenMMFreeEnergy[_d].so
#   libOpenMMFreeEnergy_static[_d].a
#----------------------------------------------------

#INCLUDE(Dart)

SET(LOG TRUE)
# ----------------------------------------------------------------------------
IF(LOG)
   SET(LOG_FILE "CMakeLog.txt" )
   FILE( WRITE ${LOG_FILE} "In OpenMM Cmake\n")
ENDIF(LOG)

IF(LOG)
   MACRO(LOG_DIR LOG_FILE DIR_LIST )
    
      FILE( APPEND ${LOG_FILE} "\n${DIR_LIST}\n")
      FOREACH(currentFile ${ARGN})
         FILE( APPEND ${LOG_FILE} "   ${currentFile}\n" )
      ENDFOREACH(currentFile)
    
    ENDMACRO(LOG_DIR) 
ENDIF(LOG)

# ----------------------------------------------------------------------------

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_FREE_ENERGY_PLUGIN_SOURCE_SUBDIRS . openmmapi platforms/reference)

# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMM_FREE_ENERGY_LIBRARY_NAME OpenMMFreeEnergy)
SET(OPENMM_FREE_ENERGY_MAJOR_VERSION 1)
SET(OPENMM_FREE_ENERGY_MINOR_VERSION 0)
SET(OPENMM_FREE_ENERGY_BUILD_VERSION 0)

ADD_DEFINITIONS(-DOPENMM_FREE_ENERGY_LIBRARY_NAME=${OPENMM_FREE_ENERGY_LIBRARY_NAME}
                -DOPENMM_FREE_ENERGY_MAJOR_VERSION=${OPENMM_FREE_ENERGY_MAJOR_VERSION}
                -DOPENMM_FREE_ENERGY_MINOR_VERSION=${OPENMM_FREE_ENERGY_MINOR_VERSION}
		          -DOPENMM_FREE_ENERGY_BUILD_VERSION=${OPENMM_FREE_ENERGY_BUILD_VERSION})

SET(SHARED_FREE_ENERGY_TARGET ${OPENMM_FREE_ENERGY_LIBRARY_NAME})
SET(STATIC_FREE_ENERGY_TARGET ${OPENMM_FREE_ENERGY_LIBRARY_NAME}_static)

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_FREE_ENERGY_TARGET ${SHARED_FREE_ENERGY_TARGET}_d)
    SET(STATIC_FREE_ENERGY_TARGET ${STATIC_FREE_ENERGY_TARGET}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

# These are all the places to search for header files which are
# to be part of the API.
SET(API_FREE_ENERGY_INCLUDE_DIRS) # start empty
FOREACH(subdir ${OPENMM_FREE_ENERGY_PLUGIN_SOURCE_SUBDIRS})
    # append
    SET(API_FREE_ENERGY_INCLUDE_DIRS ${API_FREE_ENERGY_INCLUDE_DIRS}
                                     ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include 
                                     ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/internal)
ENDFOREACH(subdir)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
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
SET(SOURCE_FREE_ENERGY_FILES) # empty
SET(SOURCE_FREE_ENERGY_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_FREE_ENERGY_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FREE_ENERGY_FILES         ${SOURCE_FREE_ENERGY_FILES}         ${src_files})   #append
    SET(SOURCE_FREE_ENERGY_INCLUDE_FILES ${SOURCE_FREE_ENERGY_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)
INCLUDE_DIRECTORIES(BEFORE ${OPENMM_DIR}/platforms/reference/src)
INCLUDE_DIRECTORIES(BEFORE ${OPENMM_DIR}/platforms/reference/src/SimTKReference)

# ----------------------------------------------------------------------------
IF(LOG)
    FILE( APPEND ${LOG_FILE} "SOURCE_FREE_ENERGY_FILES=" ${SOURCE_FREE_ENERGY_FILES} "\n" )
    FILE( APPEND ${LOG_FILE} "SOURCE_FREE_ENERGY_INCLUDE_FILES=" ${SOURCE_FREE_ENERGY_INCLUDE_FILES} "\n" )
    FILE( APPEND ${LOG_FILE} "OPENMM_PATH="            ${OPENMM_PATH} "\n" )
    FILE( APPEND ${LOG_FILE} "CUDA_PATH="              ${CUDA_PATH} "\n" )
ENDIF(LOG)

# ----------------------------------------------------------------------------

# If API_FREE_ENERGY wrappers are being generated, and add them to the build.
#SET(OPENMM_BUILD_API_FREE_ENERGY_WRAPPERS OFF CACHE BOOL "Build wrappers for C and Fortran")
#IF(OPENMM_BUILD_API_FREE_ENERGY_WRAPPERS)
#    ADD_SUBDIRECTORY(wrappers)
#    SET(SOURCE_FREE_ENERGY_FILES ${SOURCE_FREE_ENERGY_FILES} wrappers/OpenMMCWrapper.cpp  wrappers/OpenMMFortranWrapper.cpp)
#    SET_SOURCE_FREE_ENERGY_FILES_PROPERTIES(wrappers/OpenMMCWrapper.cpp wrappers/OpenMMFortranWrapper.cpp PROPERTIES GENERATED TRUE)
#ENDIF(OPENMM_BUILD_API_FREE_ENERGY_WRAPPERS)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

ADD_LIBRARY(${SHARED_FREE_ENERGY_TARGET} SHARED ${SOURCE_FREE_ENERGY_FILES} ${SOURCE_FREE_ENERGY_INCLUDE_FILES} ${API_FREE_ENERGY_ABS_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_FREE_ENERGY_TARGET} PROPERTIES COMPILE_FLAGS "-DOPENMM_BUILDING_SHARED_LIBRARY -DLEPTON_BUILDING_SHARED_LIBRARY")

ADD_LIBRARY(${STATIC_FREE_ENERGY_TARGET} STATIC ${SOURCE_FREE_ENERGY_FILES} ${SOURCE_FREE_ENERGY_INCLUDE_FILES} ${API_FREE_ENERGY_ABS_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${STATIC_FREE_ENERGY_TARGET} PROPERTIES COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES -DOPENMM_BUILDING_STATIC_LIBRARY -DLEPTON_USE_STATIC_LIBRARIES -DLEPTON_BUILDING_STATIC_LIBRARY")

#IF(OPENMM_BUILD_API_FREE_ENERGY_WRAPPERS)
#    ADD_DEPENDENCIES(${SHARED_FREE_ENERGY_TARGET} ApiWrappers)
#    ADD_DEPENDENCIES(${STATIC_FREE_ENERGY_TARGET} ApiWrappers)
#ENDIF(OPENMM_BUILD_API_FREE_ENERGY_WRAPPERS)

# On Linux need to link to libdl
FIND_LIBRARY(DL_LIBRARY dl)
IF(DL_LIBRARY)
  TARGET_LINK_LIBRARIES(${SHARED_FREE_ENERGY_TARGET} ${DL_LIBRARY})
  TARGET_LINK_LIBRARIES(${STATIC_FREE_ENERGY_TARGET} ${DL_LIBRARY})
ENDIF(DL_LIBRARY)

SET( OpenMMLib OpenMM )
# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(OpenMMLib ${OpenMMLib}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
TARGET_LINK_LIBRARIES( ${SHARED_FREE_ENERGY_TARGET} ${SHARED_TARGET} )
TARGET_LINK_LIBRARIES( ${STATIC_FREE_ENERGY_TARGET} ${STATIC_TARGET} )

ADD_SUBDIRECTORY(platforms/reference/tests)

# Which hardware platforms to build
IF(CUDA_FOUND)
    SET(OPENMM_BUILD_FREE_ENERGY_CUDA_LIB ON CACHE BOOL "Build OpenMMFreeEnergyCuda library for Nvidia GPUs")
ELSE(CUDA_FOUND)
    SET(OPENMM_BUILD_FREE_ENERGY_CUDA_LIB OFF CACHE BOOL "Build OpenMMFreeEnergyCuda library for Nvidia GPUs")
ENDIF(CUDA_FOUND)

SET(OPENMM_BUILD_FREE_ENERGY_PATH)
SET(OPENMM_BUILD_FREE_ENERGY_CUDA_PATH)
IF(OPENMM_BUILD_FREE_ENERGY_CUDA_LIB)
    ADD_SUBDIRECTORY(platforms/cuda)
    SET(OPENMM_BUILD_FREE_ENERGY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/platforms/cuda)
    SET(OPENMM_BUILD_FREE_ENERGY_CUDA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/platforms/cuda)
    SET(OPENMM_FREE_ENERGY_CUDA_SOURCE_SUBDIRS . openmmapi olla platforms/cuda)
ENDIF(OPENMM_BUILD_FREE_ENERGY_CUDA_LIB)

#SET(OPENMM_BUILD_BROOK_LIB OFF CACHE BOOL "Build OpenMMBrook library for ATI GPUs")
#IF(OPENMM_BUILD_BROOK_LIB)
#    ADD_SUBDIRECTORY(platforms/brook)
#ENDIF(OPENMM_BUILD_BROOK_LIB)
#
#SET(OPENMM_BUILD_OPENCL_LIB OFF CACHE BOOL "Build OpenMMOpenCL library for Nvidia GPUs")
#IF(OPENMM_BUILD_OPENCL_LIB)
#    ADD_SUBDIRECTORY(platforms/opencl)
#ENDIF(OPENMM_BUILD_OPENCL_LIB)

INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_FREE_ENERGY_TARGET})
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${STATIC_FREE_ENERGY_TARGET})
FILE(GLOB CORE_HEADERS     include/*.h          */include/*.h)
FILE(GLOB TOP_HEADERS      include/openmm/*.h          */include/openmm/*.h)
FILE(GLOB INTERNAL_HEADERS include/openmm/internal/*.h */include/openmm/internal/*.h )
INSTALL_FILES(/include                 FILES ${CORE_HEADERS})
INSTALL_FILES(/include/openmm          FILES ${TOP_HEADERS})
INSTALL_FILES(/include/openmm/internal FILES ${INTERNAL_HEADERS})

#
# Allow automated build and dashboard.
#
#INCLUDE (Dart)

#IF (UNIX AND NOT CYGWIN AND NOT APPLE)
#  IF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
#    ADD_DEFINITIONS(-fprofile-arcs -ftest-coverage)
#    LINK_LIBRARIES(gcov)
#  ENDIF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
#ENDIF (UNIX AND NOT CYGWIN AND NOT APPLE)


#
# Testing
#

ENABLE_TESTING()

IF (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
ELSE (EXECUTABLE_OUTPUT_PATH)
  SET (TEST_PATH .)
ENDIF (EXECUTABLE_OUTPUT_PATH)

#INCLUDE(ApiDoxygen.cmake)

#ADD_SUBDIRECTORY(tests)
#ADD_SUBDIRECTORY(examples)
