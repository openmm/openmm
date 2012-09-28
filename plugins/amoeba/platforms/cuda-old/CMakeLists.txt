#---------------------------------------------------
# OpenMMAmoeba CUDA Platform
#
# Creates OpenMM library, base name=OpenMMAmoebaCuda.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   OpenMMAmoebaCuda[_d].dll
#   OpenMMAmoebaCuda[_d].lib
#   OpenMMAmoebaCuda_static[_d].lib
# Unix:
#   libOpenMMAmoebaCuda[_d].so
#   libOpenMMAmoebaCuda_static[_d].a
#----------------------------------------------------

# ----------------------------------------------------------------------------
# logging

SET(LOG FALSE)
   
IF(LOG)
   SET(LOG_FILE "CMakeLog.txt" )
   FILE( WRITE ${LOG_FILE} "In plugins/amoeba/platforms/cuda\n")
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
SET(OPENMM_AMOEBA_SOURCE_SUBDIRS .)

# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMCUDA_LIBRARY_NAME OpenMMAmoebaCuda)

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
SET(API_AMOEBA_INCLUDE_DIRS) # start empty
FOREACH(subdir ${OPENMM_AMOEBA_SOURCE_SUBDIRS})
    # append
    SET(API_AMOEBA_INCLUDE_DIRS ${API_AMOEBA_INCLUDE_DIRS}
                                     ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/internal)
ENDFOREACH(subdir)

## ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "OPENMM_SOURCE_SUBDIRS"              ${OPENMM_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "OPENMM_AMOEBA_SOURCE_SUBDIRS"       ${OPENMM_AMOEBA_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "API_AMOEBA_INCLUDE_DIRS"            ${API_AMOEBA_INCLUDE_DIRS} )
   LOG_DIR( ${LOG_FILE} "CMAKE_CURRENT_SOURCE_DIR"           ${CMAKE_CURRENT_SOURCE_DIR} )
ENDIF(LOG)
## ----------------------------------------------------------------------------

# We'll need both *relative* path names, starting with their API_AMOEBA_INCLUDE_DIRS,
# and absolute pathnames.
SET(API_AMOEBA_REL_INCLUDE_FILES)   # start these out empty
SET(API_AMOEBA_ABS_INCLUDE_FILES)

FOREACH(dir ${API_AMOEBA_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_AMOEBA_ABS_INCLUDE_FILES ${API_AMOEBA_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_AMOEBA_REL_INCLUDE_FILES ${API_AMOEBA_REL_INCLUDE_FILES} ${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

## ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "API_AMOEBA_REL_INCLUDE_FILES" ${API_AMOEBA_REL_INCLUDE_FILES} )
   LOG_DIR( ${LOG_FILE} "OPENMM_DIR" ${OPENMM_DIR} )
ENDIF(LOG)
## ----------------------------------------------------------------------------

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_AMOEBA_SOURCE_SUBDIRS})
    FILE(GLOB_RECURSE src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.c)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    IF(LOG)
        LOG_DIR( ${LOG_FILE} "Adding include dir: " ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include )
    ENDIF(LOG)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)
INCLUDE_DIRECTORIES(${OPENMM_DIR}/platforms/cuda/include)
INCLUDE_DIRECTORIES(${OPENMM_DIR}/platforms/cuda/src)
INCLUDE_DIRECTORIES(${OPENMM_DIR}/platforms/cuda/src/kernels)

IF(LOG)
   LOG_DIR( ${LOG_FILE} "CMAKE_CURRENT_SOURCE_DIR" ${CMAKE_CURRENT_SOURCE_DIR} )
   FILE( APPEND ${LOG_FILE} "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}\n" )
   LOG_DIR( ${LOG_FILE} "SOURCE_FILES" ${SOURCE_FILES} )
   LOG_DIR( ${LOG_FILE} "SOURCE_INCLUDE_FILES" ${SOURCE_INCLUDE_FILES} )
ENDIF(LOG)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)

# SET(FINDCUDA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cuda-cmake)

# SUBDIRS (sharedTarget staticTarget)
SUBDIRS (sharedTarget)

