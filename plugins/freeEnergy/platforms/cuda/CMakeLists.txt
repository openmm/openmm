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

# ----------------------------------------------------------------------------
# logging

SET(LOG TRUE)
   
IF(LOG)
   SET(LOG_FILE "CMakeLog.txt" )
   FILE( WRITE ${LOG_FILE} "In plugins/freeEnergy/platforms/cuda\n")
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
IF( (NOT DEFINED CUDA_HAVE_GPU) OR (CUDA_HAVE_GPU) )
    SUBDIRS (tests)
ENDIF( (NOT DEFINED CUDA_HAVE_GPU) OR (CUDA_HAVE_GPU) )

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

## ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "OPENMM_SOURCE_SUBDIRS"              ${OPENMM_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "OPENMM_FREE_ENERGY_SOURCE_SUBDIRS"  ${OPENMM_FREE_ENERGY_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "API_FREE_ENERGY_INCLUDE_DIRS"       ${API_FREE_ENERGY_INCLUDE_DIRS} )
   LOG_DIR( ${LOG_FILE} "CMAKE_CURRENT_SOURCE_DIR"           ${CMAKE_CURRENT_SOURCE_DIR} )
ENDIF(LOG)
## ----------------------------------------------------------------------------

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

## ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "API_FREE_ENERGY_REL_INCLUDE_FILES" ${API_FREE_ENERGY_REL_INCLUDE_FILES} )
   LOG_DIR( ${LOG_FILE} "OPENMM_DIR" ${OPENMM_DIR} )
ENDIF(LOG)
## ----------------------------------------------------------------------------

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_FREE_ENERGY_SOURCE_SUBDIRS})
    FILE(GLOB_RECURSE src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.c)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    LOG_DIR( ${LOG_FILE} "Adding include dir: " ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include )
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

SET(FINDCUDA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cuda-cmake)

SUBDIRS (sharedTarget staticTarget)
