#----------------------------------------------------
# OpenMM Brook Platform
#
# Creates OpenMM library, base name=OpenMM_Brook.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   OpenMM_Brook[_d].dll
#   OpenMM_Brook[_d].lib
#   OpenMM_BROOK_static[_d].lib
# Unix:
#   libOpenMM_Brook[_d].so
#   libOpenMM_BROOK_static[_d].a
#----------------------------------------------------

SUBDIRS (tests)

# logging

SET(LOG TRUE)

IF(LOG)
   SET(LOG_FILE "CMakeLog.txt" )
   FILE( WRITE ${LOG_FILE} "In Brook Cmake\n")
ENDIF(LOG)

# ----------------------------------------------------------------------------

IF(LOG)
   MACRO(LOG_DIR LOG_FILE DIR_LIST DIRS1 )
   
      FILE( APPEND ${LOG_FILE} "\n${DIR_LIST}\n")
      FOREACH(currentFile ${ARGN})
         FILE( APPEND ${LOG_FILE} "   ${currentFile}\n" )
      ENDFOREACH(currentFile)
   
   ENDMACRO(LOG_DIR) 
   
   LOG_DIR( ${LOG_FILE} "API_INCLUDE_DIRS" ${API_INCLUDE_DIRS} )
ENDIF(LOG)

# ----------------------------------------------------------------------------

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)

# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OpenMM_BROOK_LIBRARY_NAME OpenMM_Brook)

SET(SHARED_BROOK_TARGET ${OpenMM_BROOK_LIBRARY_NAME})
SET(STATIC_BROOK_TARGET ${OpenMM_BROOK_LIBRARY_NAME}_static)

# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.

IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_BROOK_TARGET ${SHARED_BROOK_TARGET}_d)
    SET(STATIC_BROOK_TARGET ${STATIC_BROOK_TARGET}_d)
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

# ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "API_ABS_INCLUDE_FILES" ${API_ABS_INCLUDE_FILES} )
   LOG_DIR( ${LOG_FILE} "API_REL_INCLUDE_FILES" ${API_REL_INCLUDE_FILES} )
   LOG_DIR( ${LOG_FILE} "CMAKE_CURRENT_SOURCE_DIR" ${CMAKE_CURRENT_SOURCE_DIR} )
ENDIF(LOG)
# ----------------------------------------------------------------------------

# collect cpp source files

SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)
# SET( CMAKE_CURRENT_SOURCE_DIR /home/friedrim/src/openmm/trunk/OpenMM/platforms/brook )
FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# ----------------------------------------------------------------------------
IF(LOG)
  LOG_DIR( ${LOG_FILE} "BROOK_SOURCE_FILES" ${SOURCE_FILES} )
  LOG_DIR( ${LOG_FILE} "SOURCE_INCLUDE_FILES" ${SOURCE_INCLUDE_FILES} )
ENDIF(LOG)
# ----------------------------------------------------------------------------

# Brook setup
 
# MESSAGE("CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}" )
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/brook-cmake/FindBrook.cmake)

INCLUDE_DIRECTORIES(${BROOK_INCLUDE_DIR})
LINK_DIRECTORIES(${${BROOK_brook_LIBRARY}})

# get *br files

FILE(GLOB BROOK_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/*.br)
FILE(GLOB BROOK_INCLUDE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/*.h)
 
FILE( APPEND ${LOG_FILE} "BROOK_SRC_FILES=${BROOK_SRC_FILES}\n" )

# ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "Brook src" ${BROOK_SRC_FILES} )
   LOG_DIR( ${LOG_FILE} "Brook include" ${BROOK_INCLUDE_FILES} )
ENDIF(LOG)
# ----------------------------------------------------------------------------

# create Brook custom rules

SET(BROOK_CPP_FILES)
FOREACH(brookFile ${BROOK_SRC_FILES})
  BROOK_FILE( ${brookFile} )
ENDFOREACH(brookFile)

# ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "Brook cpp" ${BROOK_CPP_FILES} )
ENDIF(LOG)
# ----------------------------------------------------------------------------
  
# BROOK_INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  
ADD_LIBRARY(${SHARED_BROOK_TARGET} SHARED ${BROOK_CPP_FILES} ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} )
ADD_LIBRARY(${STATIC_BROOK_TARGET} STATIC ${BROOK_CPP_FILES} ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} )

# ----------------------------------------------------------------------------
IF(LOG)
   FILE( APPEND ${LOG_FILE} "\nSTATIC_BROOK_TARGET=${STATIC_BROOK_TARGET} OPENMM_LIBRARY_NAME=${OPENMM_LIBRARY_NAME}\n" )
   FILE( APPEND ${LOG_FILE} "SHARED_BROOK_TARGET=${SHARED_BROOK_TARGET} SHARED_TARGET=${SHARED_TARGET} STATIC_TARGET=${STATIC_TARGET}\n" )
   FILE( APPEND ${LOG_FILE} "PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR} BROOK_LIB_PATH=${BROOK_LIB_PATH}\n" )
ENDIF(LOG)
# ----------------------------------------------------------------------------

TARGET_LINK_LIBRARIES(${SHARED_BROOK_TARGET} debug ${OPENMM_LIBRARY_NAME}_d optimized ${OPENMM_LIBRARY_NAME})
TARGET_LINK_LIBRARIES(${SHARED_BROOK_TARGET} debug brook_d optimized brook )
TARGET_LINK_LIBRARIES(${SHARED_BROOK_TARGET} ${SHARED_TARGET} )

LINK_DIRECTORIES(${SHARED_BROOK_TARGET} ${BROOK_LIB_PATH})
LINK_DIRECTORIES(${SHARED_BROOK_TARGET} ${PROJECT_BINARY_DIR})

TARGET_LINK_LIBRARIES(${STATIC_BROOK_TARGET} debug ${OPENMM_LIBRARY_NAME}_static_d optimized ${OPENMM_LIBRARY_NAME}_static)
TARGET_LINK_LIBRARIES(${STATIC_BROOK_TARGET} debug brook_d optimized brook)
TARGET_LINK_LIBRARIES(${STATIC_BROOK_TARGET} ${STATIC_TARGET})
LINK_DIRECTORIES(${STATIC_BROOK_TARGET} ${BROOK_LIB_PATH})

#---------------------------------------------------
