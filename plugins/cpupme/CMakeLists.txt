#---------------------------------------------------
# OpenMM CPU PME Plugin
#
# Creates plugin library, base name=OpenMMPME.
# Default libraries are shared & optimized.
#
# Windows:
#   OpenMMPME[_d].dll
#   OpenMMPME[_d].lib
# Unix:
#   libOpenMMPME[_d].so
#----------------------------------------------------

IF (APPLE)
    SET (CMAKE_OSX_DEPLOYMENT_TARGET "10.6")
ENDIF (APPLE)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)


# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMPME_LIBRARY_NAME OpenMMPME)

SET(SHARED_TARGET ${OPENMMPME_LIBRARY_NAME})


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

# Find the include files.
SET(API_INCLUDE_FILES)
FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_INCLUDE_FILES ${API_INCLUDE_FILES} ${fullpaths})
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
SET_SOURCE_FILES_PROPERTIES(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS "-msse4.1")


# Include FFTW related files.
INCLUDE(FindFFTW)
INCLUDE_DIRECTORIES(${FFTW_INCLUDES})

# Build the plugin library.
ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})

IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(MAIN_OPENMM_LIB ${OPENMM_LIBRARY_NAME}_d)
ELSE (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(MAIN_OPENMM_LIB ${OPENMM_LIBRARY_NAME})
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${MAIN_OPENMM_LIB} ${PTHREADS_LIB} ${FFTW_LIBRARY} ${FFTW_THREADS_LIBRARY})
SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES COMPILE_FLAGS "-DOPENMM_PME_BUILDING_SHARED_LIBRARY")

INSTALL_TARGETS(/lib/plugins RUNTIME_DIRECTORY /lib/plugins ${SHARED_TARGET})

SUBDIRS(tests)
