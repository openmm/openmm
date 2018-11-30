#---------------------------------------------------
# OpenMM CPU PME Plugin
#
# Creates OpenMMPME plugin library.
#
# Windows:
#   OpenMMPME.dll
#   OpenMMPME.lib
#   OpenMMPME_static.lib
# Unix:
#   libOpenMMPME.so
#   libOpenMMPME_static.a
#----------------------------------------------------

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)


# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMPME_LIBRARY_NAME OpenMMPME)

SET(SHARED_TARGET ${OPENMMPME_LIBRARY_NAME})
SET(STATIC_TARGET ${OPENMMPME_LIBRARY_NAME}_static)


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
IF(NOT MSVC)
    IF(X86)
        SET_SOURCE_FILES_PROPERTIES(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS "-msse4.1")
    ELSE()
        SET_SOURCE_FILES_PROPERTIES(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS "")
    ENDIF()
ENDIF()

# Include FFTW related files.
INCLUDE_DIRECTORIES(${FFTW_INCLUDES})

# Build the shared plugin library.
IF (OPENMM_BUILD_SHARED_LIB)
    ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})

    TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${OPENMM_LIBRARY_NAME} ${PTHREADS_LIB} ${FFTW_LIBRARY})
    IF (FFTW_THREADS_LIBRARY)
        TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${FFTW_THREADS_LIBRARY})
    ENDIF (FFTW_THREADS_LIBRARY)
    SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DOPENMM_PME_BUILDING_SHARED_LIBRARY")

    INSTALL_TARGETS(/lib/plugins RUNTIME_DIRECTORY /lib/plugins ${SHARED_TARGET})
ENDIF (OPENMM_BUILD_SHARED_LIB)

# Build the static plugin library.
IF(OPENMM_BUILD_STATIC_LIB)
    ADD_LIBRARY(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_INCLUDE_FILES})

    TARGET_LINK_LIBRARIES(${STATIC_TARGET} ${OPENMM_LIBRARY_NAME}_static ${PTHREADS_LIB} ${FFTW_LIBRARY})
    IF (FFTW_THREADS_LIBRARY)
        TARGET_LINK_LIBRARIES(${STATIC_TARGET} ${FFTW_THREADS_LIBRARY})
    ENDIF (FFTW_THREADS_LIBRARY)
    SET_TARGET_PROPERTIES(${STATIC_TARGET} PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DOPENMM_PME_BUILDING_STATIC_LIBRARY")

    INSTALL_TARGETS(/lib/plugins RUNTIME_DIRECTORY /lib/plugins ${STATIC_TARGET})
ENDIF(OPENMM_BUILD_STATIC_LIB)

IF(BUILD_TESTING)
    SUBDIRS(tests)
ENDIF(BUILD_TESTING)
