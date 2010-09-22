#---------------------------------------------------
# OpenMM Serialization Library
#
# Creates OpenMM serializatin library, base name=OpenMMSerialization.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   OpenMMSerialization[_d].dll
#   OpenMMSerialization[_d].lib
# Unix:
#   libOpenMMSerialization[_d].so
#----------------------------------------------------

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)


# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMM_SERIALIZATION_LIBRARY_NAME OpenMMSerialization)


# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.
IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(OPENMM_SERIALIZATION_LIBRARY_NAME ${OPENMM_SERIALIZATION_LIBRARY_NAME}_d)
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

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    FILE(GLOB_RECURSE src_files   ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.c)
    FILE(GLOB incl_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.hpp)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    INCLUDE_DIRECTORIES(BEFORE  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Create the library
ADD_LIBRARY(${OPENMM_SERIALIZATION_LIBRARY_NAME} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
TARGET_LINK_LIBRARIES(${OPENMM_SERIALIZATION_LIBRARY_NAME} debug ${OPENMM_LIBRARY_NAME}_d optimized ${OPENMM_LIBRARY_NAME})
SET_TARGET_PROPERTIES(${OPENMM_SERIALIZATION_LIBRARY_NAME} PROPERTIES COMPILE_FLAGS "-DOPENMM_BUILDING_SHARED_LIBRARY -DTIXML_USE_STL -DIEEE_8087")
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${OPENMM_SERIALIZATION_LIBRARY_NAME})

ADD_SUBDIRECTORY(tests)
