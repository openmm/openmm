#---------------------------------------------------
# OpenMM 
#
# Creates OpenMM library, base name=OpenMM.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d).
#
# Windows:
#   OpenMM[_d].dll
#   OpenMM[_d].lib
#   OpenMM_static[_d].lib
# Unix:
#   libOpenMM[_d].so
#   libOpenMM_static[_d].a
#----------------------------------------------------

# Don't create a new project name if this is part of a mega-build from the
# parent directory
IF( NOT PROJECT_NAME )
    PROJECT (OpenMM)
ENDIF( NOT PROJECT_NAME )

# On Linux it appears that cmake 2.4 does not work with Cuda cmake rules
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#SET(CMAKE_VERBOSE_MAKEFILE 1)

# CMake 2.4 on dart server only needs to know about Dart
INCLUDE(Dart)
MARK_AS_ADVANCED(DART_ROOT)

# Avoid dart server problem with cmake 2.4 --cmb
set(cmv "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")
IF(NOT cmv EQUAL "2.4") # This whole file...

# We have custom cmake modules for FindOpenMM and running python tests
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

# Older cmake versions do not have cmake_policy command
# Cmake 2.4.5, the default cmake on red hat linux, has the
# cmake_policy command, but it does not work
# "if(POLICY ..." does not work with cmake 2.4.[56] on red hat
# (cmake 2.4.7 is OK)
if(COMMAND cmake_policy)
    if(CMAKE_MAJOR_VERSION GREATER 2 OR CMAKE_MINOR_VERSION GREATER 5)
        cmake_policy(SET CMP0003 NEW) 
        cmake_policy(SET CMP0005 NEW) 
    endif(CMAKE_MAJOR_VERSION GREATER 2 OR CMAKE_MINOR_VERSION GREATER 5)
    if(CMAKE_MAJOR_VERSION GREATER 2 OR CMAKE_MINOR_VERSION GREATER 6 OR CMAKE_PATCH_VERSION GREATER 2)
        cmake_policy(SET CMP0011 NEW) 
    endif(CMAKE_MAJOR_VERSION GREATER 2 OR CMAKE_MINOR_VERSION GREATER 6 OR CMAKE_PATCH_VERSION GREATER 2)
endif(COMMAND cmake_policy)

# Where to install
IF(${CMAKE_C_COMPILER} MATCHES "gcc")
    IF(NOT OPENMM_INSTALL_PREFIX)
        SET(OPENMM_INSTALL_PREFIX "/usr/local/openmm" CACHE PATH "Where to install OpenMM")
    ENDIF(NOT OPENMM_INSTALL_PREFIX)
ELSE(${CMAKE_C_COMPILER} MATCHES "gcc")
    IF(NOT OPENMM_INSTALL_PREFIX)
        SET(OPENMM_INSTALL_PREFIX "$ENV{ProgramFiles}/OpenMM" CACHE PATH "Where to install OpenMM")
    ENDIF(NOT OPENMM_INSTALL_PREFIX)
ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")
MARK_AS_ADVANCED(OPENMM_INSTALL_PREFIX)

# It seems that on linux and mac, everything is trying to be installed in /usr/local/openmm
# But if every install target is prefixed with /openmm/, on Windows the install files
# end up in C:/Program Files/OpenMM/openmm/ which is ugly.
# Better to set CMAKE_INSTALL_PREFIX to /usr/local/openmm and leave /openmm/ off the 
# install target names.  Plus, the user now has the opportunity to install into /usr/local/
# if she so chooses.  --cmb
IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
  SET(CMAKE_INSTALL_PREFIX ${OPENMM_INSTALL_PREFIX} CACHE PATH "Where to install OpenMM" FORCE)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS . openmmapi olla libraries/jama libraries/quern libraries/lepton libraries/sfmt platforms/reference libraries/validate)

# The build system will set ARCH64 for 64 bit builds, which require
# use of the lib64/ library directories rather than lib/.
#SET( ARCH64 OFF CACHE BOOL "ON for 64bit builds, OFF for 32bit builds")
#MARK_AS_ADVANCED( ARCH64  )
#IF (ARCH64)
#    SET(LIB64 64)
#ELSE (ARCH64)
#    SET(LIB64) # nothing
#ENDIF (ARCH64)

IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  SET( LIB64 64 )
ELSE( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  SET( LIB64  )
ENDIF( CMAKE_SIZEOF_VOID_P EQUAL 8 )

# Build 32 bit binaries, since CUDA doesn't currently work with 64 bit
IF (APPLE)
   SET (CMAKE_CXX_FLAGS "-arch i386 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5")
   SET (CMAKE_C_FLAGS "-arch i386 -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5")
ENDIF (APPLE)

IF(UNIX AND NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release build" FORCE)
ENDIF (UNIX AND NOT CMAKE_BUILD_TYPE)

IF (NOT CMAKE_CXX_FLAGS_DEBUG)
    SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "To use when CMAKE_BUILD_TYPE=Debug" FORCE)
ENDIF (NOT CMAKE_CXX_FLAGS_DEBUG)

IF (NOT CMAKE_CXX_FLAGS_RELEASE)
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING 
          "To use when CMAKE_BUILD_TYPE=Release" FORCE)
ENDIF (NOT CMAKE_CXX_FLAGS_RELEASE)



# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMM_LIBRARY_NAME OpenMM)
SET(OPENMM_MAJOR_VERSION 1)
SET(OPENMM_MINOR_VERSION 0)
SET(OPENMM_BUILD_VERSION 0)

SET(OPENMM_COPYRIGHT_YEARS "2008")

# underbar separated list of dotted authors, no spaces or commas
SET(OPENMM_AUTHORS "Peter.Eastman")

# Get the subversion revision number if we can
# It's possible that WIN32 installs use svnversion through cygwin
# so we'll try for both svnversion.exe and svnversion. Note that
# this will result in warnings if all you have is Tortoise without
# Cygwin, and your "about" string will say "unknown" rather than
# providing the SVN version of the source.

FIND_PROGRAM (SVNVERSION_PROGRAM svnversion
    PATHS "C:/cygwin/bin"
)
IF (SVNVERSION_PROGRAM)
    EXEC_PROGRAM (${SVNVERSION_PROGRAM}
                  ARGS \"${CMAKE_CURRENT_SOURCE_DIR}\"
                  OUTPUT_VARIABLE OPENMM_SVN_REVISION)
ELSE (SVNVERSION_PROGRAM)
    MESSAGE (STATUS 
      "Could not find 'svnversion' executable; 'about' will be wrong. (Cygwin provides one on Windows.)"
    )
    SET (OPENMM_SVN_REVISION unknown) 
ENDIF (SVNVERSION_PROGRAM)

# Remove colon from build version, for easier placement in directory names
STRING(REPLACE ":" "_" OPENMM_SVN_REVISION ${OPENMM_SVN_REVISION})
SET(OPENMM_SVN_REVISION "${OPENMM_SVN_REVISION}" CACHE STRING "subversion repository revision of OpenMM" FORCE)

ADD_DEFINITIONS(-DOPENMM_LIBRARY_NAME=${OPENMM_LIBRARY_NAME}
                -DOPENMM_MAJOR_VERSION=${OPENMM_MAJOR_VERSION}
                -DOPENMM_MINOR_VERSION=${OPENMM_MINOR_VERSION}
		-DOPENMM_BUILD_VERSION=${OPENMM_BUILD_VERSION})

# CMake quotes automatically when building Visual Studio projects but we need
# to add them ourselves for Linux or Cygwin. Two cases to avoid duplicate quotes
# in Visual Studio which end up in the binary.

IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
   SET(NEED_QUOTES FALSE)
ELSE (${CMAKE_GENERATOR} MATCHES "Visual Studio")
   SET(NEED_QUOTES TRUE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

##TODO: doesn't work without quotes in nightly build
SET(NEED_QUOTES TRUE)

IF(NEED_QUOTES)
   ADD_DEFINITIONS(-DOPENMM_SVN_REVISION="${OPENMM_SVN_REVISION}"
                   -DOPENMM_COPYRIGHT_YEARS="${OPENMM_COPYRIGHT_YEARS}"
                   -DOPENMM_AUTHORS="${OPENMM_AUTHORS}")
ELSE(NEED_QUOTES)
   ADD_DEFINITIONS(-DOPENMM_SVN_REVISION=${OPENMM_SVN_REVISION}
                   -DOPENMM_COPYRIGHT_YEARS=${OPENMM_COPYRIGHT_YEARS}
                   -DOPENMM_AUTHORS=${OPENMM_AUTHORS})
ENDIF(NEED_QUOTES)

# -DOPENMM_TYPE has to be defined in the target subdirectories.
# -Dsimbody_EXPORTS defined automatically when Windows DLL build is being done.

# Report the version number to the CMake UI
SET(OPENMM_VERSION
    "${OPENMM_MAJOR_VERSION}.${OPENMM_MINOR_VERSION}.${OPENMM_BUILD_VERSION}"
     CACHE STRING "This is the version of OpenMM which will be built." FORCE)
# We don't actually use this version for anything yet
MARK_AS_ADVANCED(OPENMM_VERSION)

SET(SHARED_TARGET ${OPENMM_LIBRARY_NAME})
SET(STATIC_TARGET ${OPENMM_LIBRARY_NAME}_static)


## If no one says otherwise, change the executable path to drop into the same binary
## location as the DLLs so that the test cases will use the just-build DLLs.
IF(NOT EXECUTABLE_OUTPUT_PATH)
  SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR} 
      CACHE INTERNAL "Single output directory for building all executables.")
ENDIF(NOT EXECUTABLE_OUTPUT_PATH)
IF(NOT LIBRARY_OUTPUT_PATH)
  SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR} 
      CACHE INTERNAL "Single output directory for building all libraries.")
ENDIF(NOT LIBRARY_OUTPUT_PATH)
SET(${PROJECT_NAME}_EXECUTABLE_DIR ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})
SET(${PROJECT_NAME}_LIBRARY_DIR    ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})

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

# used by plugin
SET(OPENMM_DIR ${CMAKE_CURRENT_SOURCE_DIR})

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
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in OpenMM/include if
    ## OpenMM was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

# If API wrappers are being generated, and add them to the build.
# Java and gccxml are required for API wrappers
FIND_PACKAGE(Java)
MARK_AS_ADVANCED(CLEAR JAVA_RUNTIME)
# One particular location to try first for nightly builds
FIND_PROGRAM(GCCXML_PATH gccxml PATH
    "C:/Program Files/gccxml_sherm/bin")
# In case first attempt fails
FIND_PROGRAM(GCCXML_PATH gccxml PATH
    /usr/local/bin
    "C:/Program Files/gccxml 0.9/bin"
)
IF(GCCXML_PATH AND JAVA_RUNTIME AND NOT cmv EQUAL "2.4")
    SET(OPENMM_BUILD_API_WRAPPERS ON CACHE BOOL "Build wrappers for C and Fortran")
ELSE(GCCXML_PATH AND JAVA_RUNTIME AND NOT cmv EQUAL "2.4")
    SET(OPENMM_BUILD_API_WRAPPERS OFF CACHE BOOL "Build wrappers for C and Fortran")
ENDIF(GCCXML_PATH AND JAVA_RUNTIME AND NOT cmv EQUAL "2.4")
IF(OPENMM_BUILD_API_WRAPPERS)
    ADD_SUBDIRECTORY(wrappers)
    SET(SOURCE_FILES ${SOURCE_FILES} wrappers/OpenMMCWrapper.cpp  wrappers/OpenMMFortranWrapper.cpp)
    SET_SOURCE_FILES_PROPERTIES(wrappers/OpenMMCWrapper.cpp wrappers/OpenMMFortranWrapper.cpp PROPERTIES GENERATED TRUE)
ENDIF(OPENMM_BUILD_API_WRAPPERS)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES COMPILE_FLAGS "-DOPENMM_BUILDING_SHARED_LIBRARY -DLEPTON_BUILDING_SHARED_LIBRARY -DOPENMM_VALIDATE_BUILDING_SHARED_LIBRARY")

SET(OPENMM_BUILD_STATIC_LIB OFF CACHE BOOL "Whether to build static OpenMM libraries")
IF(OPENMM_BUILD_STATIC_LIB)
    ADD_LIBRARY(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})
    SET_TARGET_PROPERTIES(${STATIC_TARGET} PROPERTIES COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES -DOPENMM_BUILDING_STATIC_LIBRARY -DLEPTON_USE_STATIC_LIBRARIES -DLEPTON_BUILDING_STATIC_LIBRARY -DOPENMMM_VALIDATE_BUILDING_STATIC_LIBRARY -DOPENMM_VALIDATE_BUILDING_STATIC_LIBRARY")
ENDIF(OPENMM_BUILD_STATIC_LIB)

IF(OPENMM_BUILD_API_WRAPPERS)
    ADD_DEPENDENCIES(${SHARED_TARGET} ApiWrappers)
    IF(OPENMM_BUILD_STATIC_LIB)
        ADD_DEPENDENCIES(${STATIC_TARGET} ApiWrappers)
    ENDIF(OPENMM_BUILD_STATIC_LIB)
ENDIF(OPENMM_BUILD_API_WRAPPERS)

# On Linux need to link to libdl
FIND_LIBRARY(DL_LIBRARY dl)
IF(DL_LIBRARY)
  TARGET_LINK_LIBRARIES(${SHARED_TARGET} ${DL_LIBRARY})
  IF(OPENMM_BUILD_STATIC_LIB)
    TARGET_LINK_LIBRARIES(${STATIC_TARGET} ${DL_LIBRARY})
  ENDIF(OPENMM_BUILD_STATIC_LIB)
ENDIF(DL_LIBRARY)
IF(WIN32)
    MARK_AS_ADVANCED(DL_LIBRARY)
ENDIF(WIN32)

ADD_SUBDIRECTORY(platforms/reference/tests)

# Which hardware platforms to build

# A bit of tedium because we are using custom FindCUDA files that happen to work...
SET(FINDCUDA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindCUDA")
SET(CUDA_BUILD_CUBIN OFF)
FIND_PACKAGE(CUDA QUIET)
IF(CUDA_FOUND)
    IF(NOT CUDA_NVCC_FLAGS)
        SET(FLAGS "")
        # Note that cmake will insert semicolons between these item automatically...
        SET(FLAGS ${FLAGS} -gencode arch=compute_11,code=sm_11)
        SET(FLAGS ${FLAGS} -gencode arch=compute_13,code=sm_13)
        SET(FLAGS ${FLAGS} -gencode arch=compute_20,code=sm_20)
        SET(FLAGS ${FLAGS} -use_fast_math)
        IF(MSVC)
            # Unfortunately the variables CUDA_NVCC_FLAGS_RELEASE and CUDA_NVCC_FLAGS_DEBUG
            # appear to be unused, at least in CMake 2.6
            # Release /MD linkage
            set(FLAGS ${FLAGS} "-Xcompiler \"/MD\"")
        ENDIF(MSVC)
        SET(CUDA_NVCC_FLAGS  "${FLAGS}"
            CACHE STRING "Semicolon delimit multiple arguments")
    ENDIF(NOT CUDA_NVCC_FLAGS)
    SET(OPENMM_BUILD_CUDA_LIB ON CACHE BOOL "Build OpenMMCuda library for Nvidia GPUs")

    # I wish I was not hardcoding /MD here
    # Avoid strange windows link error with cuda free energy 
    # plugin by specifying /MD to CUDA_NVCC_FLAGS
    MARK_AS_ADVANCED(CLEAR CUDA_NVCC_FLAGS)
    IF(MSVC)
        # Unfortunately the variables CUDA_NVCC_FLAGS_RELEASE and CUDA_NVCC_FLAGS_DEBUG
        # appear to be unused, at least in CMake 2.6
        # Release /MD linkage
        SET(HAS_NVCC_FLAG FALSE)
        IF(CUDA_NVCC_FLAGS MATCHES "-Xcompiler")
            SET(HAS_NVCC_FLAG TRUE)
        ENDIF(CUDA_NVCC_FLAGS MATCHES "-Xcompiler")
        IF(NOT HAS_NVCC_FLAG)
            set(new_flags "-Xcompiler \"/MD\"")
            if(CUDA_NVCC_FLAGS)
                set(new_flags "${CUDA_NVCC_FLAGS};${new_flags}")
            endif(CUDA_NVCC_FLAGS)
            SET(CUDA_NVCC_FLAGS "${new_flags}" 
                CACHE STRING "Semicolon delimit multiple arguments" FORCE)
        ENDIF(NOT HAS_NVCC_FLAG)
    ENDIF(MSVC)
ELSE(CUDA_FOUND)
    SET(OPENMM_BUILD_CUDA_LIB OFF CACHE BOOL "Build OpenMMCuda library for Nvidia GPUs")
ENDIF(CUDA_FOUND)
IF(OPENMM_BUILD_CUDA_LIB)
    ADD_SUBDIRECTORY(platforms/cuda)
ENDIF(OPENMM_BUILD_CUDA_LIB)
MARK_AS_ADVANCED(CUDA_VERBOSE_BUILD)
MARK_AS_ADVANCED(CUDA_BUILD_CUBIN)
MARK_AS_ADVANCED(CUDA_BUILD_EMULATION)

FIND_PACKAGE(OpenCL QUIET)
# Plus the tesla/linux is taking too long on the tests
IF(OPENCL_FOUND)
    SET(OPENMM_BUILD_OPENCL_LIB ON CACHE BOOL "Build OpenMMOpenCL library for Nvidia GPUs")
ELSE(OPENCL_FOUND)
    SET(OPENMM_BUILD_OPENCL_LIB OFF CACHE BOOL "Build OpenMMOpenCL library for Nvidia GPUs")
ENDIF(OPENCL_FOUND)
IF(OPENMM_BUILD_OPENCL_LIB)
    ADD_SUBDIRECTORY(platforms/opencl)
ENDIF(OPENMM_BUILD_OPENCL_LIB)

# FreeEnergy plugin

SET(OPENMM_BUILD_FREE_ENERGY_PLUGIN ON CACHE BOOL "Build FreeEnergy plugin")
SET(OPENMM_BUILD_FREE_ENERGY_PATH)
IF(OPENMM_BUILD_FREE_ENERGY_PLUGIN)
   SET(OPENMM_BUILD_FREE_ENERGY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/plugins/freeEnergy)
#  ADD_DEFINITIONS(-DINCLUDE_FREE_ENERGY_PLUGIN)
   ADD_SUBDIRECTORY(plugins/freeEnergy)
ENDIF(OPENMM_BUILD_FREE_ENERGY_PLUGIN)

INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_TARGET})
IF(OPENMM_BUILD_STATIC_LIB)
  INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${STATIC_TARGET})
ENDIF(OPENMM_BUILD_STATIC_LIB)
FILE(GLOB CORE_HEADERS     include/*.h          */include/*.h)
FILE(GLOB TOP_HEADERS      include/openmm/*.h          */include/openmm/*.h)
FILE(GLOB INTERNAL_HEADERS include/openmm/internal/*.h */include/openmm/internal/*.h )
INSTALL_FILES(/include                 FILES ${CORE_HEADERS})
INSTALL_FILES(/include/openmm          FILES ${TOP_HEADERS})
INSTALL_FILES(/include/openmm/internal FILES ${INTERNAL_HEADERS})

#
# Allow automated build and dashboard.
#
# redundant INCLUDE(Dart) statement
# INCLUDE (Dart)

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

# Build and install documentation

find_package(Doxygen QUIET)
mark_as_advanced(CLEAR DOXYGEN_EXECUTABLE)
IF(DOXYGEN_EXECUTABLE)
    SET(DOXY_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in 
          ${DOXY_CONFIG}
          @ONLY )
    FILE(GLOB_RECURSE OPENMM_INCLUDES "openmm/include/*.h")
    FILE(GLOB_RECURSE OLLA_INCLUDES "olla/include/*.h")
    ADD_CUSTOM_COMMAND(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG}
        DEPENDS ${OPENMM_INCLUDES} ${OLLA_INCLUDES}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Generating API documentation using Doxygen")
    ADD_CUSTOM_TARGET(DoxygenApiDocs 
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
        COMMENT "Generating API documentation using Doxygen"
        SOURCES 
            "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in"
            ${OPENMM_INCLUDES}
            ${OLLA_INCLUDES}
    ) 
    FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/html/")
    INSTALL(DIRECTORY "${PROJECT_BINARY_DIR}/html/"
            DESTINATION "docs/api/")
    INSTALL(FILES "docs/API Reference.html"
            DESTINATION "docs/")
  set(OPENMM_GENERATE_API_DOCS OFF CACHE BOOL "Whether to create API documentation using Doxygen")
  if(OPENMM_GENERATE_API_DOCS)
    ADD_CUSTOM_TARGET(GenerateDoxygenApiDocs ALL
       DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
       COMMENT "Generating API documentation using Doxygen"
       SOURCES 
            "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in"
            ${OPENMM_INCLUDES}
            ${OLLA_INCLUDES}
    )
  endif(OPENMM_GENERATE_API_DOCS)
ENDIF(DOXYGEN_EXECUTABLE)

install(FILES docs/UsersGuide/OpenMMUsersGuide.pdf
    DESTINATION docs/)

FILE(GLOB LICENSE_FILES "docs/licenses/*.txt")
install(FILES ${LICENSE_FILES}
    DESTINATION licenses/)

ADD_SUBDIRECTORY(tests)
ADD_SUBDIRECTORY(examples)

ENDIF(NOT cmv EQUAL "2.4") # This whole file...
