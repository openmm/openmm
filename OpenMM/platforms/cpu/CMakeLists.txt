#---------------------------------------------------
# OpenMM CPU Platform
#
# Creates OpenMMCPU library.
#
# Windows:
#   OpenMMCPU.dll
#   OpenMMCPU.lib
#   OpenMMCPU_static.lib
# Unix:
#   libOpenMMCPU.so
#   libOpenMMCPU_static.a
#----------------------------------------------------

SET(OPENMM_BUILD_CPU_TESTS TRUE CACHE BOOL "Whether to build CPU platform test cases")
MARK_AS_ADVANCED(OPENMM_BUILD_CPU_TESTS)
IF(BUILD_TESTING AND OPENMM_BUILD_CPU_TESTS)
    SUBDIRS(tests)
ENDIF(BUILD_TESTING AND OPENMM_BUILD_CPU_TESTS)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(OPENMM_SOURCE_SUBDIRS .)

IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_COMPILER_IS_CLANGXX 1)
ENDIF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")

IF(CMAKE_COMPILER_IS_CLANGXX)
  EXECUTE_PROCESS( COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE __clang_full_version_string )
  string (REGEX REPLACE ".*clang version ([0-9]+\\.[0-9]+).*" "\\1" CLANG_VERSION_STRING ${__clang_full_version_string})
  unset(__clang_full_version_string)
ENDIF(CMAKE_COMPILER_IS_CLANGXX)

IF(CMAKE_COMPILER_IS_CLANGXX AND CLANG_VERSION_STRING STREQUAL "3.6")
 message(FATAL_ERROR  "The OpenMM CPU platform cannot be built with your current compiler, clang-3.6, due to bugs \
in the compiler's AVX support. Either downgrade to clang-3.5, upgrade to clang-3.7 or later, switch to an \
alternative compiler such as GCC, or turn off building of the CPU platform by unsetting the CMake variable \
OPENMM_BUILD_CPU_LIB.")
ENDIF(CMAKE_COMPILER_IS_CLANGXX AND CLANG_VERSION_STRING STREQUAL "3.6")


# Collect up information about the version of the OpenMM library we're building
# and make it available to the code so it can be built into the binaries.

SET(OPENMMCPU_LIBRARY_NAME OpenMMCPU)

SET(SHARED_TARGET ${OPENMMCPU_LIBRARY_NAME})
SET(STATIC_TARGET ${OPENMMCPU_LIBRARY_NAME}_static)


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
    IF(MSVC)
         FILE(GLOB_RECURSE kernel_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/kernels/*.cu)
         SET(SOURCE_FILES ${SOURCE_FILES} ${kernel_files})
    ENDIF(MSVC)
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Install headers

FILE(GLOB CORE_HEADERS include/*.h)
INSTALL_FILES(/include/openmm/cpu FILES ${CORE_HEADERS})

IF(OPENMM_BUILD_SHARED_LIB)
    SUBDIRS (sharedTarget)
ENDIF(OPENMM_BUILD_SHARED_LIB)
IF(OPENMM_BUILD_STATIC_LIB)
    SUBDIRS (staticTarget)
ENDIF(OPENMM_BUILD_STATIC_LIB)
