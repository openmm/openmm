#---------------------------------------------------
# OpenMM UMA Plugin
#
# Creates OpenMMUMA plugin library for Universal
# Molecular Atomistic (UMA) potential integration
# using libtorch C++ API.
#----------------------------------------------------

# The source is organized into subdirectories
SET(OPENMM_UMA_PLUGIN_SOURCE_SUBDIRS . openmmapi)

SET(OPENMM_UMA_LIBRARY_NAME OpenMMUMA)
SET(SHARED_UMA_TARGET ${OPENMM_UMA_LIBRARY_NAME})
SET(STATIC_UMA_TARGET ${OPENMM_UMA_LIBRARY_NAME}_static)

# Find required dependencies
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# API include directories
SET(API_UMA_INCLUDE_DIRS)
FOREACH(subdir ${OPENMM_UMA_PLUGIN_SOURCE_SUBDIRS})
    SET(API_UMA_INCLUDE_DIRS ${API_UMA_INCLUDE_DIRS}
                             ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include
                             ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/include/internal)
ENDFOREACH(subdir)

# Collect header files
SET(API_UMA_REL_INCLUDE_FILES)
SET(API_UMA_ABS_INCLUDE_FILES)

FOREACH(dir ${API_UMA_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)
    SET(API_UMA_ABS_INCLUDE_FILES ${API_UMA_ABS_INCLUDE_FILES} ${fullpaths})
    
    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_UMA_REL_INCLUDE_FILES ${API_UMA_REL_INCLUDE_FILES} ${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# Collect source files
SET(API_UMA_SOURCE_FILES)
SET(API_UMA_SOURCE_INCLUDE_FILES)

FOREACH(subdir ${OPENMM_UMA_PLUGIN_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.cpp)
    FILE(GLOB incl_files ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/src/*.h)
    SET(API_UMA_SOURCE_FILES         ${API_UMA_SOURCE_FILES}         ${src_files})
    SET(API_UMA_SOURCE_INCLUDE_FILES ${API_UMA_SOURCE_INCLUDE_FILES} ${incl_files})
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/openmmapi/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/reference/include)
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/platforms/reference/src)
INCLUDE_DIRECTORIES(${TORCH_INCLUDE_DIRS})

# Create the library
ADD_LIBRARY(${SHARED_UMA_TARGET} SHARED ${API_UMA_SOURCE_FILES} ${API_UMA_SOURCE_INCLUDE_FILES} ${API_UMA_ABS_INCLUDE_FILES})

TARGET_LINK_LIBRARIES(${SHARED_UMA_TARGET} ${OPENMM_LIBRARY_NAME} ${TORCH_LIBRARIES})
TARGET_LINK_LIBRARIES(${SHARED_UMA_TARGET} ${PTHREADS_LIB})

SET_TARGET_PROPERTIES(${SHARED_UMA_TARGET} PROPERTIES
    COMPILE_FLAGS "-DOPENMM_UMA_BUILDING_SHARED_LIBRARY ${EXTRA_COMPILE_FLAGS}"
    LINK_FLAGS "${EXTRA_LINK_FLAGS}")

INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${SHARED_UMA_TARGET})

# Install headers
FILE(GLOB CORE_HEADERS     include/*.h          */include/*.h)
FILE(GLOB TOP_HEADERS      include/openmm/*.h)
FILE(GLOB INTERNAL_HEADERS include/openmm/internal/*.h */include/openmm/internal/*.h)
INSTALL_FILES(/include                 FILES ${CORE_HEADERS})
INSTALL_FILES(/include/openmm          FILES ${TOP_HEADERS})
INSTALL_FILES(/include/openmm/internal FILES ${INTERNAL_HEADERS})

# Build subdirectories
ADD_SUBDIRECTORY(platforms)

# Enable testing if requested
IF(BUILD_TESTING AND OPENMM_BUILD_C++_TESTS)
    ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTING AND OPENMM_BUILD_C++_TESTS)

# Python wrapper
IF(OPENMM_BUILD_PYTHON_WRAPPERS)
    ADD_SUBDIRECTORY(wrappers/python)
ENDIF(OPENMM_BUILD_PYTHON_WRAPPERS)
