#
# Test each supported platform with TestPlatform*.cpp programs
#

ENABLE_TESTING()

# For a particular platform, for example "ReferencePlatform", there must be:
#  1) A C++ object with that name derived from OpenMM::Platform
#  2) A header file with name like "ReferencePlatform.h"
SET(SUPPORTED_PLATFORMS "ReferencePlatform" CACHE STRING "Hardware platforms supported on this computer")

# Automatically create tests using files named "TestPlatform*.cpp.in"
FILE(GLOB PLATFORM_GENERIC_PROGS "TestPlatform*.cpp.in")

FOREACH(PLATFORM_GENERIC_PROG ${PLATFORM_GENERIC_PROGS})
    # Create platform-specific test for each supported platform
    FOREACH(OPENMM_PLATFORM_NAME ${SUPPORTED_PLATFORMS})

        STRING(REGEX REPLACE "TestPlatform([^\\/]*)\\.cpp.in" "\\1" TEST_NAME ${PLATFORM_GENERIC_PROG})
    	GET_FILENAME_COMPONENT(TEST_NAME ${TEST_NAME} NAME_WE)
        SET(PLATFORM_SPECIFIC_PROG "${CMAKE_CURRENT_BINARY_DIR}/Test${OPENMM_PLATFORM_NAME}${TEST_NAME}.cpp")
        CONFIGURE_FILE(${PLATFORM_GENERIC_PROG} ${PLATFORM_SPECIFIC_PROG} @ONLY)
        
    	GET_FILENAME_COMPONENT(TEST_EXECUTABLE ${PLATFORM_SPECIFIC_PROG} NAME_WE)
        
	    # Link with shared library
	    ADD_EXECUTABLE(${TEST_EXECUTABLE} ${PLATFORM_SPECIFIC_PROG})
	    TARGET_LINK_LIBRARIES(${TEST_EXECUTABLE} ${SHARED_TARGET})
	    ADD_TEST(${TEST_EXECUTABLE} ${EXECUTABLE_OUTPUT_PATH}/${TEST_EXECUTABLE})
        
    ENDFOREACH(OPENMM_PLATFORM_NAME ${SUPPORTED_PLATFORMS})
ENDFOREACH(PLATFORM_GENERIC_PROG ${PLATFORM_GENERIC_PROGS})

