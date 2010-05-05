#
# Testing
#

ENABLE_TESTING()

# INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../FindOpenCL.cmake)
INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})

# Automatically create tests using files named "Test*.cpp"
FILE(GLOB TEST_PROGS "*Test*.cpp")
FOREACH(TEST_PROG ${TEST_PROGS})
    GET_FILENAME_COMPONENT(TEST_ROOT ${TEST_PROG} NAME_WE)

    # Link with shared library
    ADD_EXECUTABLE(${TEST_ROOT} ${TEST_PROG})
    TARGET_LINK_LIBRARIES(${TEST_ROOT} ${SHARED_TARGET})
    ADD_TEST(${TEST_ROOT} ${EXECUTABLE_OUTPUT_PATH}/${TEST_ROOT})

    # Link with static library
#     SET(TEST_STATIC ${TEST_ROOT}Static)
#     ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
#     SET_TARGET_PROPERTIES(${TEST_STATIC}
#                 PROPERTIES
#                 COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
#                 )
#     TARGET_LINK_LIBRARIES(${TEST_STATIC} ${STATIC_TARGET})
#     ADD_TEST(${TEST_STATIC} ${EXECUTABLE_OUTPUT_PATH}/${TEST_STATIC})

ENDFOREACH(TEST_PROG ${TEST_PROGS})

# TestOpenCLUsingParameterFile customized w/ command-line argument (input file name used in test)

#ADD_EXECUTABLE(TestOpenCLUsingParameterFile TstOpenCLUsingParameterFile.cpp)
#TARGET_LINK_LIBRARIES(TestOpenCLUsingParameterFile ${SHARED_TARGET})
#ADD_TEST(TestOpenCLUsingParameterFile "${EXECUTABLE_OUTPUT_PATH}/TestOpenCLUsingParameterFile" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/lambdaSdObcParameters.txt")
#ADD_TEST(TestOpenCLUsingParameterFile "${EXECUTABLE_OUTPUT_PATH}/TestOpenCLUsingParameterFile" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/bptiMdRfNoPbcParameters.txt")

SET(TEST_ROOT TestOpenCLUsingParameterFile)
SET(TEST_PROG TstOpenCLUsingParameterFile.cpp)
SET(TEST_STATIC ${TEST_ROOT}Static)
SET(INCLUDE_OPENCL_STATIC 1)
#IF(INCLUDE_OPENCL_STATIC)
#   ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
#   SET_TARGET_PROPERTIES(${TEST_STATIC}
#                         PROPERTIES
#                         COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
#                        )
#   TARGET_LINK_LIBRARIES(${TEST_STATIC} ${STATIC_TARGET} ${STATIC_BROOK_TARGET})
#   ADD_TEST(${TEST_STATIC} "${EXECUTABLE_OUTPUT_PATH}/TestOpenCLUsingParameterFileStatic" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/lambdaSdObcParameters.txt")
#   ADD_TEST(${TEST_STATIC} "${EXECUTABLE_OUTPUT_PATH}/TestOpenCLUsingParameterFileStatic" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/bptiMdRfNoPbcParameters.txt")
#  ADD_TEST(${TEST_STATIC} "${EXECUTABLE_OUTPUT_PATH}/TestOpenCLUsingParameterFileStatic" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/bptiMdRfPbcParameters.txt" " +checkEnergyForceConsistent -checkForces" )
#ENDIF(INCLUDE_OPENCL_STATIC)


