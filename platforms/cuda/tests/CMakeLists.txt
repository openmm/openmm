#
# Testing
#

ENABLE_TESTING()

INCLUDE(${CMAKE_SOURCE_DIR}/platforms/cuda/cuda-cmake/FindCuda.cmake)
INCLUDE_DIRECTORIES(${CUDA_INCLUDE})

# Automatically create tests using files named "Test*.cpp"
FILE(GLOB TEST_PROGS "*Test*.cpp")
# TestCudaRandom has never worked on windows, so let's stop polluting the dashboard with it.
IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
    FILE(GLOB TEST_CUDA_RANDOM_PROG "*TestCudaRandom.cpp")
    LIST(REMOVE_ITEM TEST_PROGS "${TEST_CUDA_RANDOM_PROG}")
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
FOREACH(TEST_PROG ${TEST_PROGS})
    GET_FILENAME_COMPONENT(TEST_ROOT ${TEST_PROG} NAME_WE)

    # Link with shared library
    CUDA_ADD_EXECUTABLE(${TEST_ROOT} ${TEST_PROG})
    TARGET_LINK_LIBRARIES(${TEST_ROOT} ${SHARED_TARGET})
    ADD_TEST(${TEST_ROOT} ${EXECUTABLE_OUTPUT_PATH}/${TEST_ROOT})

    # Link with static library
#     SET(TEST_STATIC ${TEST_ROOT}Static)
#     CUDA_ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
#     SET_TARGET_PROPERTIES(${TEST_STATIC}
#                 PROPERTIES
#                 COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
#                 )
#     TARGET_LINK_LIBRARIES(${TEST_STATIC} ${STATIC_TARGET})
#     ADD_TEST(${TEST_STATIC} ${EXECUTABLE_OUTPUT_PATH}/${TEST_STATIC})

ENDFOREACH(TEST_PROG ${TEST_PROGS})

# TestCudaUsingParameterFile customized w/ command-line argument (input file name used in test) 

ADD_EXECUTABLE(TestCudaUsingParameterFile TstCudaUsingParameterFile.cpp)
TARGET_LINK_LIBRARIES(TestCudaUsingParameterFile ${SHARED_TARGET})
ADD_TEST(TestCudaUsingParameterFile "${EXECUTABLE_OUTPUT_PATH}/TestCudaUsingParameterFile" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/lambdaSdObcParameters.txt")
ADD_TEST(TestCudaUsingParameterFile "${EXECUTABLE_OUTPUT_PATH}/TestCudaUsingParameterFile" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/bptiMdRfNoPbcParameters.txt")

SET(TEST_ROOT TestCudaUsingParameterFile)
SET(TEST_PROG TstCudaUsingParameterFile.cpp)
SET(TEST_STATIC ${TEST_ROOT}Static)
SET(INCLUDE_CUDA_STATIC 1)
#IF(INCLUDE_CUDA_STATIC)
#   ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
#   SET_TARGET_PROPERTIES(${TEST_STATIC}
#                         PROPERTIES
#                         COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
#                        )
#   TARGET_LINK_LIBRARIES(${TEST_STATIC} ${STATIC_TARGET} ${STATIC_BROOK_TARGET})
#   ADD_TEST(${TEST_STATIC} "${EXECUTABLE_OUTPUT_PATH}/TestCudaUsingParameterFileStatic" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/lambdaSdObcParameters.txt")
#   ADD_TEST(${TEST_STATIC} "${EXECUTABLE_OUTPUT_PATH}/TestCudaUsingParameterFileStatic" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/bptiMdRfNoPbcParameters.txt")
#  ADD_TEST(${TEST_STATIC} "${EXECUTABLE_OUTPUT_PATH}/TestCudaUsingParameterFileStatic" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/bptiMdRfPbcParameters.txt" " +checkEnergyForceConsistent -checkForces" )
#ENDIF(INCLUDE_CUDA_STATIC)


