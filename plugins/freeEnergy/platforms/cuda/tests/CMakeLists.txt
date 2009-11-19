#
# Testing
#

ENABLE_TESTING()

INCLUDE(${CMAKE_SOURCE_DIR}/platforms/cuda/cuda-cmake/FindCuda.cmake)
INCLUDE_DIRECTORIES(${CUDA_INCLUDE})
INCLUDE_DIRECTORIES(${OPENMM_DIR}/platforms/cuda/include)
INCLUDE_DIRECTORIES(${OPENMM_DIR}/platforms/cuda/src)
INCLUDE_DIRECTORIES(${OPENMM_DIR}/platforms/cuda/src/kernels)
Set( SHARED_OPENMM_TARGET OpenMMFreeEnergy)
Set( STATIC_OPENMM_TARGET OpenMMFreeEnergy_static)
Set( SHARED_CUDA_TARGET OpenMMCuda)
Set( STATIC_CUDA_TARGET OpenMMCuda_static)
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_CUDA_TARGET   ${SHARED_CUDA_TARGET}_d)
    SET(SHARED_OPENMM_TARGET ${SHARED_OPENMM_TARGET}_d)
    SET(STATIC_CUDA_TARGET ${STATIC_CUDA_TARGET}_d)
    Set(STATIC_OPENMM_TARGET ${STATIC_OPENMM_TARGET}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
#LINK_DIRECTORIES


# Automatically create tests using files named "Test*.cpp"
FILE(GLOB TEST_PROGS "*Test*.cpp")
FOREACH(TEST_PROG ${TEST_PROGS})
    GET_FILENAME_COMPONENT(TEST_ROOT ${TEST_PROG} NAME_WE)

    # Link with shared library
    CUDA_ADD_EXECUTABLE(${TEST_ROOT} ${TEST_PROG})
    TARGET_LINK_LIBRARIES(${TEST_ROOT} ${SHARED_TARGET} ${SHARED_OPENMM_TARGET} ${SHARED_CUDA_TARGET})
    ADD_TEST(${TEST_ROOT} ${EXECUTABLE_OUTPUT_PATH}/${TEST_ROOT})



    # Link with static library
     SET(TEST_STATIC ${TEST_ROOT}Static)
     CUDA_ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
     SET_TARGET_PROPERTIES(${TEST_STATIC}
                 PROPERTIES
                 COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
                 )
     TARGET_LINK_LIBRARIES(${TEST_STATIC} ${STATIC_TARGET} ${STATIC_OPENMM_TARGET} ${STATIC_CUDA_TARGET})
     ADD_TEST(${TEST_STATIC} ${EXECUTABLE_OUTPUT_PATH}/${TEST_STATIC})

ENDFOREACH(TEST_PROG ${TEST_PROGS})

# TestCudaUsingParameterFile customized w/ command-line argument (input file name used in test) 

#ADD_EXECUTABLE(TestCudaUsingParameterFile TstCudaUsingParameterFile.cpp)
#TARGET_LINK_LIBRARIES(TestCudaUsingParameterFile ${SHARED_TARGET})
#ADD_TEST(TestCudaUsingParameterFile "${EXECUTABLE_OUTPUT_PATH}/TestCudaUsingParameterFile" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/lambdaSdObcParameters.txt")
#ADD_TEST(TestCudaUsingParameterFile "${EXECUTABLE_OUTPUT_PATH}/TestCudaUsingParameterFile" "-parameterFileName" "${CMAKE_CURRENT_SOURCE_DIR}/bptiMdRfNoPbcParameters.txt")
#
#SET(TEST_ROOT TestCudaUsingParameterFile)
#SET(TEST_PROG TstCudaUsingParameterFile.cpp)
#SET(TEST_STATIC ${TEST_ROOT}Static)
#SET(INCLUDE_CUDA_STATIC 1)
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


