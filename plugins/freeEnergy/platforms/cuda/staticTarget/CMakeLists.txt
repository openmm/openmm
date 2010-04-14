#
# Include CUDA related files.
#
SET(OPENMM_BUILD_FREE_ENERGY_PATH ${CMAKE_SOURCE_DIR}/plugins/freeEnergy)

SET(CUDA_NVCC_BUILD_FLAGS "-DCUDPP_STATIC_LIB" )
# INCLUDE(${FINDCUDA_DIR}/FindCuda.cmake)
INCLUDE_DIRECTORIES(${CUDA_INCLUDE})
LINK_DIRECTORIES(${CUDA_TARGET_LINK})
FOREACH(subdir ${OPENMM_SOURCE_SUBDIRS})
    FILE(GLOB src_files ${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/${subdir}/src/*.cu  ${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/${subdir}/src/*/*.cu)
    FOREACH(file ${src_files})
        FILE(RELATIVE_PATH file ${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda ${file})
        SET(SOURCE_FILES         ${SOURCE_FILES}         ${file})   #append
    ENDFOREACH(file)
    CUDA_INCLUDE_DIRECTORIES(BEFORE ${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/${subdir}/include)
    CUDA_INCLUDE_DIRECTORIES(BEFORE ${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/${subdir}/src)
ENDFOREACH(subdir)
#CUDA_INCLUDE_DIRECTORIES(BEFORE ${OPENMM_BUILD_FREE_ENERGY_PATH}/jama/include)
CUDA_INCLUDE_DIRECTORIES(${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/../src
                         ${OPENMM_DIR}/platforms/cuda/src
                         ${OPENMM_DIR}/platforms/cuda/include
                         ${OPENMM_DIR}/platforms/cuda/src/kernels
                         ${OPENMM_DIR}/openmmapi/include ) 

CUDA_ADD_LIBRARY(${STATIC_TARGET} STATIC ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

# required for getting OPENMM_EXPORT to be set correctly in 'class OPENMM_EXPORT CudaStreamFactory', ...
# see OpenMM/openmmapi/include/internal/windowsExport.h for details

SET(CUDA_STATIC_COMPILE_FLAG "-DOPENMMCUDAFREEENERGY_BUILDING_STATIC_LIBRARY -DOPENMMCUDAFREEENERGY_USE_STATIC_LIBRARIES -DCUDPP_STATIC_LIB")
SET_TARGET_PROPERTIES(${STATIC_TARGET} PROPERTIES COMPILE_FLAGS ${CUDA_STATIC_COMPILE_FLAG})

TARGET_LINK_LIBRARIES(${STATIC_TARGET} optimized ${OPENMM_LIBRARY_NAME}_static)
TARGET_LINK_LIBRARIES(${STATIC_TARGET} debug ${OPENMM_LIBRARY_NAME}_static_d optimized ${OPENMM_LIBRARY_NAME}_static ${CUFFT_TARGET_LINK})

INSTALL_TARGETS(/lib/plugins RUNTIME_DIRECTORY /lib/plugins ${STATIC_TARGET})
