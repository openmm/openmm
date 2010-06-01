#
# Include CUDA related files.
#
# ----------------------------------------------------------------------------
# logging

SET(LOG TRUE)
   
IF(LOG)
   SET(LOG_FILE "CMakeLog.txt" )
   FILE( WRITE ${LOG_FILE} "In freeEnergy/platforms/cuda/sharedTarget Cmake\n")
#  FILE( APPEND ${LOG_FILE} "BROOK_LIB_PATH=${BROOK_LIB_PATH}\n")
ENDIF(LOG)

IF(LOG)
   MACRO(LOG_DIR LOG_FILE DIR_LIST )
   
      FILE( APPEND ${LOG_FILE} "\n${DIR_LIST}\n")
      FOREACH(currentFile ${ARGN})
         FILE( APPEND ${LOG_FILE} "   ${currentFile}\n" )
      ENDFOREACH(currentFile)
   
   ENDMACRO(LOG_DIR) 
ENDIF(LOG)

# ----------------------------------------------------------------------------

SET(OPENMM_BUILD_FREE_ENERGY_PATH ${CMAKE_SOURCE_DIR}/plugins/freeEnergy) 

# ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "Pre OPENMM_SOURCE_SUBDIRS" ${OPENMM_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "Pre OPENMM_FREE_ENERGY_SOURCE_SUBDIRS " ${OPENMM_FREE_ENERGY_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "Pre SOURCE_FILES" ${SOURCE_FILES} )
ENDIF(LOG)
## ----------------------------------------------------------------------------

# message("CUDA_NVCC_FLAGS = ${CUDA_NVCC_FLAGS}")
# INCLUDE(${FINDCUDA_DIR}/FindCuda.cmake)
# message("CUDA_NVCC_FLAGS = ${CUDA_NVCC_FLAGS}")
INCLUDE_DIRECTORIES(${CUDA_INCLUDE})
LINK_DIRECTORIES(${CUDA_TARGET_LINK})
FOREACH(subdir ${OPENMM_FREE_ENERGY_SOURCE_SUBDIRS})
    FILE(GLOB src_files ${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/${subdir}/src/*.cu ${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/src/*/*.cu)
    SET(SOURCE_FILES ${SOURCE_FILES} ${src_files})
    CUDA_INCLUDE_DIRECTORIES(BEFORE ${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/../${subdir}/include)
ENDFOREACH(subdir)

# ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "OPENMM_BUILD_FREE_ENERGY_PATH" ${OPENMM_BUILD_FREE_ENERGY_PATH} )
   FILE( APPEND ${LOG_FILE} "OPENMM_BUILD_FREE_ENERGY_PATH=${OPENMM_BUILD_FREE_ENERGY_PATH}\n")
   LOG_DIR( ${LOG_FILE} "OPENMM_SOURCE_SUBDIRS" ${OPENMM_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "CMAKE_SOURCE_DIR" ${CMAKE_SOURCE_DIR} )
   LOG_DIR( ${LOG_FILE} "CUDA_INCLUDE" ${CUDA_INCLUDE} )
   LOG_DIR( ${LOG_FILE} "CUDA_TARGET_LINK" ${CUDA_TARGET_LINK} )
   LOG_DIR( ${LOG_FILE} "SHARED_TARGET" ${SHARED_TARGET} )
   LOG_DIR( ${LOG_FILE} "OPENMM_DIR" ${OPENMM_DIR} )
   LOG_DIR( ${LOG_FILE} "SOURCE_FILES" ${SOURCE_FILES} )
ENDIF(LOG)
## ----------------------------------------------------------------------------

# CUDA_INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/jama/include)

CUDA_INCLUDE_DIRECTORIES(${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/../src
                         ${OPENMM_BUILD_FREE_ENERGY_PATH}/platforms/cuda/include
                         ${OPENMM_DIR}/platforms/cuda/src
                         ${OPENMM_DIR}/platforms/cuda/include
                         ${OPENMM_DIR}/platforms/cuda/src/kernels
                         ${OPENMM_DIR}/openmmapi/include ) 

CUDA_ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

TARGET_LINK_LIBRARIES(${SHARED_TARGET} debug ${OPENMM_LIBRARY_NAME}_d optimized ${OPENMM_LIBRARY_NAME} )
TARGET_LINK_LIBRARIES(${SHARED_TARGET} debug ${OPENMM_LIBRARY_NAME}Cuda_d optimized ${OPENMM_LIBRARY_NAME}Cuda )
TARGET_LINK_LIBRARIES(${SHARED_TARGET} debug ${OPENMM_FREE_ENERGY_LIBRARY_NAME}_d optimized ${OPENMM_FREE_ENERGY_LIBRARY_NAME} )
SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES COMPILE_FLAGS "-DOPENMMCUDAFREEENERGY_BUILDING_SHARED_LIBRARY")

INSTALL_TARGETS(/lib/plugins RUNTIME_DIRECTORY /lib/plugins ${SHARED_TARGET})
