#
# Include REFERENCE related files.
#
# ----------------------------------------------------------------------------
# logging

SET(LOG FALSE)
   
IF(LOG)
   SET(LOG_FILE "CMakeLog.txt" )
   FILE( WRITE ${LOG_FILE} "In amoeba/platforms/reference/sharedTarget Cmake\n")
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

SET(OPENMM_BUILD_AMOEBA_PATH ${CMAKE_SOURCE_DIR}/plugins/amoeba) 

# ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "Pre OPENMM_SOURCE_SUBDIRS" ${OPENMM_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "Pre OPENMM_AMOEBA_SOURCE_SUBDIRS " ${OPENMM_AMOEBA_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "Pre SOURCE_FILES" ${SOURCE_FILES} )
ENDIF(LOG)
## ----------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${REFERENCE_INCLUDE})
LINK_DIRECTORIES(${REFERENCE_TARGET_LINK})
FOREACH(subdir ${OPENMM_AMOEBA_SOURCE_SUBDIRS})
    FILE(GLOB src_files ${OPENMM_BUILD_AMOEBA_PATH}/platforms/reference/${subdir}/src/*.cu ${OPENMM_BUILD_AMOEBA_PATH}/platforms/reference/src/*/*.cu)
    SET(SOURCE_FILES ${SOURCE_FILES} ${src_files})
    INCLUDE_DIRECTORIES(BEFORE ${OPENMM_BUILD_AMOEBA_PATH}/platforms/reference/../${subdir}/include)
ENDFOREACH(subdir)

# ----------------------------------------------------------------------------
IF(LOG)
   LOG_DIR( ${LOG_FILE} "OPENMM_BUILD_AMOEBA_PATH" ${OPENMM_BUILD_AMOEBA_PATH} )
   FILE( APPEND ${LOG_FILE} "OPENMM_BUILD_AMOEBA_PATH=${OPENMM_BUILD_AMOEBA_PATH}\n")
   LOG_DIR( ${LOG_FILE} "OPENMM_SOURCE_SUBDIRS" ${OPENMM_SOURCE_SUBDIRS} )
   LOG_DIR( ${LOG_FILE} "CMAKE_SOURCE_DIR" ${CMAKE_SOURCE_DIR} )
   LOG_DIR( ${LOG_FILE} "REFERENCE_INCLUDE" ${REFERENCE_INCLUDE} )
   LOG_DIR( ${LOG_FILE} "REFERENCE_TARGET_LINK" ${REFERENCE_TARGET_LINK} )
   LOG_DIR( ${LOG_FILE} "SHARED_TARGET" ${SHARED_TARGET} )
   LOG_DIR( ${LOG_FILE} "OPENMM_DIR" ${OPENMM_DIR} )
   LOG_DIR( ${LOG_FILE} "SOURCE_FILES" ${SOURCE_FILES} )
ENDIF(LOG)
## ----------------------------------------------------------------------------

# REFERENCE_INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR}/jama/include)

INCLUDE_DIRECTORIES(${OPENMM_BUILD_AMOEBA_PATH}/platforms/reference/../src
                         ${OPENMM_BUILD_AMOEBA_PATH}/platforms/reference/include
                         ${OPENMM_DIR}/platforms/reference/src
                         ${OPENMM_DIR}/platforms/reference/include
                         ${OPENMM_DIR}/platforms/reference/src/kernels
                         ${OPENMM_DIR}/openmmapi/include ) 

ADD_LIBRARY(${SHARED_TARGET} SHARED ${SOURCE_FILES} ${SOURCE_INCLUDE_FILES} ${API_ABS_INCLUDE_FILES})

TARGET_LINK_LIBRARIES(${SHARED_TARGET} debug ${OPENMM_LIBRARY_NAME}_d optimized ${OPENMM_LIBRARY_NAME} )
TARGET_LINK_LIBRARIES(${SHARED_TARGET} debug ${OPENMM_LIBRARY_NAME}Cuda_d optimized ${OPENMM_LIBRARY_NAME}Cuda )
TARGET_LINK_LIBRARIES(${SHARED_TARGET} debug ${OPENMM_AMOEBA_LIBRARY_NAME}_d optimized ${OPENMM_AMOEBA_LIBRARY_NAME} )
SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES COMPILE_FLAGS "-DOPENMMREFERENCE_BUILDING_SHARED_LIBRARY -DOPENMMREFERENCEAMOEBA_BUILDING_SHARED_LIBRARY")

INSTALL_TARGETS(/lib/plugins RUNTIME_DIRECTORY /lib/plugins ${SHARED_TARGET})
