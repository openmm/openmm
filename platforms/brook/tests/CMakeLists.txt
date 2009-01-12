#
# Testing
#

ENABLE_TESTING()

# ----------------------------------------------------------------------------
   
# logging

SET(LOG TRUE)
   
IF(LOG)
   SET(LOG_FILE "CMakeLog.txt" )
   FILE( WRITE ${LOG_FILE} "In Brook Test Cmake\n")
#  FILE( APPEND ${LOG_FILE} "BROOK_LIB_PATH=${BROOK_LIB_PATH}\n")
ENDIF(LOG)

# ----------------------------------------------------------------------------

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../brook-cmake/FindBrook.cmake)

SET(BROOK_LIB brook)
  
SET(OpenMM_BROOK_LIBRARY_NAME OpenMM_Brook)
  
SET(SHARED_BROOK_TARGET ${OpenMM_BROOK_LIBRARY_NAME})
SET(STATIC_BROOK_TARGET ${OpenMM_BROOK_LIBRARY_NAME}_static)

# Automatically create tests using files named "Test*.cpp"
FILE(GLOB TEST_PROGS "*Test*.cpp")

FOREACH(TEST_PROG ${TEST_PROGS})

    GET_FILENAME_COMPONENT(TEST_ROOT ${TEST_PROG} NAME_WE)

   # Link with shared library

   ADD_EXECUTABLE(${TEST_ROOT} ${TEST_PROG})
   TARGET_LINK_LIBRARIES(${TEST_ROOT} ${SHARED_TARGET})
   #ADD_TEST(${TEST_ROOT} ${EXECUTABLE_OUTPUT_PATH}/${TEST_ROOT})

   # ----------------------------------------------------------------------------
   IF(LOG)
      FILE( APPEND ${LOG_FILE} "TARGET_LINK_LIBRARIES: ${TEST_PROG} TARGET=${SHARED_TARGET} BROOK_TARGET=${ROOK_TARGET} BROOK_LIB=${BROOK_LIB}\n")
   ENDIF(LOG)
   # ----------------------------------------------------------------------------

   SET( CMAKE_EXE_LINKER_FLAGS            "/NODEFAULTLIB:\"LIBCMT.lib\"")
   SET( CMAKE_EXE_LINKER_FLAGS_DEBUG      "/NODEFAULTLIB:\"LIBCMTD.lib\"")

   ADD_DEFINITIONS( -D_WIN32 )

    # Link with static library

    SET(TEST_STATIC ${TEST_ROOT}Static)
    ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
#   SET_TARGET_PROPERTIES(${TEST_STATIC}
#                PROPERTIES
#                COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
#                )

   TARGET_LINK_LIBRARIES(${TEST_STATIC} ${STATIC_TARGET} ${STATIC_BROOK_TARGET} ${BROOK_LIB})
   # ADD_TEST(${TEST_STATIC} ${EXECUTABLE_OUTPUT_PATH}/${TEST_STATIC})

   # ----------------------------------------------------------------------------
   IF(LOG)
      FILE( APPEND ${LOG_FILE} "TARGET_LINK_LIBRARIES: ${TEST_STATIC} STATIC_TARGET=${STATIC_TARGET} STATIC_BROOK_TARGET=${STATIC_BROOK_TARGET} BROOK_LIB=${BROOK_LIB}\n")
   ENDIF(LOG)
   # ----------------------------------------------------------------------------

#  ADD_TEST(${TEST_STATIC} ${EXECUTABLE_OUTPUT_PATH}/${TEST_STATIC})

ENDFOREACH(TEST_PROG ${TEST_PROGS})

# ----------------------------------------------------------------------------

IF(LOG)
   FILE( APPEND ${LOG_FILE} "Leaving Brook Test Cmake\n")
ENDIF(LOG)

