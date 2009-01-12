# - Look for the BrookGPU streaming extension to C language 
# 
# BROOK_FILE :   .BR -> .CPP 
# BROOK_INCLUDE_DIR :   Include directory for Brook.hpp 
# BROOK_C[XX]FLAGS :   Flags needed to compile the produced CPP file 
# BROOK_LIBRARIES :   All needed libraries, including OpenGL for OGL backend 
# 
# Used internally : 
# BROOK_CC :   Location of BRCC 
# BROOK_xxx_LIBRARY :   Location of the various libraries used by brook 
# 

# ----------------------------------------------------------------------------

FIND_PATH(BROOK_INCLUDE_DIR brook 
   $ENV{BROOKDIR}/include 
   $ENV{BROOKROOT}/sdk/include 
   /usr/include/ 
   /usr/local/include/  
   ~/src/cvs/brook/include 
   ~/src/brook/include 
   ~/brook/include) 

SET(BROOK_CXXFLAGS "-I${BROOK_INCLUDE_DIR}") 
SET(BROOK_CFLAGS "${BROOK_CXXFLAGS}") 

FIND_PROGRAM(BROOK_CC brcc 
   $ENV{BROOKDIR}/bin 
   $ENV{BROOKROOT}/sdk/bin 
   /usr/bin/ 
   /usr/local/bin/ 
   ~/src/cvs/brook/bin 
   ~/src/brook/bin 
   ~/brook/bin) 

# Search for all libraries 
# - both BASE and RUNTIME TARGETS 

# ----------------------------------------------------------------------------

FIND_LIBRARY(BROOK_brook_LIBRARY 
   NAMES  
      brook
      brook_d
   PATHS 
      $ENV{BROOKDIR}/lib 
      $ENV{BROOKDIR}/bin 
      $ENV{BROOKROOT}/sdk/lib 
      /usr/lib 
      /usr/local/lib 
      ~/src/brook/bin 
      ~/brook/bin) 

# if found, add to list 

IF (BROOK_brook_LIBRARY) 
   SET(BROOK_LIBRARIES ${BROOK_LIBRARIES} ${BROOK_brook_LIBRARY}) 
ENDIF (BROOK_brook_LIBRARY) 

# all individual libs are advanced settings 

MARK_AS_ADVANCED(BROOK_brook_LIBRARY) 

# ----------------------------------------------------------------------------

IF(LOG)
   FILE( APPEND ${LOG_FILE} "\nIn FindBrook.cmake\n" )
   FILE( APPEND ${LOG_FILE} "BROOK_INCLUDE_DIR=${BROOK_INCLUDE_DIR}\n" )
   FILE( APPEND ${LOG_FILE} "BROOK_CC=${BROOK_CC}\n" )
   FILE( APPEND ${LOG_FILE} "BROOK_brook_LIBRARY=${BROOK_brook_LIBRARY}\n" )
   FILE( APPEND ${LOG_FILE} "BROOKROOT=<$ENV{BROOKROOT}>\n" )
   FILE( APPEND ${LOG_FILE} "sub_lib=<${sub_lib}>\n" )
ENDIF(LOG)

# ----------------------------------------------------------------------------

# check if includes and main lib are here 

IF (BROOK_INCLUDE_DIR AND BROOK_brook_LIBRARY AND BROOK_CC) 

   SET(BROOK_FOUND TRUE) 

   GET_FILENAME_COMPONENT(BROOK_LIB_PATH "${BROOK_brook_LIBRARY}" PATH) 
   # ----------------------------------------------------------------------------
   IF(LOG)
      FILE( APPEND ${LOG_FILE} "BROOK_LIB_PATH=${BROOK_LIB_PATH}\n" )
   ENDIF(LOG)
   # ----------------------------------------------------------------------------

   # Implementation to allow interpreting/compiling Brook files 

   MACRO(BROOK_FILE FILENAME) 

      IF(LOG)
         FILE( APPEND ${LOG_FILE} "1 In BROOK_FILE: ${FILENAME}\n" )
      ENDIF(LOG)

      # split input names 

      GET_FILENAME_COMPONENT(PATH "${FILENAME}" PATH) 
      GET_FILENAME_COMPONENT(HEAD "${FILENAME}" NAME_WE) # without trailing ".BR" 

      # File names 

      SET(OUTPATH "${CMAKE_CURRENT_BINARY_DIR}/src/gpu") 
      SET(BROOK_PREFIX "${OUTPATH}/${HEAD}") 
      SET(OUTFILE "${BROOK_PREFIX}.cpp") # file produced by Brook 
      # SET(INFILE "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}") # canonical input name 
      SET(INFILE "${FILENAME}") # canonical input name 

      # ----------------------------------------------------------------------------
      IF(LOG)
         FILE( APPEND ${LOG_FILE} "In BROOK_FILE: ${FILENAME} FILENAME\n" )
         FILE( APPEND ${LOG_FILE} "   Path=${PATH} HEAD=${HEAD}\n" )
         FILE( APPEND ${LOG_FILE} "   OUTPATH=${OUTPATH}\n   BROOK_PREFIX=${BROOK_PREFIX}\n" )
         FILE( APPEND ${LOG_FILE} "   OUTFILE=${OUTFILE}\n   INFILE=${INFILE}\n\n" )
      ENDIF(LOG)
      # ----------------------------------------------------------------------------

      # create output path, if it does not exist 

      IF(NOT EXISTS "${OUTPATH}") 
         FILE(MAKE_DIRECTORY "${OUTPATH}") 
      ENDIF(NOT EXISTS "${OUTPATH}") 

      # Run Brook 

      ADD_CUSTOM_COMMAND( 
         OUTPUT    "${OUTFILE}" 
         COMMAND   "${BROOK_CC}" 
         ARGS      "-o${BROOK_PREFIX}" 
                   "${INFILE}" 
         DEPENDS   "${INFILE}" ) 

      # Flag file as generated 

      SET_SOURCE_FILES_PROPERTIES("${OUTFILE}" PROPERTIES GENERATED TRUE) 

      # accumulate Brook cpp files

      SET(BROOK_CPP_FILES ${BROOK_CPP_FILES} ${OUTFILE})

   ENDMACRO(BROOK_FILE) 
ELSE (BROOK_INCLUDE_DIR AND BROOK_brook_LIBRARY AND BROOK_CC) 
   SET(BROOK_FOUND FALSE) 
ENDIF (BROOK_INCLUDE_DIR AND BROOK_brook_LIBRARY AND BROOK_CC) 

# Some verbosity 

IF (NOT BROOK_FOUND) 
   IF (BROOK_FIND_REQUIRED) 
      MESSAGE(FATAL_ERROR "Could not find BROOK" ) 
   ENDIF (BROOK_FIND_REQUIRED) 
ENDIF (NOT BROOK_FOUND) 

# ----------------------------------------------------------------------------
IF(LOG)
   FILE( APPEND ${LOG_FILE} "BROOK_FOUND=<${BROOK_FOUND}>\n" )
   FILE( APPEND ${LOG_FILE} "\nLeaving FindBrook.cmake\n" )
ENDIF(LOG)
# ----------------------------------------------------------------------------
