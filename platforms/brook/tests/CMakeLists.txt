# # Testing
#

ENABLE_TESTING()

SET(BROOK_SDK ../../../brook/sdk)
SET(BROOK_SRC ../src )
SET(BROOK_LIB brook_d)

# SET(BROOK_GPU  ../../../platforms/brook/src/gpu)
SET(BROOK_GPU  C:\\cygwin\\home\\friedrim\\src\\openmm\\trunk\\OpenMM\\platforms\\brook\\src\\gpu)
SET(BROOK_GPU_LIB gmxgpu_d)

# Automatically create tests using files named "Test*.cpp"
FILE(GLOB TEST_PROGS "*Test*.cpp")
FOREACH(TEST_PROG ${TEST_PROGS})
    GET_FILENAME_COMPONENT(TEST_ROOT ${TEST_PROG} NAME_WE)

    # Link with shared library
#     ADD_EXECUTABLE(${TEST_ROOT} ${TEST_PROG})
#     LINK_DIRECTORIES( ${TEST_ROOT} ${BROOK_SDK} )
#     TARGET_LINK_LIBRARIES(${TEST_ROOT} ${SHARED_TARGET} ${BROOK_LIB} )
#     ADD_TEST(${TEST_ROOT} ${EXECUTABLE_OUTPUT_PATH}/${TEST_ROOT})


    # Link with static library
    SET(TEST_STATIC ${TEST_ROOT}Static)
    LINK_DIRECTORIES(${TEST_STATIC} ${BROOK_SDK}/lib ${BROOK_GPU})
    # STRING_APPEND( CMAKE_EXE_LINKER_FLAGS_DEBUG      "/NODEFAULTLIB:\"LIBCMT.lib\"")
    SET( CMAKE_EXE_LINKER_FLAGS_DEBUG      "/NODEFAULTLIB:\"LIBCMT.lib\"")
    ADD_EXECUTABLE(${TEST_STATIC} ${TEST_PROG})
    INCLUDE_DIRECTORIES(BEFORE ${BROOK_SRC} ${BROOK_GPU})
    ADD_DEFINITIONS(-D_WIN32 )
    SET_TARGET_PROPERTIES(${TEST_STATIC}
                 PROPERTIES
                 COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
                 )
    IF (CMAKE_BUILD_TYPE MATCHES Debug)
       TARGET_LINK_LIBRARIES(${TEST_STATIC} ${STATIC_TARGET} ${BROOK_GPU_LIB}_d ${BROOK_LIB}_d)
    ELSE(CMAKE_BUILD_TYPE MATCHES Debug)
       TARGET_LINK_LIBRARIES(${TEST_STATIC} ${STATIC_TARGET} ${BROOK_GPU_LIB} ${BROOK_LIB})
    ENDIF (CMAKE_BUILD_TYPE MATCHES Debug)

    ADD_TEST(${TEST_STATIC} ${EXECUTABLE_OUTPUT_PATH}/${TEST_STATIC})

ENDFOREACH(TEST_PROG ${TEST_PROGS})

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
# BROOK_FILE code based on : 
#   http://www.vtk.org/Wiki/CMakeUserUseBison 
# 
# CUDA has a nice CMake file too, BTW : 
# 

#FIND_PATH(BROOK_INCLUDE_DIR "$ENV{BROOKROOT}/sdk/include" ) 
SET(BROOK_INCLUDE_DIR "$ENV{BROOKROOT}\\sdk\\include" ) 

SET(BROOK_CXXFLAGS "-I${BROOK_INCLUDE_DIR}") 
SET(BROOK_CFLAGS "${BROOK_CXXFLAGS}") 

FIND_PROGRAM(BROOK_CC brcc $ENV{BROOKROOT}/sdk/bin ) 

# FIND_LIBRARY(BROOK_LIBRARIES NAMES brook PATHS ENV BROOKROOT PATH_SUFFIXES sdk/lib )
FIND_LIBRARY(BROOK_LIBRARIES brook PATHS $ENV{BROOKROOT} "$ENV{BROOKROOT}/sdk/lib" ENV BROOKROOT PATH_SUFFIXES sdk/lib )

MESSAGE(STATUS "\n********* In GpuBrook ${BROOK_CC} <${BROOK_LIBRARIES}> incl=<${BROOK_INCLUDE_DIR}>  *********\n" )

# Search for all libraries 
# - both BASE and RUNTIME TARGETS 
# FOREACH(sub_lib brook ogl gpu-dx9 gpu-ctm gpu cpu GLEE) 
#FOREACH(sub_lib brook gpu ) 
#
#   # Look for lib 
#
#   FIND_LIBRARY(BROOK_${sub_lib}_LIBRARY 
#      NAMES  
#         ${sub_lib} 
#      PATHS 
#         "$ENV{BROOKROOT}"
#         "$ENV{BROOKROOT}/sdk" 
#         "$ENV{BROOKROOT}/sdk/lib"
#      PATH_SUFFIXES
#         sdk
#         sdk/lib
#      ) 
#
#   # if found, add to list 
#   IF (BROOK_${sub_lib}_LIBRARY) 
#      SET(BROOK_LIBRARIES ${BROOK_LIBRARIES} ${BROOK_${sub_lib}_LIBRARY}) 
#   ELSE (BROOK_${sub_lib}_LIBRARY)
#      MESSAGE(STATUS "\n********* X  $ENV{BROOKROOT} <${sub_lib}> ${BROOK_${sub_lib}_LIBRARY}   *********\n" )
#   ENDIF (BROOK_${sub_lib}_LIBRARY) 
#
#   # all individual libs are advanced settings 
#   MARK_AS_ADVANCED(BROOK_${sub_lib}_LIBRARY) 
#
#ENDFOREACH(sub_lib) 

# SET(BROOK_brook_LIBRARY $ENV{BROOKROOT}/sdk/lib/brook/lib) 
# MARK_AS_ADVANCED(BROOK_brook_LIBRARY) 
# SET(BROOK_LIBRARIES ${BROOK_LIBRARIES} ${BROOK_brook_LIBRARY}) 

MESSAGE(STATUS "\n********* In GpuBrook ${BROOK_LIBRARIES}  *********\n" )

# Runtime target "OGL" requires OpenGL and GLU. 
# IF (BROOK_ogl_LIBRARY)    
#    # OpenGL module 
#    IF(NOT OPENGL_FOUND) 
#       #INCLUDE(FindOpenGL) 
#       FIND_PACKAGE(OpenGL REQUIRED) 
#    ENDIF(NOT OPENGL_FOUND) 
#    # Add all OpenGL includes (which, in turn, may add X11 includes) 
#    IF (OPENGL_FOUND) 
#            SET (BROOK_LIBRARIES ${BROOK_LIBRARIES} ${OPENGL_LIBRARIES}) 
#    ENDIF (OPENGL_FOUND) 
# 
#    # NOTE: Old CMake used to have a separate FindGLU.cmake file 
# ENDIF (BROOK_ogl_LIBRARY) 

# Runtime target "CPU" can use OpenMP where available 
#... not implemented yet... 

# Runtime targets "GPU-DX9" and "GPU-CTM" need to be implemented too... 
# GPU-DX9 : needs Microsoft DirectX9 
# GPU-CTM : needs ATI's CTM (Note: will get deprecated by Brook+/CAL) 

# check if includes and main lib are here 
IF (BROOK_INCLUDE_DIR AND BROOK_brook_LIBRARY AND BROOK_CC) 
   SET(BROOK_FOUND TRUE) 

   # Implementation to allow interpreting/compiling Brook files 
   MACRO(BROOK_FILE FILENAME) 
      # split input names 
      GET_FILENAME_COMPONENT(PATH "${FILENAME}" PATH) 
      GET_FILENAME_COMPONENT(HEAD "${FILENAME}" NAME_WE) # without trailing ".BR" 


      # File names 
      SET(OUTPATH "${CMAKE_CURRENT_BINARY_DIR}/${PATH}") 
      SET(BROOK_PREFIX "${OUTPATH}/${HEAD}") 
      SET(OUTFILE "${BROOK_PREFIX}.cpp") # file produced by Brook 
      SET(INFILE "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}") # canonical input name 

      # check output path 
      IF(NOT EXISTS "${OUTPATH}") 
         FILE(MAKE_DIRECTORY "${OUTPATH}") 
      ENDIF(NOT EXISTS "${OUTPATH}") 

      # Run Brook 
      ADD_CUSTOM_COMMAND( 
         OUTPUT   "${OUTFILE}" 
         COMMAND   "${BROOK_CC}" 
         ARGS   "-o${BROOK_PREFIX}" 
            "${INFILE}" 
         DEPENDS   "${INFILE}" ) 

      # Flag file as generated 
      SET_SOURCE_FILES_PROPERTIES("${OUTFILE}" PROPERTIES GENERATED TRUE) 

   ENDMACRO(BROOK_FILE) 
    
ENDIF (BROOK_INCLUDE_DIR AND BROOK_brook_LIBRARY AND BROOK_CC) 

# Some verbosity 
IF (BROOK_FOUND) 
   IF (NOT BROOK_FIND_QUIETLY) 
      MESSAGE(STATUS "Found BROOK: ${BROOK_CC} " ) 
   ENDIF (NOT BROOK_FIND_QUIETLY) 
ELSE (BROOK_FOUND) 
   IF (BROOK_FIND_REQUIRED) 
      MESSAGE(FATAL_ERROR "Could not find BROOK" ) 
   ENDIF (BROOK_FIND_REQUIRED) 
ENDIF (BROOK_FOUND) 

# include_directories (${HELLO_SOURCE_DIR}/Hello)
# 
# # Make sure the linker can find the Hello library once it is built.
# link_directories (${HELLO_BINARY_DIR}/Hello)
# 
# # Add executable called "helloDemo" that is built from the source files
# # "demo.cxx" and "demo_b.cxx".  The extensions are automatically found.
# add_executable (helloDemo demo.cxx demo_b.cxx)
# 
# # Link the executable to the Hello library.
# target_link_libraries (helloDemo Hello)
# 
