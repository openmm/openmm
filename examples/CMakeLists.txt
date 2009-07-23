# Generate and install examples.
#
# This is boilerplate code for generating a set of executables, one per
# .cpp file in an "examples" subdirectory. 
#
# For IDEs that can deal with PROJECT_LABEL properties (at least
# Visual Studio) the projects for building each of these adhoc
# executables will be labeled "Example - TheExampleName" if a file
# TheExampleName.cpp is found in this directory.
#
# We check the BUILD_TESTING_{SHARED,STATIC} variables to determine
# whether to build dynamically linked, statically linked, or both
# versions of the executable.


SET(OpenMM_CWRAPPER "OpenMMCWrapper")
SET(OpenMM_FWRAPPER "OpenMMFortranWrapper")
SET(OpenMM_FMODULE  "OpenMMFortranModule")

# CWrapper is always a static library but you need a different
# one if you want to link with the OpenMM static library.
SET(CWRAPPER_FOR_SHARED ${OpenMM_CWRAPPER})
SET(CWRAPPER_FOR_STATIC ${OpenMM_CWRAPPER}_static)

# Visual Studio will add "_d" to target names automatically
# but on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(CWRAPPER_FOR_SHARED ${CWRAPPER_FOR_SHARED}_d)
    SET(CWRAPPER_FOR_STATIC ${CWRAPPER_FOR_STATIC}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

SET(CPP_EXAMPLES HelloArgon HelloSodiumChloride HelloEthane HelloWaterBox)
SET(C_EXAMPLES HelloArgonInC HelloSodiumChlorideInC)
SET(F_EXAMPLES HelloArgonInFortran HelloSodiumChlorideInFortran)

SET(BUILD_TESTING_SHARED 1)
SET(BUILD_TESTING_STATIC 1)

ADD_LIBRARY(${CWRAPPER_FOR_SHARED} 
		    ${OpenMM_CWRAPPER}.cpp ${OpenMM_FWRAPPER}.cpp ${OpenMM_CWRAPPER}.h)
SET_TARGET_PROPERTIES(${CWRAPPER_FOR_SHARED}
		PROPERTIES
		  PROJECT_LABEL "C Bindings")
TARGET_LINK_LIBRARIES(${CWRAPPER_FOR_SHARED} ${SHARED_TARGET})
		
ADD_LIBRARY(${CWRAPPER_FOR_STATIC}
		    ${OpenMM_CWRAPPER}.cpp ${OpenMM_FWRAPPER}.cpp ${OpenMM_CWRAPPER}.h)
SET_TARGET_PROPERTIES(${CWRAPPER_FOR_STATIC}
		PROPERTIES
		  COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
		  PROJECT_LABEL "C Bindings for static OpenMM")
		  
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${CWRAPPER_FOR_SHARED})
INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${CWRAPPER_FOR_STATIC})


FOREACH(EX_ROOT ${CPP_EXAMPLES})
    IF (BUILD_TESTING_SHARED)
        # Link with shared library
        ADD_EXECUTABLE(${EX_ROOT} ${EX_ROOT}.cpp)
        SET_TARGET_PROPERTIES(${EX_ROOT}
		PROPERTIES
	      PROJECT_LABEL "Example - ${EX_ROOT}")
        TARGET_LINK_LIBRARIES(${EX_ROOT} ${SHARED_TARGET})
    ENDIF (BUILD_TESTING_SHARED)

    IF (BUILD_TESTING_STATIC)
        # Link with static library
        SET(EX_STATIC ${EX_ROOT}Static)
        ADD_EXECUTABLE(${EX_STATIC} ${EX_ROOT}.cpp)
        SET_TARGET_PROPERTIES(${EX_STATIC}
		PROPERTIES
		COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
		PROJECT_LABEL "Example - ${EX_STATIC}")
        TARGET_LINK_LIBRARIES(${EX_STATIC} ${STATIC_TARGET})
    ENDIF (BUILD_TESTING_STATIC)

    INSTALL(FILES ${EX_ROOT}.cpp DESTINATION examples)

ENDFOREACH(EX_ROOT ${CPP_EXAMPLES})

FOREACH(EX_ROOT ${C_EXAMPLES})
    IF (BUILD_TESTING_SHARED)
        # Link with shared library
		# We need at least one .cpp here to get CMake to include
		# C++ libraries on the link line.
        ADD_EXECUTABLE(${EX_ROOT} ${EX_ROOT}.c Empty.cpp)
        SET_TARGET_PROPERTIES(${EX_ROOT}
		PROPERTIES
	      PROJECT_LABEL "Example C - ${EX_ROOT}")
        TARGET_LINK_LIBRARIES(${EX_ROOT} ${CWRAPPER_FOR_SHARED} ${SHARED_TARGET})
    ENDIF (BUILD_TESTING_SHARED)

    IF (BUILD_TESTING_STATIC)
        # Link with static library
        SET(EX_STATIC ${EX_ROOT}Static)
		# We need at least one .cpp here to get CMake to include
		# C++ libraries on the static link line.
        ADD_EXECUTABLE(${EX_STATIC} ${EX_ROOT}.c Empty.cpp)
        SET_TARGET_PROPERTIES(${EX_STATIC}
		PROPERTIES
		COMPILE_FLAGS "-DOPENMM_USE_STATIC_LIBRARIES"
		PROJECT_LABEL "Example C - ${EX_STATIC}")
        TARGET_LINK_LIBRARIES(${EX_STATIC} ${CWRAPPER_FOR_STATIC} ${STATIC_TARGET})
    ENDIF (BUILD_TESTING_STATIC)

    INSTALL(FILES ${EX_ROOT}.c DESTINATION examples)

ENDFOREACH(EX_ROOT ${C_EXAMPLES})

FOREACH(EX_ROOT ${F_EXAMPLES})
    INSTALL(FILES ${EX_ROOT}.f90 DESTINATION examples)
ENDFOREACH(EX_ROOT ${F_EXAMPLES})

INSTALL(FILES ${OpenMM_CWRAPPER}.h ${OpenMM_CWRAPPER}.cpp ${OpenMM_FWRAPPER}.cpp
		DESTINATION examples)
INSTALL(FILES ${OpenMM_FMODULE}.f90 DESTINATION examples)

INSTALL(FILES README.txt DESTINATION examples)
INSTALL(FILES Makefile NMakefile DESTINATION examples)

INSTALL(FILES HelloArgonVS8/HelloArgon.vcproj 
		      HelloArgonVS8/HelloArgon.sln 
	    DESTINATION examples/VisualStudio2005)
INSTALL(FILES HelloArgonInCVS8/HelloArgonInC.vcproj 
		      HelloArgonInCVS8/HelloArgonInC.sln 
	    DESTINATION examples/VisualStudio2005)
INSTALL(FILES HelloOpenMM.ppt 
        DESTINATION examples)
