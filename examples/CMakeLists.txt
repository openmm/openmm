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

SET(CPP_EXAMPLES HelloArgon HelloSodiumChloride HelloEthane HelloWaterBox)
SET(C_EXAMPLES HelloArgonInC HelloSodiumChlorideInC)
SET(F_EXAMPLES HelloArgonInFortran HelloSodiumChlorideInFortran)

SET(BUILD_TESTING_SHARED 1)
SET(BUILD_TESTING_STATIC OFF)
IF(OPENMM_BUILD_STATIC_LIB)
    SET(BUILD_TESTING_STATIC ON)
ENDIF(OPENMM_BUILD_STATIC_LIB)
		
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

# Only build wrapper examples if wrappers have been built
IF(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS)
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_BINARY_DIR}/wrappers)

    FOREACH(EX_ROOT ${C_EXAMPLES})
        IF (BUILD_TESTING_SHARED)
            # Link with shared library
            # We need at least one .cpp here to get CMake to include
            # C++ libraries on the link line.
            ADD_EXECUTABLE(${EX_ROOT} ${EX_ROOT}.c Empty.cpp)
            SET_TARGET_PROPERTIES(${EX_ROOT}
            PROPERTIES
              PROJECT_LABEL "Example C - ${EX_ROOT}")
            TARGET_LINK_LIBRARIES(${EX_ROOT} ${SHARED_TARGET})
            ADD_DEPENDENCIES(${EX_ROOT} ApiWrappers)
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
            TARGET_LINK_LIBRARIES(${EX_STATIC} ${STATIC_TARGET})
            ADD_DEPENDENCIES(${EX_STATIC} ApiWrappers)
        ENDIF (BUILD_TESTING_STATIC)
  ENDFOREACH(EX_ROOT ${C_EXAMPLES})
ENDIF(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS)

FOREACH(EX_ROOT ${C_EXAMPLES})
  INSTALL(FILES ${EX_ROOT}.c DESTINATION examples)
ENDFOREACH(EX_ROOT ${C_EXAMPLES})

FOREACH(EX_ROOT ${F_EXAMPLES})
    INSTALL(FILES ${EX_ROOT}.f90 DESTINATION examples)
ENDFOREACH(EX_ROOT ${F_EXAMPLES})

INSTALL(FILES simulateAmber.py simulatePdb.py simulateGromacs.py testInstallation.py argon-chemical-potential.py input.inpcrd input.prmtop input.pdb input.gro input.top
        DESTINATION examples)

INSTALL(FILES VisualStudio/HelloArgon.vcproj 
              VisualStudio/HelloArgon.sln 
              VisualStudio/HelloArgonInC.sln 
              VisualStudio/HelloArgonInC.vcproj 
              VisualStudio/HelloArgonInFortran.sln 
              VisualStudio/HelloArgonInFortran.vfproj 
        DESTINATION examples/VisualStudio)

INSTALL(FILES README.txt DESTINATION examples)
INSTALL(FILES Makefile NMakefile DESTINATION examples)

INSTALL(FILES MakefileNotes.txt Empty.cpp DESTINATION examples)
