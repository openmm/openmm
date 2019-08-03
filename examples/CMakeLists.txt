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

FOREACH(EX_ROOT ${CPP_EXAMPLES})
    IF (OPENMM_BUILD_SHARED_LIB)
        # Link with shared library
        ADD_EXECUTABLE(${EX_ROOT} ${EX_ROOT}.cpp)
        SET_TARGET_PROPERTIES(${EX_ROOT}
            PROPERTIES
            PROJECT_LABEL "Example - ${EX_ROOT}"
            LINK_FLAGS "${EXTRA_LINK_FLAGS}"
            COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
        TARGET_LINK_LIBRARIES(${EX_ROOT} ${SHARED_TARGET})
    ENDIF (OPENMM_BUILD_SHARED_LIB)

    IF (OPENMM_BUILD_STATIC_LIB)
        # Link with static library
        SET(EX_STATIC ${EX_ROOT}Static)
        ADD_EXECUTABLE(${EX_STATIC} ${EX_ROOT}.cpp)
        SET_TARGET_PROPERTIES(${EX_STATIC}
            PROPERTIES
            PROJECT_LABEL "Example - ${EX_STATIC}"
            LINK_FLAGS "${EXTRA_LINK_FLAGS}"
            COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DOPENMM_USE_STATIC_LIBRARIES")
        TARGET_LINK_LIBRARIES(${EX_STATIC} ${STATIC_TARGET})
    ENDIF (OPENMM_BUILD_STATIC_LIB)

    INSTALL(FILES ${EX_ROOT}.cpp DESTINATION examples)

ENDFOREACH(EX_ROOT ${CPP_EXAMPLES})

# Only build wrapper examples if wrappers have been built
IF(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS)
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_BINARY_DIR}/wrappers)

    FOREACH(EX_ROOT ${C_EXAMPLES})
        IF (OPENMM_BUILD_SHARED_LIB)
            # Link with shared library
            # We need at least one .cpp here to get CMake to include
            # C++ libraries on the link line.
            ADD_EXECUTABLE(${EX_ROOT} ${EX_ROOT}.c Empty.cpp)
            SET_TARGET_PROPERTIES(${EX_ROOT}
                PROPERTIES
                PROJECT_LABEL "Example C - ${EX_ROOT}"
                LINK_FLAGS "${EXTRA_LINK_FLAGS}"
                COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
            TARGET_LINK_LIBRARIES(${EX_ROOT} ${SHARED_TARGET})
            ADD_DEPENDENCIES(${EX_ROOT} ApiWrappers)
        ENDIF (OPENMM_BUILD_SHARED_LIB)

        IF (OPENMM_BUILD_STATIC_LIB)
            # Link with static library
            SET(EX_STATIC ${EX_ROOT}Static)
            # We need at least one .cpp here to get CMake to include
            # C++ libraries on the static link line.
            ADD_EXECUTABLE(${EX_STATIC} ${EX_ROOT}.c Empty.cpp)
            SET_TARGET_PROPERTIES(${EX_STATIC}
                PROPERTIES
                PROJECT_LABEL "Example C - ${EX_STATIC}"
                LINK_FLAGS "${EXTRA_LINK_FLAGS}"
                COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DOPENMM_USE_STATIC_LIBRARIES")
            TARGET_LINK_LIBRARIES(${EX_STATIC} ${STATIC_TARGET})
            ADD_DEPENDENCIES(${EX_STATIC} ApiWrappers)
        ENDIF (OPENMM_BUILD_STATIC_LIB)
  ENDFOREACH(EX_ROOT ${C_EXAMPLES})
ENDIF(OPENMM_BUILD_C_AND_FORTRAN_WRAPPERS)

FOREACH(EX_ROOT ${C_EXAMPLES})
  INSTALL(FILES ${EX_ROOT}.c DESTINATION examples)
ENDFOREACH(EX_ROOT ${C_EXAMPLES})

FOREACH(EX_ROOT ${F_EXAMPLES})
    INSTALL(FILES ${EX_ROOT}.f90 DESTINATION examples)
ENDFOREACH(EX_ROOT ${F_EXAMPLES})

INSTALL(FILES simulateAmber.py simulatePdb.py simulateGromacs.py benchmark.py argon-chemical-potential.py input.inpcrd input.prmtop input.pdb input.gro input.top 5dfr_minimized.pdb 5dfr_solv-cube_equil.pdb apoa1.pdb
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
