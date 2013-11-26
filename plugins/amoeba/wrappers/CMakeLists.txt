# WRAPPER_DOXYGEN_DIR is a workspace directory where wrapper files will be created
set(WRAPPER_DOXYGEN_DIR "${CMAKE_CURRENT_BINARY_DIR}/doxygen")
file(MAKE_DIRECTORY "${WRAPPER_DOXYGEN_DIR}")

# Step 1 - Create Doxyfile to point to OpenMM headers
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${WRAPPER_DOXYGEN_DIR}/Doxyfile
    @ONLY
)

# Step 2 - Run doxygen to analyze the headers
add_custom_command(
    OUTPUT "${WRAPPER_DOXYGEN_DIR}/xml/index.xml"
    COMMAND "${DOXYGEN_EXECUTABLE}"
    DEPENDS "${WRAPPER_DOXYGEN_DIR}/Doxyfile"
    WORKING_DIRECTORY "${WRAPPER_DOXYGEN_DIR}"
    COMMENT "Parsing OpenMM header files with Doxygen..."
)

# Step 3 - Generate the wrappers
ADD_CUSTOM_COMMAND(OUTPUT AmoebaOpenMMCWrapper.h AmoebaOpenMMCWrapper.cpp AmoebaOpenMMFortranModule.f90 AmoebaOpenMMFortranWrapper.cpp
     COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/generateAmoebaWrappers.py" "${WRAPPER_DOXYGEN_DIR}/xml" "${CMAKE_CURRENT_BINARY_DIR}"
     DEPENDS "${WRAPPER_DOXYGEN_DIR}/xml/index.xml" "${CMAKE_CURRENT_SOURCE_DIR}/generateAmoebaWrappers.py")

ADD_CUSTOM_TARGET(AmoebaApiWrappers DEPENDS AmoebaOpenMMCWrapper.h AmoebaOpenMMCWrapper.cpp AmoebaOpenMMFortranModule.f90 AmoebaOpenMMFortranWrapper.cpp)

INSTALL_FILES(/include FILES AmoebaOpenMMCWrapper.h AmoebaOpenMMFortranModule.f90)
