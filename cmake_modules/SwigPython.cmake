# OPENMM_PYTHON_STAGING_DIR is a staging area for python, swig, and C files in the python package we are making.
set(OPENMM_PYTHON_STAGING_DIR "${CMAKE_BINARY_DIR}/python" CACHE PATH "Temporary staging area for Python API wrappers")
mark_as_advanced(OPENMM_PYTHON_STAGING_DIR)

# Create the directory where configuration files will go
set(SWIG_CONFIG_DIRECTORY ${OPENMM_PYTHON_STAGING_DIR}/src/swig_doxygen/config)
file(MAKE_DIRECTORY ${SWIG_CONFIG_DIRECTORY})

# Each plugin may define its own SWIG configuration file for building the Python API.  It
# should call ADD_SWIG_CONFIG_FILE() to add its file to the list and copy the file to the
# staging directory.
add_custom_target(CopySwigConfigFiles)
set(OPENMM_SWIG_CONFIG_FILES "")
macro(ADD_SWIG_CONFIG_FILE FILE)
    set(OPENMM_SWIG_CONFIG_FILES ${OPENMM_SWIG_CONFIG_FILES} ${FILE})
    set(OPENMM_SWIG_CONFIG_FILES ${OPENMM_SWIG_CONFIG_FILES} PARENT_SCOPE)
    file(COPY ${FILE} DESTINATION ${SWIG_CONFIG_DIRECTORY})
endmacro(ADD_SWIG_CONFIG_FILE FILE)
