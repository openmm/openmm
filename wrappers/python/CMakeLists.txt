if(POLICY CMP0026)
    cmake_policy(SET CMP0026 OLD)
endif()
#############################################
### Copy all source files to staging area ###
#############################################

# OPENMM_PYTHON_STAGING_DIR is a staging area for python, swig, and C files in the python package we are making.
set(OPENMM_PYTHON_STAGING_DIR "${CMAKE_BINARY_DIR}/python"
    CACHE PATH "Temporary staging area for Python API wrappers")
mark_as_advanced(OPENMM_PYTHON_STAGING_DIR)

# Create package directory structure
file(MAKE_DIRECTORY ${OPENMM_PYTHON_STAGING_DIR}/simtk/openmm)
file(MAKE_DIRECTORY ${OPENMM_PYTHON_STAGING_DIR}/simtk/unit)
file(MAKE_DIRECTORY ${OPENMM_PYTHON_STAGING_DIR}/src/swig_doxygen/swig_lib/python)
file(MAKE_DIRECTORY ${OPENMM_PYTHON_STAGING_DIR}/tests)

##############################################################################
### Identify files that need to be copied from source area to staging area ###
##############################################################################
set(STAGING_OUTPUT_FILES "") # Will contain all required package files
file(GLOB STAGING_INPUT_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/README.txt"
)

###############################################################################
# Collect the necessary information for building and compiling the Python
# package:
#  - WRAPPER_INCLUDE_DIRS: All folders containing header files that the
#    wrappers must be compiled against.
#  - ALL_PLUGIN_TARGET_LOCATIONS: All of the plugin libraries that should be
#    copied into the Python installation
#  - ALL_LIBRARIES: The OpenMM libraries that the wrappers should be linked
#    against
#  - ALL_LIBRARY_LOCATIONS: The paths to these libraries
###############################################################################

set(WRAPPER_BASE_SUBDIRS . openmmapi olla serialization plugins/amoeba/openmmapi plugins/rpmd/openmmapi plugins/drude/openmmapi)
set(WRAPPER_INCLUDE_DIRS) # start empty
foreach(subdir ${WRAPPER_BASE_SUBDIRS})
    # append
    SET(WRAPPER_INCLUDE_DIRS ${WRAPPER_INCLUDE_DIRS}
                         ${CMAKE_SOURCE_DIR}/${subdir}/include
                         ${CMAKE_SOURCE_DIR}/${subdir}/include/openmm
                         ${CMAKE_SOURCE_DIR}/${subdir}/include/openmm/internal)
endforeach(subdir)

get_property(ALL_PLUGIN_TARGETS GLOBAL PROPERTY AllPluginTargets)
set(ALL_PLUGIN_TARGET_LOCATIONS)
foreach(PLUGIN_TARGET ${ALL_PLUGIN_TARGETS})
    get_property(PLUGIN_TARGET_LOCATION TARGET ${PLUGIN_TARGET} PROPERTY LOCATION)
    set(ALL_PLUGIN_TARGET_LOCATIONS "${ALL_PLUGIN_TARGET_LOCATIONS}" "${PLUGIN_TARGET_LOCATION}")
endforeach()

set(ALL_LIBRARIES "OpenMM" "OpenMMAmoeba" "OpenMMRPMD" "OpenMMDrude")
set(ALL_LIBRARY_LOCATIONS)
foreach(LIBRARY_TARGET ${ALL_LIBRARIES})
    get_property(TARGET_LOCATION TARGET ${LIBRARY_TARGET} PROPERTY LOCATION)
    set(ALL_LIBRARY_LOCATIONS "${ALL_LIBRARY_LOCATIONS}" "${TARGET_LOCATION}")
endforeach()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py
    ${OPENMM_PYTHON_STAGING_DIR}/setup.py @ONLY)

###############################################################################
# Collect the current git revision and whether or not we are on a git tag
# (used to indicate releases), and use this information to write
# ${OPENMM_PYTHON_STAGING_DIR}/simtk/openmm/version.py
###############################################################################

execute_process(
  COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_REVISION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

execute_process(
  COMMAND git describe --tags --exact-match
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_TAG_EXACT_MATCH_STR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

if(NOT GIT_REVISION)
  set(GIT_REVISION "Unknown")
endif()
if (NOT GIT_TAG_EXACT_MATCH_STR)
    SET(GIT_TAG_EXACT_MATCH "False")
else()
    SET(GIT_TAG_EXACT_MATCH "True")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.py.in
    ${OPENMM_PYTHON_STAGING_DIR}/simtk/openmm/version.py @ONLY)


set(SUBDIRS src simtk tests)
foreach(SUBDIR ${SUBDIRS})
    file(GLOB_RECURSE STAGING_INPUT_FILES1 RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*README.txt"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.i"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.ini"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.sh"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.xml"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.pdb"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.pdbx"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.prmtop"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.prm"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.inpcrd"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.crd"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.gro"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.parm7"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.rst7"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.ncrst"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.dms"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.top"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.par"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*.str"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/*psf"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/charmm22.*"
        "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/par*.inp"
    )
    foreach(file ${STAGING_INPUT_FILES1})
        set(STAGING_INPUT_FILES ${STAGING_INPUT_FILES} "${file}")
    endforeach(file ${STAGING_INPUT_FILES1})
endforeach(SUBDIR ${SUBDIRS})

# message("STAGING_INPUT_FILES = " ${STAGING_INPUT_FILES})
foreach(INIT_FILE ${STAGING_INPUT_FILES})
    set(infile "${CMAKE_CURRENT_SOURCE_DIR}/${INIT_FILE}")
    set(outfile "${OPENMM_PYTHON_STAGING_DIR}/${INIT_FILE}")
    add_custom_command(
        OUTPUT "${outfile}"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${infile}" "${outfile}"
        DEPENDS "${infile}"
        COMMENT "CMake-copying file ${infile} to ${outfile}")
    set(STAGING_OUTPUT_FILES ${STAGING_OUTPUT_FILES} "${outfile}")
    # message("STAGING_OUTPUT_FILE = " ${outfile})
endforeach(INIT_FILE ${STAGING_INPUT_FILES})

########################################################
### Rebuild the Swig wrappers if the API has changed ###
########################################################
# Look in ~/bin/swig first, more likely to be latest version on build machines
find_program(SWIG_EXECUTABLE swig PATHS
    "$ENV{HOME}/bin" "/Users/builder/bin" "/home/builder/bin"
    NO_DEFAULT_PATH)
find_program(SWIG_EXECUTABLE swig)
MARK_AS_ADVANCED(SWIG_VERSION)
if(SWIG_EXECUTABLE)
    execute_process(COMMAND ${SWIG_EXECUTABLE} -version
        OUTPUT_VARIABLE SWIG_VERSION_output ERROR_VARIABLE SWIG_VERSION_output)
    string(REGEX REPLACE ".*SWIG Version[^0-9.]*\([0-9.]+\).*" "\\1"
        SWIG_VERSION_output "${SWIG_VERSION_output}")
    set(SWIG_VERSION ${SWIG_VERSION_output} CACHE STRING "Swig version" FORCE)
else(SWIG_EXECUTABLE)
    set(SWIG_VERSION "0.0.0" CACHE STRING "Swig version" FORCE)
endif(SWIG_EXECUTABLE)
# Enforce swig version
if (PYTHON_VERSION_MAJOR EQUAL 3)
    string(COMPARE LESS "${SWIG_VERSION}" "2.0.4" SWIG_VERSION_ERROR)
    if(SWIG_VERSION_ERROR)
        message("Swig version must be 2.0.4 or greater for Python 3! (You have ${SWIG_VERSION})")
    endif(SWIG_VERSION_ERROR)
else(PYTHON_VERSION_MAJOR EQUAL 3)
    string(COMPARE LESS "${SWIG_VERSION}" "1.3.39" SWIG_VERSION_ERROR)
    if(SWIG_VERSION_ERROR)
        message("Swig version must be 1.3.39 or greater for Python 2! (You have ${SWIG_VERSION})")
    endif(SWIG_VERSION_ERROR)
endif(PYTHON_VERSION_MAJOR EQUAL 3)

find_package(Doxygen REQUIRED)
mark_as_advanced(CLEAR DOXYGEN_EXECUTABLE)

# SWIG_OPENMM_DIR is package area where swig files will be created/used
set(SWIG_OPENMM_DIR "${OPENMM_PYTHON_STAGING_DIR}/src/swig_doxygen")
file(MAKE_DIRECTORY ${SWIG_OPENMM_DIR}/swig_lib/python)

### Create OpenMMSwigHeaders.i ###
# Step 1 - Create Doxyfile to point to OpenMM headers
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/swig_doxygen/doxygen/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/src/swig_doxygen/doxygen/Doxyfile
    @ONLY
)

# Step 2 - Run doxygen in non-package area
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/src/swig_doxygen/doxygen/xml/index.xml"
    COMMAND "${DOXYGEN_EXECUTABLE}"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/src/swig_doxygen/doxygen/Doxyfile"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/src/swig_doxygen/doxygen"
    COMMENT "Parsing OpenMM header files with doxygen..."
)

# OpenMM REFERENCE platform first

# Files which will be created by python script and used by swig (dynamically generated)
set(SWIG_INPUT_FILES_REF
    "${SWIG_OPENMM_DIR}/swig_lib/python/pythonprepend.i"
    "${SWIG_OPENMM_DIR}/swig_lib/python/pythonappend.i"
    "${SWIG_OPENMM_DIR}/OpenMM_docstring.i"
    "${SWIG_OPENMM_DIR}/OpenMMSwigHeaders.i"
    # "${SWIG_OPENMM_DIR}/OpenMM.i"
)

# Swig input files which should have been copied from source tree
# file(GLOB SWIG_INPUT_FILES2 "${SWIG_OPENMM_DIR}/swig_lib/python/*.i")
set(SWIG_INPUT_FILES2
    "${SWIG_OPENMM_DIR}/swig_lib/python/exceptions.i"
    "${SWIG_OPENMM_DIR}/swig_lib/python/extend.i"
    "${SWIG_OPENMM_DIR}/swig_lib/python/features.i"
    "${SWIG_OPENMM_DIR}/swig_lib/python/header.i"
    "${SWIG_OPENMM_DIR}/swig_lib/python/pythoncode.i"
    "${SWIG_OPENMM_DIR}/swig_lib/python/typemaps.i"
)

# Create input files for swig
add_custom_command(
    OUTPUT ${SWIG_INPUT_FILES_REF}
    COMMAND ${PYTHON_EXECUTABLE} "${SWIG_OPENMM_DIR}/swigInputBuilder.py"
        -i "${CMAKE_CURRENT_BINARY_DIR}/src/swig_doxygen/doxygen/xml"
        -c swigInputConfig.py
        -d OpenMM_docstring.i
        -o OpenMMSwigHeaders.i
        -a swig_lib/python/pythonprepend.i
        -z swig_lib/python/pythonappend.i
        -v "${SWIG_VERSION}"
    WORKING_DIRECTORY "${SWIG_OPENMM_DIR}"
    DEPENDS
        "${SWIG_OPENMM_DIR}/swigInputConfig.py"
        "${SWIG_OPENMM_DIR}/swigInputBuilder.py"
        "${CMAKE_CURRENT_BINARY_DIR}/src/swig_doxygen/doxygen/xml/index.xml"
    COMMENT "Creating OpenMM Python swig input files..."
)

add_custom_target(OpenMMSwigHeaders ALL DEPENDS "${SWIG_OPENMM_DIR}/OpenMMSwigHeaders.i")
INSTALL_FILES(/include FILES "${SWIG_OPENMM_DIR}/OpenMMSwigHeaders.i")

#~ swig -python -c++ \
#~      -outdir $PYTHON_PACKAGE_DIR \
#~      -o OpenMMSwig.cxx \
#~      OpenMM.i

# Run swig
add_custom_command(
    OUTPUT "${SWIG_OPENMM_DIR}/OpenMMSwig.cxx" "${OPENMM_PYTHON_STAGING_DIR}/simtk/openmm/openmm.py"
    COMMAND ${SWIG_EXECUTABLE}
        -python -c++
        -outdir "${OPENMM_PYTHON_STAGING_DIR}/simtk/openmm"
        -o OpenMMSwig.cxx
        OpenMM.i
    WORKING_DIRECTORY "${SWIG_OPENMM_DIR}"
    DEPENDS
        "${SWIG_OPENMM_DIR}/OpenMM.i"
        ${SWIG_INPUT_FILES_REF}
        ${SWIG_INPUT_FILES2}
    COMMENT "Creating OpenMM Python module sources with swig..."
)
add_custom_target(RunSwig DEPENDS
    "${SWIG_OPENMM_DIR}/OpenMMSwig.cxx"
    "${OPENMM_PYTHON_STAGING_DIR}/simtk/openmm/openmm.py")

set(STAGING_OUTPUT_FILES ${STAGING_OUTPUT_FILES}
    "${OPENMM_PYTHON_STAGING_DIR}/src/swig_doxygen/OpenMMSwig.cxx"
    "${OPENMM_PYTHON_STAGING_DIR}/simtk/openmm/openmm.py"
    "${OPENMM_PYTHON_STAGING_DIR}/simtk/openmm/version.py"
    "${OPENMM_PYTHON_STAGING_DIR}/setup.py"
)


###########################################################################
### Run python setup.py
###########################################################################

# Binary distribution
add_custom_target(
    PythonBdist
    COMMAND "${PYTHON_EXECUTABLE}" setup.py bdist_wheel
    DEPENDS ${STAGING_OUTPUT_FILES} ${ALL_PLUGIN_TARGETS}
    WORKING_DIRECTORY ${OPENMM_PYTHON_STAGING_DIR}
    COMMENT "Installing OpenMM Python binary module..."
)


# Install Python package (to system location)
add_custom_target(
    PythonInstall
    COMMAND "${PYTHON_EXECUTABLE}" setup.py install
    # COMMAND "touch ${CMAKE_BINARY_DIR}/.PythonInstall"
    DEPENDS ${STAGING_OUTPUT_FILES} ${ALL_PLUGIN_TARGETS}
    WORKING_DIRECTORY ${OPENMM_PYTHON_STAGING_DIR}
    COMMENT "Installing OpenMM Python binary module..."
)

install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build \"${CMAKE_BINARY_DIR}\" --target PythonInstall)")
