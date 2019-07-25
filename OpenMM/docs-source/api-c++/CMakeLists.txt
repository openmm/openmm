set(STAGING_OUTPUT_FILES "") # Will contain all required package files
file(GLOB STAGING_INPUT_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "*.rst"
    "*.rst.jinja2"
    "*.py"
    "_static/logo.png"
)

set(WRAPPER_DOXYGEN_DIR "${CMAKE_CURRENT_BINARY_DIR}/doxygen")
file(MAKE_DIRECTORY "${WRAPPER_DOXYGEN_DIR}")

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${WRAPPER_DOXYGEN_DIR}/Doxyfile
    @ONLY
)

add_custom_command(
    OUTPUT "${WRAPPER_DOXYGEN_DIR}/xml/index.xml"
    COMMAND "${DOXYGEN_EXECUTABLE}"
    DEPENDS "${WRAPPER_DOXYGEN_DIR}/Doxyfile"
    WORKING_DIRECTORY "${WRAPPER_DOXYGEN_DIR}"
    COMMENT "Parsing OpenMM header files with Doxygen..."
)

foreach(INIT_FILE ${STAGING_INPUT_FILES})
    set(infile "${CMAKE_CURRENT_SOURCE_DIR}/${INIT_FILE}")
    set(outfile "${CMAKE_CURRENT_BINARY_DIR}/${INIT_FILE}")
    add_custom_command(
        OUTPUT "${outfile}"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${infile}" "${outfile}"
        DEPENDS "${infile}"
        COMMENT "CMake-copying file ${infile} to ${outfile}")
    set(STAGING_OUTPUT_FILES ${STAGING_OUTPUT_FILES} "${outfile}")
endforeach(INIT_FILE ${STAGING_INPUT_FILES})

add_custom_command(
    OUTPUT  "${CMAKE_CURRENT_BINARY_DIR}/library.rst"
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/render.py"
            "${WRAPPER_DOXYGEN_DIR}/xml/"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/render.py"
            "${CMAKE_CURRENT_BINARY_DIR}/library.rst.jinja2"
            "${WRAPPER_DOXYGEN_DIR}/xml/index.xml"
    )

add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/api-c++/index.html"
    COMMAND "${PYTHON_EXECUTABLE}" -m sphinx . "${CMAKE_BINARY_DIR}/api-c++"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
            "${CMAKE_CURRENT_BINARY_DIR}/library.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/index.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/_static/logo.png"
            "${WRAPPER_DOXYGEN_DIR}/xml/index.xml"
)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/api-c++/)
add_custom_target(C++ApiDocs
        DEPENDS ${CMAKE_BINARY_DIR}/api-c++/index.html)
INSTALL(DIRECTORY "${CMAKE_BINARY_DIR}/api-c++/"
        DESTINATION "docs/api-c++/")
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../C++ API Reference.html"
        DESTINATION "docs/")
