set(STAGING_OUTPUT_FILES "") # Will contain all required package files
file(GLOB STAGING_INPUT_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "index.rst"
    "core.rst.jinja2"
    "integrators.rst.jinja2"
    "forces.rst.jinja2"
    "extras.rst.jinja2"
    "render.py"
    "_static/logo.png"
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


set(DOXYGEN_XML_DIRECTORY "${CMAKE_BINARY_DIR}/wrappers/python/src/swig_doxygen/doxygen/xml/")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/conf.py" "${CMAKE_CURRENT_BINARY_DIR}/conf.py")


add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/core.rst"
           "${CMAKE_CURRENT_BINARY_DIR}/integrators.rst"
           "${CMAKE_CURRENT_BINARY_DIR}/forces.rst"
           "${CMAKE_CURRENT_BINARY_DIR}/extras.rst"
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/render.py"
            "${CMAKE_BINARY_DIR}/wrappers/python/src/swig_doxygen/doxygen/xml/"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/render.py"
            "${CMAKE_CURRENT_BINARY_DIR}/core.rst.jinja2"
            "${CMAKE_CURRENT_BINARY_DIR}/integrators.rst.jinja2"
            "${CMAKE_CURRENT_BINARY_DIR}/forces.rst.jinja2"
            "${CMAKE_CURRENT_BINARY_DIR}/extras.rst.jinja2"
            "${CMAKE_BINARY_DIR}/wrappers/python/src/swig_doxygen/doxygen/xml/index.xml"
            RunSwig
)

add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/api-c++/index.html"
    COMMAND "${PYTHON_EXECUTABLE}" -m sphinx . "${CMAKE_BINARY_DIR}/api-c++"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
            "${CMAKE_CURRENT_BINARY_DIR}/core.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/forces.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/extras.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/integrators.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/index.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/_static/logo.png"
)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/api-c++/)
add_custom_target(C++ApiDocs
        DEPENDS ${CMAKE_BINARY_DIR}/api-c++/index.html)
INSTALL(DIRECTORY "${CMAKE_BINARY_DIR}/api-c++/"
        DESTINATION "docs/api-c++/")
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../C++ API Reference.html"
        DESTINATION "docs/")
