set(STAGING_OUTPUT_FILES "") # Will contain all required package files
file(GLOB STAGING_INPUT_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "app.rst.jinja2"
    "library.rst.jinja2"
    "index.rst"
    "render.py"
    "conf.py"
    "process-docstring.py"
    "_static/logo.png"
    "_templates/class.rst"
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
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/app.rst"
           "${CMAKE_CURRENT_BINARY_DIR}/library.rst"
    COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/render.py"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/app.rst.jinja2"
            "${CMAKE_CURRENT_BINARY_DIR}/library.rst.jinja2"
            "${CMAKE_CURRENT_BINARY_DIR}/render.py"
)

add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/api-python/index.html"
    COMMAND "${PYTHON_EXECUTABLE}" -m sphinx . "${CMAKE_BINARY_DIR}/api-python"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
            "${CMAKE_CURRENT_BINARY_DIR}/process-docstring.py"
            "${CMAKE_CURRENT_BINARY_DIR}/app.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/library.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/index.rst"
            "${CMAKE_CURRENT_BINARY_DIR}/_static/logo.png"
            "${CMAKE_CURRENT_BINARY_DIR}/_templates/class.rst"
            # PythonInstall
)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/api-python/)
add_custom_target(PythonApiDocs DEPENDS ${CMAKE_BINARY_DIR}/api-python/index.html)
INSTALL(DIRECTORY "${CMAKE_BINARY_DIR}/api-python/"
        DESTINATION "docs/api-python/")
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../Python API Reference.html"
        DESTINATION "docs/")
