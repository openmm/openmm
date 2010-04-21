find_package(Doxygen)

set(OPENMM_GENERATE_API_DOCS OFF CACHE BOOL "Whether to create API documentation using Doxygen")
IF(DOXYGEN_FOUND)
    SET(DOXY_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")

    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in 
          ${DOXY_CONFIG}
          @ONLY )

    ADD_CUSTOM_COMMAND(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Generating API documentation using Doxygen"
        SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in") 
    ADD_CUSTOM_TARGET(DoxygenApiDocs 
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Generating API documentation using Doxygen"
        SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in") 
    FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/html/")

    IF(OPENMM_GENERATE_API_DOCS)
        INSTALL(DIRECTORY "${PROJECT_BINARY_DIR}/html/"
                DESTINATION "docs/api/")
    ENDIF(OPENMM_GENERATE_API_DOCS)
ELSE(DOXYGEN_FOUND)
ENDIF(DOXYGEN_FOUND)

