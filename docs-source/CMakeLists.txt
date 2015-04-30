#
# Build and install API documentation
#

find_package(Doxygen QUIET)
mark_as_advanced(CLEAR DOXYGEN_EXECUTABLE)
IF(DOXYGEN_EXECUTABLE)
    # Generate C++ API documentation

    SET(DOXY_CONFIG_C++ "${CMAKE_BINARY_DIR}/DoxyfileC++")
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/DoxyfileC++.in 
          ${DOXY_CONFIG_C++}
          @ONLY )
    FILE(GLOB_RECURSE OPENMM_INCLUDES "openmm/include/*.h")
    FILE(GLOB_RECURSE OLLA_INCLUDES "olla/include/*.h")
    ADD_CUSTOM_COMMAND(
        OUTPUT "${CMAKE_BINARY_DIR}/api-c++/index.html"
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG_C++}
        DEPENDS ${OPENMM_INCLUDES} ${OLLA_INCLUDES}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        COMMENT "Generating C++ API documentation using Doxygen")
    ADD_CUSTOM_TARGET(C++ApiDocs 
        DEPENDS "${CMAKE_BINARY_DIR}/api-c++/index.html"
        COMMENT "Generating C++ API documentation using Doxygen"
        SOURCES 
            "${CMAKE_CURRENT_SOURCE_DIR}/DoxyfileC++.in"
            ${OPENMM_INCLUDES}
            ${OLLA_INCLUDES}
    ) 
    FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/api-c++/")
    INSTALL(DIRECTORY "${PROJECT_BINARY_DIR}/api-c++/"
            DESTINATION "docs/api-c++/")
    INSTALL(FILES "C++ API Reference.html"
            DESTINATION "docs/")
    ADD_CUSTOM_TARGET(DoxygenApiDocs
        DEPENDS "${CMAKE_BINARY_DIR}/api-c++/index.html"
        COMMENT "Generating C++ API documentation using Doxygen"
        SOURCES 
        "${CMAKE_CURRENT_SOURCE_DIR}/DoxyfileC++.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/DoxyfilePython.in"
        ${OPENMM_INCLUDES}
        ${OLLA_INCLUDES}
    )
    set(OPENMM_GENERATE_API_DOCS OFF CACHE BOOL "Whether to create API documentation using Doxygen")
    IF (OPENMM_GENERATE_API_DOCS)
        SET_TARGET_PROPERTIES(DoxygenApiDocs PROPERTIES EXCLUDE_FROM_ALL FALSE)
    ENDIF (OPENMM_GENERATE_API_DOCS)

    # Generate Python API documentation

    IF (OPENMM_BUILD_PYTHON_WRAPPERS)
        SET(DOXY_CONFIG_PYTHON "${CMAKE_BINARY_DIR}/DoxyfilePython")
        CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/DoxyfilePython.in 
            ${DOXY_CONFIG_PYTHON}
            @ONLY )
        ADD_CUSTOM_COMMAND(
            OUTPUT "${CMAKE_BINARY_DIR}/api-python/index.html"
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG_PYTHON}
            DEPENDS RunSwig
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
            COMMENT "Generating Python API documentation using Doxygen")
        ADD_CUSTOM_TARGET(PythonApiDocs 
            DEPENDS "${CMAKE_BINARY_DIR}/api-python/index.html"
            COMMENT "Generating Python API documentation using Doxygen"
            SOURCES 
                "${CMAKE_CURRENT_SOURCE_DIR}/DoxyfilePython.in"
                ${OPENMM_INCLUDES}
                ${OLLA_INCLUDES}
        ) 
        FILE(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/api-python/")
        INSTALL(DIRECTORY "${PROJECT_BINARY_DIR}/api-python/"
                DESTINATION "docs/api-python/")
        INSTALL(FILES "Python API Reference.html"
                DESTINATION "docs/")
        ADD_DEPENDENCIES(DoxygenApiDocs PythonApiDocs)
    ENDIF (OPENMM_BUILD_PYTHON_WRAPPERS)
ENDIF(DOXYGEN_EXECUTABLE)


#
# Build and install the User Guide and Developer Guide
#

SET(SPHINX_BUILD_DIR "${CMAKE_BINARY_DIR}/sphinx-docs/")
FILE(MAKE_DIRECTORY "${SPHINX_BUILD_DIR}")

FILE(GLOB_RECURSE USER_GUIDE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/usersguide/*)
FILE(GLOB_RECURSE DEVELOPER_GUIDE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/developerguide/*)
ADD_CUSTOM_COMMAND(
    OUTPUT "${SPHINX_BUILD_DIR}/userguide/latex/OpenMMUsersGuide.pdf"
    COMMAND "${CMAKE_MAKE_PROGRAM}" BUILDDIR="${SPHINX_BUILD_DIR}/userguide" OPENMM_VERSION="${OPENMM_MAJOR_VERSION}.${OPENMM_MINOR_VERSION}" latexpdf
    DEPENDS ${USER_GUIDE_FILES}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/usersguide"
    COMMENT "Generating PDF user guide"
)
ADD_CUSTOM_COMMAND(
    OUTPUT "${SPHINX_BUILD_DIR}/developerguide/latex/OpenMMDeveloperGuide.pdf"
    COMMAND "${CMAKE_MAKE_PROGRAM}" BUILDDIR="${SPHINX_BUILD_DIR}/developerguide" OPENMM_VERSION="${OPENMM_MAJOR_VERSION}.${OPENMM_MINOR_VERSION}" latexpdf
    DEPENDS ${DEVELOPER_GUIDE_FILES}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/developerguide"
    COMMENT "Generating PDF developer guide"
)
ADD_CUSTOM_TARGET(sphinxpdf
    DEPENDS "${SPHINX_BUILD_DIR}/userguide/latex/OpenMMUsersGuide.pdf" "${SPHINX_BUILD_DIR}/developerguide/latex/OpenMMDeveloperGuide.pdf"
)

ADD_CUSTOM_COMMAND(
    OUTPUT "${SPHINX_BUILD_DIR}/userguide/html/index.html"
    COMMAND "${CMAKE_MAKE_PROGRAM}" BUILDDIR="${SPHINX_BUILD_DIR}/userguide" OPENMM_VERSION="${OPENMM_MAJOR_VERSION}.${OPENMM_MINOR_VERSION}" html
    DEPENDS ${USER_GUIDE_FILES}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/usersguide"
    COMMENT "Generating HTML user guide"
)
ADD_CUSTOM_COMMAND(
    OUTPUT "${SPHINX_BUILD_DIR}/developerguide/html/index.html"
    COMMAND "${CMAKE_MAKE_PROGRAM}" BUILDDIR="${SPHINX_BUILD_DIR}/developerguide" OPENMM_VERSION="${OPENMM_MAJOR_VERSION}.${OPENMM_MINOR_VERSION}" html
    DEPENDS ${DEVELOPER_GUIDE_FILES}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/developerguide"
    COMMENT "Generating HTML developer guide"
)
ADD_CUSTOM_TARGET(sphinxhtml
    DEPENDS "${SPHINX_BUILD_DIR}/userguide/html/index.html" "${SPHINX_BUILD_DIR}/developerguide/html/index.html"
)

install(FILES "${SPHINX_BUILD_DIR}/userguide/latex/OpenMMUsersGuide.pdf" "${SPHINX_BUILD_DIR}developerguide/latex/OpenMMDeveloperGuide.pdf"
    DESTINATION docs/ OPTIONAL)

FILE(GLOB LICENSE_FILES "licenses/*.txt")
install(FILES ${LICENSE_FILES}
    DESTINATION licenses/)
