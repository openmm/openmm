#---------------------------------------------------
# OpenMM Serialization Classes
#----------------------------------------------------

INSTALL_FILES(/include/openmm/serialization FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/openmm/serialization/SerializationNode.h)
INSTALL_FILES(/include/openmm/serialization FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/openmm/serialization/SerializationProxy.h)
INSTALL_FILES(/include/openmm/serialization FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/openmm/serialization/XmlSerializer.h)

SET(OPENMM_BUILD_SERIALIZATION_TESTS TRUE CACHE BOOL "Whether to build serialization test cases")
MARK_AS_ADVANCED(OPENMM_BUILD_SERIALIZATION_TESTS)
IF(BUILD_TESTING AND OPENMM_BUILD_SERIALIZATION_TESTS)
    ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTING AND OPENMM_BUILD_SERIALIZATION_TESTS)
