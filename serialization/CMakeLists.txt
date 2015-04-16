#---------------------------------------------------
# OpenMM Serialization Classes
#----------------------------------------------------

INSTALL_FILES(/include/openmm/serialization FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/openmm/serialization/SerializationNode.h)
INSTALL_FILES(/include/openmm/serialization FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/openmm/serialization/SerializationProxy.h)
INSTALL_FILES(/include/openmm/serialization FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/openmm/serialization/XmlSerializer.h)

IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTING)
