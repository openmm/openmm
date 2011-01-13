#set(GCCXML_EXTRA_ARGS "" CACHE STRING "Additional arguments to gccxml, such as '--gccxml-compiler;msvc8'")
#SET(GCCXML_ARGS) # start empty
#FOREACH(subdir ${API_INCLUDE_DIRS})
#    SET(GCCXML_ARGS ${GCCXML_ARGS} -I${subdir})
#ENDFOREACH(subdir)
#SET(GCCXML_ARGS ${GCCXML_ARGS} ${GCCXML_EXTRA_ARGS})

SET(SAXON_JAR ${CMAKE_SOURCE_DIR}/wrappers/saxonb9-1-0-7j/saxon9.jar)

# OpenMMAmoeba

ADD_CUSTOM_COMMAND(OUTPUT AmoebaOpenMMApi.xml COMMAND ${GCCXML_PATH} ${GCCXML_ARGS} -I${CMAKE_SOURCE_DIR}/openmmapi/include -I${CMAKE_SOURCE_DIR}/olla/include ${CMAKE_SOURCE_DIR}/plugins/amoeba/openmmapi/include/OpenMMAmoeba.h -fxml=AmoebaOpenMMApi.xml)

ADD_CUSTOM_COMMAND(OUTPUT AmoebaOpenMMCWrapper.h COMMAND ${JAVA_RUNTIME} -jar ${SAXON_JAR} -t -s:AmoebaOpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/plugins/amoeba/wrappers/CWrapper_Header.xslt -o:AmoebaOpenMMCWrapper.h DEPENDS AmoebaOpenMMApi.xml)
#
ADD_CUSTOM_COMMAND(OUTPUT AmoebaOpenMMCWrapper.cpp COMMAND ${JAVA_RUNTIME} -jar ${SAXON_JAR} -t -s:AmoebaOpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/plugins/amoeba/wrappers/CWrapper_Source.xslt -o:AmoebaOpenMMCWrapper.cpp DEPENDS AmoebaOpenMMApi.xml)

ADD_CUSTOM_COMMAND(OUTPUT AmoebaOpenMMFortranModule.f90 COMMAND ${JAVA_RUNTIME} -jar ${SAXON_JAR} -t -s:AmoebaOpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/plugins/amoeba/wrappers/FortranWrapper_Header.xslt -o:AmoebaOpenMMFortranModule.f90 DEPENDS AmoebaOpenMMApi.xml ${CMAKE_SOURCE_DIR}/plugins/amoeba/wrappers/FortranWrapper_Header.xslt)

ADD_CUSTOM_COMMAND(OUTPUT AmoebaOpenMMFortranWrapper.cpp COMMAND ${JAVA_RUNTIME} -jar ${SAXON_JAR} -t -s:AmoebaOpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/plugins/amoeba/wrappers/FortranWrapper_Source.xslt -o:AmoebaOpenMMFortranWrapper.cpp DEPENDS AmoebaOpenMMApi.xml)

ADD_CUSTOM_TARGET(AmoebaApiWrappers DEPENDS AmoebaOpenMMCWrapper.h AmoebaOpenMMCWrapper.cpp  AmoebaOpenMMFortranModule.f90 AmoebaOpenMMFortranWrapper.cpp)

INSTALL_FILES(/include FILES AmoebaOpenMMCWrapper.h AmoebaOpenMMFortranModule.f90 )
