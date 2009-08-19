INCLUDE(FindJava)

set(GCCXML_EXTRA_ARGS "" CACHE STRING "Additional arguments to gccxml, such as '--gccxml-compiler;msvc8'")
SET(GCCXML_ARGS) # start empty
FOREACH(subdir ${API_INCLUDE_DIRS})
    SET(GCCXML_ARGS ${GCCXML_ARGS} -I${subdir})
ENDFOREACH(subdir)
SET(GCCXML_ARGS ${GCCXML_ARGS} ${GCCXML_EXTRA_ARGS})

find_program(GCCXML_PATH gccxml PATH
    /usr/local/bin
    "C:/Program Files/gccxml 0.9/bin"
)

ADD_CUSTOM_COMMAND(OUTPUT OpenMMApi.xml COMMAND ${GCCXML_PATH} ${GCCXML_ARGS} ${CMAKE_SOURCE_DIR}/openmmapi/include/OpenMM.h -fxml=OpenMMApi.xml)
ADD_CUSTOM_COMMAND(OUTPUT OpenMMCWrapper.h COMMAND ${JAVA_RUNTIME} -jar ${CMAKE_SOURCE_DIR}/wrappers/saxonb9-1-0-7j/saxon9.jar -t -s:OpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/wrappers/CWrapper_Header.xslt -o:OpenMMCWrapper.h DEPENDS OpenMMApi.xml)
ADD_CUSTOM_COMMAND(OUTPUT OpenMMCWrapper.cpp COMMAND ${JAVA_RUNTIME} -jar ${CMAKE_SOURCE_DIR}/wrappers/saxonb9-1-0-7j/saxon9.jar -t -s:OpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/wrappers/CWrapper_Source.xslt -o:OpenMMCWrapper.cpp DEPENDS OpenMMApi.xml)
ADD_CUSTOM_COMMAND(OUTPUT OpenMMFortranModule.f90 COMMAND ${JAVA_RUNTIME} -jar ${CMAKE_SOURCE_DIR}/wrappers/saxonb9-1-0-7j/saxon9.jar -t -s:OpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/wrappers/FortranWrapper_Header.xslt -o:OpenMMFortranModule.f90 DEPENDS OpenMMApi.xml)
ADD_CUSTOM_COMMAND(OUTPUT OpenMMFortranWrapper.cpp COMMAND ${JAVA_RUNTIME} -jar ${CMAKE_SOURCE_DIR}/wrappers/saxonb9-1-0-7j/saxon9.jar -t -s:OpenMMApi.xml
     -xsl:${CMAKE_SOURCE_DIR}/wrappers/FortranWrapper_Source.xslt -o:OpenMMFortranWrapper.cpp DEPENDS OpenMMApi.xml)

ADD_CUSTOM_TARGET(ApiWrappers DEPENDS OpenMMCWrapper.h OpenMMCWrapper.cpp OpenMMFortranModule.f90 OpenMMFortranWrapper.cpp)


INSTALL_FILES(/include FILES OpenMMCWrapper.h OpenMMFortranModule.f90)
