# Build a C wrapper for OpenMM


# Should rebuild this if anything changes, but here we're
# just listing the local cpp and h files as dependencies.
ADD_LIBRARY(${OpenMM_CWRAPPER_LIB}  SHARED
	OpenMM_CWrapper.cpp OpenMM_CWrapper.h) 	
TARGET_LINK_LIBRARIES(${OpenMM_CWRAPPER_LIB} ${SHARED_TARGET})

SET_TARGET_PROPERTIES(${OpenMM_CWRAPPER_LIB} PROPERTIES
	PROJECT_LABEL "Example - C Wrapper")

INSTALL(FILES OpenMM_CWrapper.cpp 
              OpenMM_CWrapper.h 
              OpenMM_Module.f90 
        DESTINATION examples)
