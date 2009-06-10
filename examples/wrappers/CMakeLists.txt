# Build a C wrapper for OpenMM


# Should rebuild this if anything changes, but here we're
# just listing the local cpp and h files as dependencies.
ADD_LIBRARY(${OpenMM_CWRAPPER_LIB}  STATIC 
	OpenMM_CWrapper.cpp OpenMM_CWrapper.h) 	

SET_TARGET_PROPERTIES(${OpenMM_CWRAPPER_LIB} PROPERTIES
	PROJECT_LABEL "Example - C Wrapper")

# Don't install yet
#INSTALL_TARGETS(/lib${LIB64}/plugins 
#                RUNTIME_DIRECTORY /lib${LIB64}/plugins
#                ${OpenMM_CWRAPPER_LIB})







