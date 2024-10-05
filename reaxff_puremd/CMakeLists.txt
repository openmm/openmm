# CMake project for reaxff-puremd
cmake_minimum_required(VERSION 3.8.1)

project(reaxff_puremd)


message(STATUS "Configuring and building the PuReMD software for ReaxFF via QMMM (may be time-consuming)...")

set(libspuremd_la_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}spuremd.la)
set(libspuremd_la_FULL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD/sPuReMD/lib/${libspuremd_la_NAME})
if("${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "Darwin")
       set(libspuremd_FULL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD/sPuReMD/lib/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}spuremd.1${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
       set(libspuremd_FULL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD/sPuReMD/lib/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}spuremd${CMAKE_SHARED_LIBRARY_SUFFIX}.1.0.0)
endif()

# ensure that install location for libspuremd exists before executing ExternalProject build
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)

# refresh timestemps of Autotools files (as storing in a git repo messes this up)
execute_process(
	COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD/aclocal.m4 ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD/configure ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD/Makefile.am ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD/Makefile.in
	RESULT_VARIABLE PUREMD_TOUCH_RESULT)

if(NOT PUREMD_TOUCH_RESULT EQUAL "0")
	message(FATAL_ERROR "reaxff-puremd build system timestamp update failed with error code ${PUREMD_TOUCH_RESULT}")
endif()

# execute ExternalProject build during the configure phase
# since target properties (dependent library, linker flags) are
# only known via the libtool file (*.la) _after_ the build completes

# ExternalProject configure phrase
execute_process(
	COMMAND ${CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD
	RESULT_VARIABLE PUREMD_CONFIGURE_RESULT)

if(NOT PUREMD_CONFIGURE_RESULT EQUAL "0")
	message(FATAL_ERROR "reaxff-puremd configure process failed with error code ${PUREMD_CONFIGURE_RESULT}")
endif()

# ExternalProject build and install phrases
execute_process(
	COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_SOURCE_DIR}/PuReMD
	RESULT_VARIABLE PUREMD_BUILD_RESULT)

if(NOT PUREMD_BUILD_RESULT EQUAL "0")
	message(FATAL_ERROR "reaxff-puremd build process failed with error code ${PUREMD_BUILD_RESULT}")
endif()

# after building ExternalProject, pull out library and linker dependencies
# for libspuremd from the associated libtool file (*.la)
file(READ ${libspuremd_la_FULL_PATH} libspuremd_la_TEXT)
string(REGEX MATCH "dependency_libs='([-A-Za-z0-9_ \t]*)'" _ "${libspuremd_la_TEXT}")
set(libspuremd_dep_libs ${CMAKE_MATCH_1})
string(REGEX MATCH "inherited_linker_flags='([-A-Za-z0-9_ \t]*)'" _ "${libspuremd_la_TEXT}")
set(libspuremd_link_flags ${CMAKE_MATCH_1})

add_library(libspuremd SHARED IMPORTED GLOBAL)
set_target_properties(libspuremd PROPERTIES IMPORTED_LOCATION ${libspuremd_FULL_PATH})
if(libspuremd_dep_libs)
	set_target_properties(libspuremd PROPERTIES LINK_INTERFACE_LIBRARIES "${libspuremd_dep_libs}")
endif()
if(libspuremd_link_flags)
	set_target_properties(libspuremd PROPERTIES LINK_FLAGS "${libspuremd_link_flags}")
endif()
