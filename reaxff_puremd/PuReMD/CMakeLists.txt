# CMake ExternalProject for reaxff-puremd
cmake_minimum_required(VERSION 3.8.1)

project(reaxff_puremd_ExternalProject)

include(ExternalProject)

if (OPENMP)
	set(reaxff_puremd_configure_options --libdir=${CMAKE_INSTALL_PREFIX}/lib --enable-openmp=yes --enable-qmmm=yes CC=${CMAKE_C_COMPILER})
else()
	set(reaxff_puremd_configure_options --libdir=${CMAKE_INSTALL_PREFIX}/lib --enable-openmp=no --enable-serial=yes --enable-qmmm=yes CC=${CMAKE_C_COMPILER})
endif()
if("${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "Darwin")
	set(libspuremd_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}spuremd.1${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(libspuremd_LINK1 ${CMAKE_SHARED_LIBRARY_PREFIX}spuremd.1.0.0${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(libspuremd_LINK2 ${CMAKE_SHARED_LIBRARY_PREFIX}spuremd${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(libspuremd_FULL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/sPuReMD/lib/.libs/${libspuremd_NAME})
else()
	set(libspuremd_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}spuremd${CMAKE_SHARED_LIBRARY_SUFFIX}.1.0.0)
	set(libspuremd_LINK1 ${CMAKE_SHARED_LIBRARY_PREFIX}spuremd${CMAKE_SHARED_LIBRARY_SUFFIX}.1)
	set(libspuremd_LINK2 ${CMAKE_SHARED_LIBRARY_PREFIX}spuremd${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(libspuremd_FULL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/sPuReMD/lib/.libs/${libspuremd_NAME})
endif()

ExternalProject_Add(reaxff_puremd_build
	PREFIX ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	
	CONFIGURE_COMMAND ./configure ${reaxff_puremd_configure_options}
	LOG_CONFIGURE TRUE

	BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} clean all "AR_FLAGS=cr" "CFLAGS=-O2 -fPIC"
	BUILD_IN_SOURCE TRUE
	BUILD_ALWAYS TRUE
	BUILD_BYPRODUCTS ${libspuremd_FULL_PATH}
	LOG_BUILD TRUE

	INSTALL_COMMAND ${CMAKE_COMMAND} -E copy ${libspuremd_FULL_PATH} ${CMAKE_INSTALL_PREFIX}/lib/${libspuremd_NAME}
	COMMAND ${CMAKE_COMMAND} -E create_symlink ${libspuremd_NAME} ${CMAKE_INSTALL_PREFIX}/lib/${libspuremd_LINK1}
	COMMAND ${CMAKE_COMMAND} -E create_symlink ${libspuremd_LINK1} ${CMAKE_INSTALL_PREFIX}/lib/${libspuremd_LINK2}
	LOG_INSTALL TRUE

	TEST_COMMAND ""
)
