#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([PuReMD], [1.0], [ohearnku@msu.edu hma@msu.edu])
# Do not allow AC_PROG_CC to set CFLAGS (this line must be after AC_INIT but before AC_PROG_CC)
save_CFLAGS="${CFLAGS}"
: ${CFLAGS=""}
AM_INIT_AUTOMAKE([1.15 subdir-objects -Wall -Werror])
# Enable silent build rules by default.
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])], [AC_SUBST([AM_DEFAULT_VERBOSITY],[1])])
AM_PROG_AR
LT_PREREQ([2.2])
LT_INIT([dlopen])

AC_CONFIG_MACRO_DIR([m4])

# Headline formatter
AC_DEFUN([CONFIGURE_HEADLINE],
[
        echo; echo "+++ $1 +++"
])

AC_ARG_ENABLE([serial],
	      [AS_HELP_STRING([--enable-serial],
			      [enable serial support @<:@default: no@:>@])],
	      [pack_serial_enabled=${enableval}], [pack_serial_enabled=no])
AC_ARG_ENABLE([openmp],
	      [AS_HELP_STRING([--enable-openmp],
			      [enable OpenMP support @<:@default: yes@:>@])],
	      [pack_openmp_enabled=${enableval}], [pack_openmp_enabled=yes])

# Build standalone PuReMD code.
AC_ARG_ENABLE([puremd-standalone],
	      [AS_HELP_STRING([--enable-puremd-standalone],
			      [enable build for standalone PuReMD code @<:@default: yes@:>@])],
	      [puremd_standalone=${enableval}], [puremd_standalone=yes])

# Build code in QM/MM mode.
AC_ARG_ENABLE([qmmm],
	      [AS_HELP_STRING([--enable-qmmm],
			      [enable build for code in QM/MM mode @<:@default: no@:>@])],
	      [qmmm=${enableval}], [qmmm=no])

# Build LAMMPS/reaxc integration code.
AC_ARG_ENABLE([lammps-reaxc],
	      [AS_HELP_STRING([--enable-lammps-reaxc],
			      [enable build for LAMMPS/reaxc integration code @<:@default: no@:>@])],
	      [lammps_reaxc=${enableval}], [lammps_reaxc=no])

# Build dual QEq solver code.
AC_ARG_ENABLE([qeq-dual-solver],
	      [AS_HELP_STRING([--enable-qeq-dual-solver],
			      [enable dual QEq charge solver code @<:@default: yes@:>@])],
	      [qeq_dual_solver=${enableval}], [qeq_dual_solver=yes])

# Build neutral territory code.
AC_ARG_ENABLE([neutral-territory],
	      [AS_HELP_STRING([--enable-neutral-territory],
			      [enable neutral territory code @<:@default: no@:>@])],
	      [neutral_territory=${enableval}], [neutral_territory=no])

# Build neutral territory code.
AC_ARG_ENABLE([log-perf],
	      [AS_HELP_STRING([--enable-log-perf],
			      [enable performance logging code @<:@default: yes@:>@])],
	      [log_perf=${enableval}], [log_perf=yes])

# Provides debug compilation mode.
AC_ARG_ENABLE([debug],
	      [AS_HELP_STRING([--enable-debug],
			      [enable debug support @<:@default: no@:>@])],
	      [DEBUG=${enableval}], [DEBUG=no]
	      )
if test "x${DEBUG}" = "xyes"; then
#	#TODO: fix exporting to subdirs
#	# See: http://stackoverflow.com/questions/34124337/changing-flags-in-configure-ac-vs-caching-with-subprojects
	export BUILD_DEBUG="yes"
fi

# gprof flags.
AC_ARG_ENABLE([gprof],
	[AS_HELP_STRING([--enable-gprof],
		[enable support for profiling with gprof @<:@default: no@:>@])],
	[case "${enableval}" in
		yes) 
			export BUILD_GPROF="yes"
			;;
		no)
			;;
		*)
			AC_MSG_ERROR([bad value ${enableval} for --enable-gprof (only yes or no are possible)]) ;;
	esac],
	[]	      
)


# Checks for compilers.
AC_LANG([C])
AC_PROG_CC([icc gcc clang cc])
AC_PROG_CPP
CFLAGS="${save_CFLAGS}"

AC_LANG([C++])
AC_PROG_CXX([icpc g++ clang++ CC c++])
AC_PROG_CXXCPP

AC_CONFIG_HEADERS([common/include/config.h])

# sPuReMD
if test "x${pack_serial_enabled}" = "xyes" || test "x${pack_openmp_enabled}" = "xyes"; then
	if test "x${pack_serial_enabled}" = "xyes" || test "x${pack_openmp_enabled}" != "xyes"; then
		BUILD_OPENMP="no"
	else
		BUILD_OPENMP="yes"
	fi

	AC_LANG_PUSH([C])

	AC_CONFIG_SRCDIR([sPuReMD/src/traj.h])

	# Checks for libraries.
	AC_SEARCH_LIBS([exp], [m])
	AC_SEARCH_LIBS([sqrt], [m])
	AC_SEARCH_LIBS([pow], [m])
	AC_SEARCH_LIBS([acos], [m])
	AC_SEARCH_LIBS([cos], [m])
	AC_SEARCH_LIBS([sin], [m])
	AC_SEARCH_LIBS([tan], [m])
	AC_SEARCH_LIBS([fabs], [m])

	# Checks for typedefs, structures, and compiler characteristics.
	AC_C_INLINE

	# Checks for library functions.
	AC_FUNC_MALLOC
	AC_FUNC_REALLOC
	AC_FUNC_STRTOD
	AC_CHECK_FUNCS([gettimeofday memset])

	# Check for compiler vendor. If the compiler is recognized,
	#   the variable ax_cv_c_compiler_vendor is set accordingly.
	AX_COMPILER_VENDOR
	if test "x${ax_cv_c_compiler_vendor}" = "xgnu"; then
		if test "x${BUILD_DEBUG}" = "xyes"; then
			DEFAULT_CFLAGS="-Wall -O0 -g2 -D_GLIBCXX_DEBUG"
			AC_DEFINE([DEBUG], [1], [Define to 1 if you have debug support enabled.])
		else
			DEFAULT_CFLAGS="-Wall -O2"
		fi
		if test "x${BUILD_GPROF}" = "xyes"; then
			DEFAULT_CFLAGS+=" -pg"
		fi
	elif test "x${ax_cv_c_compiler_vendor}" = "xintel"; then
		if test "x${BUILD_DEBUG}" = "xyes"; then
			DEFAULT_CFLAGS="-Wall -O0"
			AC_DEFINE([DEBUG], [1], [Define to 1 if you have debug support enabled.])
		else
			DEFAULT_CFLAGS="-Wall -O2"
		fi
		if test "x${BUILD_GPROF}" = "xyes"; then
			DEFAULT_CFLAGS+=" -p"
		fi
	fi
	AC_SUBST([DFT_CFLAGS], ["${DEFAULT_CFLAGS}"])

	# Check for OpenMP support.
	if test "x${BUILD_OPENMP}" = "xyes"; then
		AX_OPENMP

		if test "x${OPENMP_CFLAGS}" = "x"; then
			AC_MSG_WARN([
-----------------------------------------------
 Unable to find OpenMP support on this system.
 Building a single-threaded version.
-----------------------------------------------])
		else
			# bug due to recent Intel compiler change (?)
			if test "x${ax_cv_c_compiler_vendor}" = "xintel"; then
				OPENMP_CFLAGS="-qopenmp"
			fi
			OPENMP_LIBS="-lgomp"
		fi

	#	AC_MSG_NOTICE([OPENMP_CFLAGS=${OPENMP_CFLAGS}])
	#	AC_MSG_NOTICE([OPENMP_LIBS=${OPENMP_LIBS}])
		AC_SUBST([OMP_CFLAGS], ["${OPENMP_CFLAGS}"])
		AC_SUBST([OMP_LIBS], ["${OPENMP_LIBS}"])
	fi

	# Check for zlib
	ZLIB_LIBS=
	AC_CHECK_HEADERS([zlib.h], [ZLIB_FOUND_HEADERS="yes"])
	if test "x${ZLIB_FOUND_HEADERS}" = "xyes"; then
		# BUG: AC_SEARCH_LIBS fails to reset LIBS when nested in conditional statements, so do it manually
		save_LIBS="${LIBS}"
		ZLIB_FOUND_LIBS="yes"
		AC_SEARCH_LIBS([gzopen], [z], [], [ZLIB_FOUND_LIBS="no"], [])
		AC_SEARCH_LIBS([gzeof], [z], [], [ZLIB_FOUND_LIBS="no"], [])
		AC_SEARCH_LIBS([gzgets], [z], [], [ZLIB_FOUND_LIBS="no"], [])
		AC_SEARCH_LIBS([gzseek], [z], [], [ZLIB_FOUND_LIBS="no"], [])
		AC_SEARCH_LIBS([gzclose], [z], [], [ZLIB_FOUND_LIBS="no"], [])
		AC_CHECK_TYPES([gzFile], [], [ZLIB_FOUND_LIBS="no"], [#include<zlib.h>])
		LIBS="${save_LIBS}"
		AS_IF([test "x${ZLIB_FOUND_LIBS}" != "xyes"],
		      [AC_MSG_ERROR([Unable to find zlib compatible library.])])
		ZLIB_LIBS="-lz"
		AC_DEFINE([HAVE_ZLIB], [1], [Define to 1 if you have zlib support enabled.])
	else
		AC_MSG_WARN([
-----------------------------------------------
 Unable to find compatible zlib on this system.
 Disabling support for dependent functionality.
-----------------------------------------------])
	fi
	AC_SUBST([Z_LIBS], ["${ZLIB_LIBS}"])

	# Check for LAPACKE
	LAPACKE_LIBS=
	AC_CHECK_HEADERS([mkl.h], [MKL_FOUND_HEADERS="yes"])
	if test "x${MKL_FOUND_HEADERS}" = "xyes"; then
		# BUG: AC_SEARCH_LIBS fails to reset LIBS when nested in conditional statements, so do it manually
		save_LIBS="${LIBS}"
		AC_SEARCH_LIBS([LAPACKE_dgels], [mkl_intel_ilp64],
			       [MKL_FOUND_LIBS="yes"], [MKL_FOUND_LIBS="no"],
			       [-lmkl_sequential -lmkl_core -lpthread -lm -ldl])
		LIBS="${save_LIBS}"
		AS_IF([test "x${MKL_FOUND_LIBS}" != "xyes"],
		      [AC_MSG_ERROR([Unable to find MKL LAPACKE library.])])
		LAPACKE_LIBS="-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
		AC_DEFINE([HAVE_LAPACKE_MKL], [1], [Define to 1 if you have MKL LAPACKE support enabled.])
	else
		AC_CHECK_HEADERS([lapacke.h], [LAPACKE_FOUND_HEADERS="yes"])
		if test "x${LAPACKE_FOUND_HEADERS}" = "xyes"; then
			# BUG: AC_SEARCH_LIBS fails to reset LIBS when nested in conditional statements, so do it manually
			save_LIBS="${LIBS}"
			AC_SEARCH_LIBS([LAPACKE_dgels], [lapacke],
				       [LAPACKE_FOUND_LIBS="yes"], [LAPACKE_FOUND_LIBS="no"],
				       [-llapack])
			LIBS="${save_LIBS}"
			AS_IF([test "x${LAPACKE_FOUND_LIBS}" != "xyes"],
			      [AC_MSG_ERROR([Unable to find LAPACKE library.])])
			LAPACKE_LIBS="-llapacke -llapack"
			AC_DEFINE([HAVE_LAPACKE], [1], [Define to 1 if you have LAPACKE support enabled.])
		else
			AC_MSG_WARN([
-----------------------------------------------
 Unable to find LAPACKE on this system.
 Disabling support for dependent methods.
-----------------------------------------------])
		fi
	fi
	AC_SUBST([L_LIBS], ["${LAPACKE_LIBS}"])

	# Check for Tensorflow
	TENSORFLOW_LIBS=
	AC_CHECK_HEADERS([tensorflow/c/c_api.h], [TENSORFLOW_FOUND_HEADERS="yes"])
	if test "x${TENSORFLOW_FOUND_HEADERS}" = "xyes"; then
		# BUG: AC_SEARCH_LIBS fails to reset LIBS when nested in conditional statements, so do it manually
		save_LIBS="${LIBS}"
		AC_SEARCH_LIBS([TF_NewStatus], [tensorflow],
			       [TENSORFLOW_FOUND_LIBS="yes"], [TENSORFLOW_FOUND_LIBS="no"],
			       [-ltensorflow_framework])
		LIBS="${save_LIBS}"
		AS_IF([test "x${TENSORFLOW_FOUND_LIBS}" != "xyes"],
		      [AC_MSG_WARN([
----------------------------------------------------
 Unable to find Tensorflow libraries on this system.
 Disabling support for dependent methods.
----------------------------------------------------])])
		if test "x${TENSORFLOW_FOUND_HEADERS}" = "xyes" && test "x${TENSORFLOW_FOUND_LIBS}" = "xyes"; then
			TENSORFLOW_LIBS="-ltensorflow"
		fi
		AC_DEFINE([HAVE_TENSORFLOW], [1], [Define to 1 if you have TENSORFLOW support enabled.])
	else
		AC_MSG_WARN([
--------------------------------------------------
 Unable to find Tensorflow headers on this system.
 Disabling support for dependent methods.
--------------------------------------------------])
	fi
	AC_SUBST([T_LIBS], ["${TENSORFLOW_LIBS}"])

	# Check for GNU Scientific Library (for energy minimization)
	GSL_LIBS=
	AC_CHECK_HEADERS([gsl/gsl_sf_bessel.h], [GSL_FOUND_HEADERS="yes"])
	if test "x${GSL_FOUND_HEADERS}" = "xyes"; then
		# BUG: AC_SEARCH_LIBS fails to reset LIBS when nested in conditional statements, so do it manually
		save_LIBS="${LIBS}"
		AC_SEARCH_LIBS([gsl_sf_bessel_J0], [gsl],
			       [GSL_FOUND_LIBS="yes"], [GSL_FOUND_LIBS="no"],
			       [-lgslcblas -m])
		LIBS="${save_LIBS}"
		AS_IF([test "x${GSL_FOUND_LIBS}" != "xyes"],
		      [AC_MSG_WARN([
----------------------------------------------------
 Unable to find GNU scientific library (GSL) on this system.
 Disabling support for dependent methods.
----------------------------------------------------])])
		if test "x${GSL_FOUND_HEADERS}" = "xyes" && test "x${GSL_FOUND_LIBS}" = "xyes"; then
			GSL_LIBS="-lgsl -lgslcblas -m"
		fi
		AC_DEFINE([HAVE_GSL], [1], [Define to 1 if you have GNU scientific library (GSL) support enabled.])
	else
		AC_MSG_WARN([
--------------------------------------------------
 Unable to find GNU scientific library (GSL) headers on this system.
 Disabling support for dependent methods.
--------------------------------------------------])
	fi
	AC_SUBST([G_LIBS], ["${GSL_LIBS}"])

	# Build code in QM/MM mode
	AS_IF([test "x${qmmm}" = "xyes"],
	      [AC_DEFINE([QMMM], [1], [Define to 1 to build PuReMD code in QMMM mode.])])

	AC_LANG_POP([C])
fi
AM_CONDITIONAL([BUILD_S_OMP], [test "x${pack_serial_enabled}" = "xyes" || test "x${pack_openmp_enabled}" = "xyes"])


# Tests using Google C++ testing framework (gtest)
AC_LANG_PUSH([C++])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
  [[#ifndef __cplusplus
    #error "broken C++"
    #endif]])],,
  [CXX=;])
if test "x${CXX}" != "x"; then
	AC_CHECK_HEADER([gtest/gtest.h], [HAVE_GTEST_H="yes"])
	if test "x${HAVE_GTEST_H}" != "xyes"; then
		AC_MSG_WARN([gtest.h not found. Unable to build tests. Continuing...])
	else
		ACX_PTHREAD([have_pthread=yes])

		# define gtest variables 
		if test "x${have_pthread}" = "xyes"; then
			GTEST_CPPFLAGS="-DGTEST_HAS_PTHREAD=1"
			GTEST_CXXFLAGS="${PTHREAD_CFLAGS}"
			GTEST_LDFLAGS=
			GTEST_LIBS="${PTHREAD_LIBS}"
		else
			GTEST_CPPFLAGS="-DGTEST_HAS_PTHREAD=0"
			GTEST_CXXFLAGS=
			GTEST_LDFLAGS=
			GTEST_LIBS=
		fi
	fi
else
	AC_MSG_WARN([C++ compiler not found. Unable to build tests. Continuing...])
fi
AM_CONDITIONAL([BUILD_TEST], [test "x${HAVE_GTEST_H}" = "xyes"])
AC_LANG_POP([C++])


AC_CONFIG_FILES([Makefile])

AC_OUTPUT
