ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS = -O2 -funroll-loops -Wall

SUBDIRS =
DIST_SUBDIRS =
EXTRA_DIST = 

# sPuReMD
if BUILD_S_OMP
lib_LTLIBRARIES = sPuReMD/lib/libspuremd.la
sPuReMD_lib_libspuremd_la_SOURCES = sPuReMD/src/allocate.c \
				    sPuReMD/src/allocate.h \
				    sPuReMD/src/analyze.c \
				    sPuReMD/src/analyze.h \
				    sPuReMD/src/bond_orders.c \
				    sPuReMD/src/bond_orders.h \
				    sPuReMD/src/bonds.c \
				    sPuReMD/src/bonds.h \
				    sPuReMD/src/box.c \
				    sPuReMD/src/box.h \
				    sPuReMD/src/charges.c \
				    sPuReMD/src/charges.h \
				    sPuReMD/src/control.c \
				    sPuReMD/src/control.h \
				    sPuReMD/src/ffield.c \
				    sPuReMD/src/ffield.h \
				    sPuReMD/src/forces.c \
				    sPuReMD/src/forces.h \
				    sPuReMD/src/hydrogen_bonds.c \
				    sPuReMD/src/hydrogen_bonds.h \
				    sPuReMD/src/geo_tools.c \
				    sPuReMD/src/geo_tools.h \
				    sPuReMD/src/grid.c \
				    sPuReMD/src/grid.h \
				    sPuReMD/src/init_md.c \
				    sPuReMD/src/init_md.h \
				    sPuReMD/src/io_tools.c \
				    sPuReMD/src/io_tools.h \
				    sPuReMD/src/integrate.c \
				    sPuReMD/src/integrate.h \
				    sPuReMD/src/lin_alg.c \
				    sPuReMD/src/lin_alg.h \
				    sPuReMD/src/list.c \
				    sPuReMD/src/list.h \
				    sPuReMD/src/lookup.c \
				    sPuReMD/src/lookup.h \
				    sPuReMD/src/multi_body.c \
				    sPuReMD/src/multi_body.h \
				    sPuReMD/src/neighbors.c \
				    sPuReMD/src/neighbors.h \
				    sPuReMD/src/nonbonded.c \
				    sPuReMD/src/nonbonded.h \
				    sPuReMD/src/random.c \
				    sPuReMD/src/random.h \
				    sPuReMD/src/reax_types.h \
				    sPuReMD/src/reset_tools.c \
				    sPuReMD/src/reset_tools.h \
				    sPuReMD/src/restart.c \
				    sPuReMD/src/restart.h \
				    sPuReMD/src/spuremd.c \
				    sPuReMD/src/spuremd.h \
				    sPuReMD/src/system_props.c \
				    sPuReMD/src/system_props.h \
				    sPuReMD/src/torsion_angles.c \
				    sPuReMD/src/torsion_angles.h \
				    sPuReMD/src/tool_box.c \
				    sPuReMD/src/tool_box.h \
				    sPuReMD/src/traj.c \
				    sPuReMD/src/traj.h \
				    sPuReMD/src/valence_angles.c \
				    sPuReMD/src/valence_angles.h \
				    sPuReMD/src/vector.c \
				    sPuReMD/src/vector.h
sPuReMD_lib_libspuremd_la_CFLAGS = @DFT_CFLAGS@ @OMP_CFLAGS@
sPuReMD_lib_libspuremd_la_CPPFLAGS = -I sPuReMD/src
sPuReMD_lib_libspuremd_la_LDFLAGS = -version-info 1:0:0
sPuReMD_lib_libspuremd_la_LIBADD = @OMP_LIBS@ @Z_LIBS@ @L_LIBS@ @T_LIBS@ @G_LIBS@

bin_PROGRAMS = sPuReMD/bin/spuremd
sPuReMD_bin_spuremd_SOURCES = sPuReMD/src/driver.c
sPuReMD_bin_spuremd_CFLAGS = @DFT_CFLAGS@ @OMP_CFLAGS@
sPuReMD_bin_spuremd_CPPFLAGS = -I sPuReMD/src
sPuReMD_bin_spuremd_LDADD = sPuReMD/lib/libspuremd.la @OMP_LIBS@ @Z_LIBS@ @L_LIBS@ @T_LIBS@ @G_LIBS@

check_PROGRAMS =
TESTS =

if BUILD_TEST
check_PROGRAMS += sPuReMD/tests/test_spuremd
check_PROGRAMS += sPuReMD/tests/test_lin_alg
check_PROGRAMS += sPuReMD/tests/test_vector
TESTS += $(check_PROGRAMS)

sPuReMD_tests_test_spuremd_SOURCES = sPuReMD/tests/test_spuremd.cpp
sPuReMD_tests_test_spuremd_CPPFLAGS = -I sPuReMD/src $(GTEST_CPPFLAGS)
sPuReMD_tests_test_spuremd_CXXFLAGS = $(GTEST_CXXFLAGS)
sPuReMD_tests_test_spuremd_LDFLAGS = $(GTEST_LDFLAGS) $(GTEST_LIBS)
sPuReMD_tests_test_spuremd_LDADD = sPuReMD/lib/libspuremd.la -lgtest

sPuReMD_tests_test_lin_alg_SOURCES = sPuReMD/tests/test_lin_alg.cpp
sPuReMD_tests_test_lin_alg_CPPFLAGS = -I sPuReMD/src $(GTEST_CPPFLAGS)
sPuReMD_tests_test_lin_alg_CXXFLAGS = $(GTEST_CXXFLAGS)
sPuReMD_tests_test_lin_alg_LDFLAGS = $(GTEST_LDFLAGS) $(GTEST_LIBS)
sPuReMD_tests_test_lin_alg_LDADD = -lgtest

sPuReMD_tests_test_vector_SOURCES = sPuReMD/tests/test_vector.cpp
sPuReMD_tests_test_vector_CPPFLAGS = -I sPuReMD/src $(GTEST_CPPFLAGS)
sPuReMD_tests_test_vector_CXXFLAGS = $(GTEST_CXXFLAGS)
sPuReMD_tests_test_vector_LDFLAGS = $(GTEST_LDFLAGS) $(GTEST_LIBS)
sPuReMD_tests_test_vector_LDADD = -lgtest
endif
endif


dist-hook:
	rm -rf `find $(distdir) -name .git`
