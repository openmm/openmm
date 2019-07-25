# ----------------------------------------------------------------------
# Makefile for OpenMM Preview Release 4 workshop "hello world" examples.
# August 18, 2009
# See https://simtk.org/home/openmm.
# ----------------------------------------------------------------------
# This assumes you have gcc compilers for whatever language you are
# using: g++ for C++ and C, gfortran for Fortran 95.
# 
# For the C and Fortran examples, we're depending on your version of
# OpenMM to have been built with the automatically-generated API
# wrappers.
#
# This has had only minimal testing, although it has been known to
# work. It is likely to work fine for C and C++. For Fortran, you
# may need to add some of the C/C++ libraries: 
#    -lc -lm -lstdc++ (or -lstdc++.6) -lgcc -lgcc_s
# but this wasn't required for these examples on Centos 5.2 using
# gcc 4.1.2.

# Check whether this is the right capitalization for your install directory.
OpenMM_INSTALL_DIR=/usr/local/openmm
CFLAGS = -g 
FFLAGS = -g -ffree-line-length-none

# Extra libraries required when gfortran links with a C++ module.
# If this doesn't work, look in /usr/lib to
# see what versions of libstdc++.so you have and try different
# versions.
FCPPLIBS = -lstdc++

LIB_DIR=$(OpenMM_INSTALL_DIR)/lib
INCLUDE_DIR=$(OpenMM_INSTALL_DIR)/include
LIBS= -lOpenMM

ALL_CPP_EXAMPLES = HelloArgon HelloSodiumChloride HelloEthane HelloWaterBox
ALL_C_EXAMPLES   = HelloArgonInC HelloSodiumChlorideInC
ALL_F95_EXAMPLES = HelloArgonInFortran HelloSodiumChlorideInFortran

ALL_PROGS = $(ALL_CPP_EXAMPLES) $(ALL_C_EXAMPLES) $(ALL_F95_EXAMPLES)

default: HelloArgon

all : $(ALL_PROGS)

# Treat all .cpp source files the same way.
.cpp : 
	g++ $(CFLAGS) -I$(INCLUDE_DIR) $< -L$(LIB_DIR) $(LIBS) -o $*

HelloArgonInC: HelloArgonInC.c 
	g++ $(CFLAGS) -I$(INCLUDE_DIR) HelloArgonInC.c \
        -L$(LIB_DIR) $(LIBS) -o HelloArgonInC

HelloSodiumChlorideInC: HelloSodiumChlorideInC.c 
	g++ $(CFLAGS) -I$(INCLUDE_DIR) HelloSodiumChlorideInC.c \
        -L$(LIB_DIR) $(LIBS) -o HelloSodiumChlorideInC

HelloArgonInFortran: HelloArgonInFortran.f90 openmm.mod 
	gfortran $(FFLAGS) HelloArgonInFortran.f90 \
	    -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) $(FCPPLIBS) -o HelloArgonInFortran

HelloSodiumChlorideInFortran: HelloSodiumChlorideInFortran.f90 openmm.mod 
	gfortran $(FFLAGS) HelloSodiumChlorideInFortran.f90 \
	    -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) $(FCPPLIBS) -o HelloSodiumChlorideInFortran


# Build Fortran 95 Module file
openmm.mod:
	gfortran -c $(FFLAGS) $(INCLUDE_DIR)/OpenMMFortranModule.f90

clean : 
	rm $(ALL_PROGS) *.o *.mod *.obj *.exe

