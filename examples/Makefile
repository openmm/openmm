# Uncomment DEBUG to use OpenMM libraries instead of release libraries
# DEBUG=_d

# Default install directory
#
# Linux
# CFLAGS = -g -m32
#
# Mac
# CFLAGS = -g -m32 -bind_at_load

# Check whether this is the right capitalization for your install directory.
OpenMM_INSTALL_DIR=/usr/local/OpenMM
CFLAGS = -g 
FFLAGS = -g

LIB_DIR=$(OpenMM_INSTALL_DIR)/lib
INCLUDE_DIR=$(OpenMM_INSTALL_DIR)/include
# assume local directory for C and Fortran wrappers
WRAPPER_DIR=.
LIBS= -lOpenMM$(DEBUG)

ALL_CPP_EXAMPLES = HelloArgon HelloSodiumChloride HelloEthane HelloWaterBox
ALL_C_EXAMPLES   = HelloArgonInC HelloSodiumChlorideInC
ALL_F95_EXAMPLES = HelloArgonInFortran HelloSodiumChlorideInFortran

ALL_PROGS = $(ALL_CPP_EXAMPLES) $(ALL_C_EXAMPLES) $(ALL_F95_EXAMPLES)

default: HelloArgon

all : $(ALL_PROGS)

# Treat all .cpp source files the same way (except the one that 
# implements the C Wrappers).
.cpp : 
	g++ $(CFLAGS) $< -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) -o $*

HelloArgonInC: HelloArgonInC.c OpenMM_CWrapper.o
	g++ $(CFLAGS) -I$(WRAPPER_DIR) HelloArgonInC.c OpenMM_CWrapper.o \
        -L$(LIB_DIR) $(LIBS) -o HelloArgonInC

HelloSodiumChlorideInC: HelloSodiumChlorideInC.c OpenMM_CWrapper.o
	g++ $(CFLAGS) -I$(WRAPPER_DIR) HelloSodiumChlorideInC.c OpenMM_CWrapper.o \
        -L$(LIB_DIR) $(LIBS) -o HelloSodiumChlorideInC

HelloArgonInFortran: HelloArgonInFortran.f90 openmm.mod
	gfortran $(FFLAGS) HelloArgonInFortran.f90 OpenMM_CWrapper.o \
	    -L$(LIB_DIR) $(LIBS) -o HelloArgonInFortran

HelloSodiumChlorideInFortran: HelloSodiumChlorideInFortran.f90 openmm.mod
	gfortran $(FFLAGS) HelloSodiumChlorideInFortran.f90 OpenMM_CWrapper.o \
	    -L$(LIB_DIR) $(LIBS) -o HelloSodiumChlorideInFortran

# Build C Wrappers (Fortran depends on this too)
OpenMM_CWrapper.o: 
	g++ -c $(CFLAGS) -I$(WRAPPER_DIR) -I$(INCLUDE_DIR) \
         $(WRAPPER_DIR)/OpenMM_CWrapper.cpp -o OpenMM_CWrapper.o

# Build Fortran 95 Module file
openmm.mod:
	gfortran -c $(FFLAGS) $(WRAPPER_DIR)/OpenMM_Module.f90

clean : 
	rm $(ALL_PROGS) *.o *.mod

