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

LIB_DIR=$(OpenMM_INSTALL_DIR)/lib
INCLUDE_DIR=$(OpenMM_INSTALL_DIR)/include
LIBS= -lOpenMM$(DEBUG)

ALL_CPP_EXAMPLES = HelloArgon HelloSodiumChloride HelloEthane HelloWaterBox
ALL_C_EXAMPLES   = HelloArgonInC HelloSodiumChlorideInC
ALL_F95_EXAMPLES = HelloArgonInFortran HelloSodiumChlorideInFortran

default: HelloArgon

all : $(ALL_PROGS)

# Treat all .cpp source files the same way
.cpp : 
	g++ $(CFLAGS) $< -I$(INCLUDE_DIR) -L$(LIB_DIR) $(LIBS) -o $*

.c : OpenMM_CWrapper
	g++ $(CFLAGS) $< OpenMM_CWrapper.o -L$(LIB_DIR) $(LIBS) -o $*

OpenMM_CWrapper: OpenMM_CWrapper.cpp
	g++ -c $(CFLAGS) $< -I$(INCLUDE_DIR) -o OpenMM_CWrapper.o

clean : 
	rm $(ALL_PROGS)

