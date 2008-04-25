# ---------------------------------------------------------------------------------------

include ../makeConfig/Top.mk

# ---------------------------------------------------------------------------------------

# target is gmxGbsa library
# ld reference is in ../gromacs/vc7.mk
# include file reference is in ../gromacs/gmxgpu.h

STATIC_LIBRARY      := gmxGbsa

# ---------------------------------------------------------------------------------------

BaseFiles           := $(GbsaBaseFiles)
IncludeDirectories  := $(SimTkDirectory) $(GromacsIncludeDirectory)

# ---------------------------------------------------------------------------------------

include $(MakeDirectory)/BaseMake.mk

# ---------------------------------------------------------------------------------------
