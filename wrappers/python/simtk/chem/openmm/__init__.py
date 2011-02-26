#
#
#

"""
Package simtk.chem.openmm

This package wraps the simtk.chem.openmm.openmm module.
When imported, it loads the swig module and then does some magic
to make the POSIX function "dlopen" work on Linux.
It also tries to load any plugin modules it can find.
"""

__author__ = "Randall J. Radmer"
__version__ = "1.0"

import os, sys, glob
if sys.platform == "win32":
    libPrefix=""
    libExt="dll"
elif sys.platform == 'darwin':
    libPrefix="lib"
    libExt="dylib"
else:
    libPrefix="lib"
    libExt="so"
    # The following is an evil incantation that is needed to permit
    # the POSIX "dlopen" function to work.  I do not understand
    # it.  If a better solution is known, please forward to the
    # PyOpenMM code maintainers.
    import ctypes
    flags = sys.getdlopenflags()
    sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL)

import simtk.chem
from simtk.chem.openmm.openmm import *

skipPluginFilenames=["OpenMMFreeEnergy"]
if len(simtk.chem.__path__) > 1:
   raise Exception("more than one item in simtk.chem.openmm.__path__")
pluginPath = os.path.join(simtk.chem.__path__[0],
                          'openmm',
                          'OpenMM',
                          'plugins')
pluginLoadingErrors={}
pluginLoadedLibNames=[]
libFilenames=glob.glob(os.path.join(pluginPath, "*.%s" % (libExt)))

# Load plugins
for filename in libFilenames:
    skipPlugin=False
    for pFilename in skipPluginFilenames:
        fullFilename = "%s%s.%s" % (libPrefix, pFilename, libExt)
        if os.path.split(filename)[-1] == fullFilename:
            skipPlugin=True
            break
    if skipPlugin: continue
    try:
        Platform.loadPluginLibrary(os.path.join(pluginPath, filename))
        pluginLoadedLibNames.append(filename)
    except:
        pluginLoadingErrors[filename]=sys.exc_info()


