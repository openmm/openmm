#
#
#

"""
Package simtk.openmm

This package wraps the simtk.openmm.openmm module.
When imported, it loads the swig module and then does some magic
to make the POSIX function "dlopen" work on Linux.
It also tries to load any plugin modules it can find.
"""

__author__ = "Randall J. Radmer"

import os, sys, glob, os.path
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


from simtk.openmm.openmm import *
from simtk.openmm.vec3 import Vec3
from simtk.openmm import version
if os.getenv('OPENMM_PLUGIN_DIR') is None and os.path.isdir(version.openmm_library_path):
    pluginLoadedLibNames = Platform.loadPluginsFromDirectory(os.path.join(version.openmm_library_path, 'plugins'))
else:
    pluginLoadedLibNames = Platform.loadPluginsFromDirectory(Platform.getDefaultPluginsDirectory())
__version__ = Platform.getOpenMMVersion()
