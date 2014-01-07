"""
Package simtk.openmm

This package wraps the simtk.openmm.openmm module.
When imported, it loads the swig module and then does some magic
to make the POSIX function "dlopen" work on Linux.
It also tries to load any plugin modules it can find.
"""

__author__ = "Randall J. Radmer"

import os, os.path
from simtk.openmm.openmm import *
from simtk.openmm.vec3 import Vec3
from simtk.openmm import version
if os.getenv('OPENMM_PLUGIN_DIR') is None and os.path.isdir(version.openmm_library_path):
    pluginLoadedLibNames = Platform.loadPluginsFromDirectory(os.path.join(version.openmm_library_path, 'plugins'))
else:
    pluginLoadedLibNames = Platform.loadPluginsFromDirectory(Platform.getDefaultPluginsDirectory())
__version__ = Platform.getOpenMMVersion()
