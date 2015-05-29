"""
Package simtk.openmm

This package wraps the simtk.openmm.openmm module.
When imported, it loads the swig module and then does some magic
to make the POSIX function "dlopen" work on Linux.
It also tries to load any plugin modules it can find.
"""

__author__ = "Randall J. Radmer"

import os, os.path
import sys
from simtk.openmm import version

if sys.platform == 'win32':
    _path = os.environ['PATH']
    os.environ['PATH'] = '%(lib)s;%(lib)s\plugins;%(path)s' % {
        'lib': version.openmm_library_path, 'path': _path}

from simtk.openmm.openmm import *
from simtk.openmm.vec3 import Vec3
from simtk.openmm.mtsintegrator import MTSIntegrator
from simtk.openmm.amd import AMDIntegrator, AMDForceGroupIntegrator, DualAMDIntegrator

if os.getenv('OPENMM_PLUGIN_DIR') is None and os.path.isdir(version.openmm_library_path):
    pluginLoadedLibNames = Platform.loadPluginsFromDirectory(os.path.join(version.openmm_library_path, 'plugins'))
else:
    pluginLoadedLibNames = Platform.loadPluginsFromDirectory(Platform.getDefaultPluginsDirectory())

if sys.platform == 'win32':
    os.environ['PATH'] = _path
    del _path
__version__ = Platform.getOpenMMVersion()
