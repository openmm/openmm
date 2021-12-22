"""OpenMM is a toolkit for molecular simulation. It can be used either as a
stand-alone application for running simulations, or as a library you call
from your own code. It provides a combination of extreme flexibility
(through custom forces and integrators), openness, and high performance
(especially on recent GPUs) that make it truly unique among simulation codes.
"""
from __future__ import absolute_import
__author__ = "Peter Eastman"

import os, os.path
import sys
from . import version

if sys.platform == 'win32':
    _path = os.environ['PATH']
    os.environ['PATH'] = '%(dll_path)s;%(lib_path)s\plugins;%(path)s' % {
        'lib_path': version.openmm_library_path,
        'dll_path': version.openmm_dll_path,
        'path': _path}

    if sys.version_info[:2] >= (3, 8):
        os.add_dll_directory(version.openmm_dll_path)
        os.add_dll_directory(os.path.join(version.openmm_library_path, plugins))

from openmm.openmm import *
from openmm.vec3 import Vec3
from openmm.mtsintegrator import MTSIntegrator, MTSLangevinIntegrator
from openmm.amd import AMDIntegrator, AMDForceGroupIntegrator, DualAMDIntegrator

if os.getenv('OPENMM_PLUGIN_DIR') is None and os.path.isdir(version.openmm_library_path):
    pluginLoadedLibNames = Platform.loadPluginsFromDirectory(os.path.join(version.openmm_library_path, 'plugins'))
else:
    pluginLoadedLibNames = Platform.loadPluginsFromDirectory(Platform.getDefaultPluginsDirectory())

if sys.platform == 'win32':
    os.environ['PATH'] = _path
    del _path
__version__ = Platform.getOpenMMVersion()

class OpenMMException(Exception):
    """This is the class used for all exceptions thrown by the C++ library."""
    pass
