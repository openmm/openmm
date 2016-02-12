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
