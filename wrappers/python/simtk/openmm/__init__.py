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

from simtk.openmm.openmm import *
from simtk.openmm.vec3 import Vec3
from simtk.openmm.mtsintegrator import MTSIntegrator
from simtk.openmm.amd import AMDIntegrator, AMDForceGroupIntegrator, DualAMDIntegrator

pluginLoadedLibNames = Platform.loadPluginsFromDirectory(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib', 'plugins'))

__version__ = Platform.getOpenMMVersion()
