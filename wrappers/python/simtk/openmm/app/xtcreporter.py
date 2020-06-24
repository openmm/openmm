"""
pdbreporter.py: Outputs simulation trajectories in PDB format

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import absolute_import
__author__ = "Jade Shi"
__version__ = "1.0"

import numpy as np
import simtk.openmm as mm
from simtk.unit import *
import mdtraj as mdt
import os

class XTCReporter(object):
    """XTCReporter outputs a series of frames from a Simulation to a XTC file.

    To use it, create a XTCReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, enforcePeriodicBox=False):
        """Create a XTCReporter.

        Parameters
        ----------
        simulation : Simuluation
            The Simulation to generate a report for
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  This is False by default, as not every simulation employs PBC.
        """
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._out = file
        self._topology = None
        self._nextModel = 0
        
    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A six element tuple. The first element is the number of steps
            until the next report. The next four elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.  The final element specifies whether
            positions should be wrapped to lie in a single periodic box.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False, self._enforcePeriodicBox)

    def report(self, simulation, state=None):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : NOT USED, this is only here to be consistent with the PDBReporter's report function.
        """
        if self._nextModel == 0:
            self._topology = mdt.Trajectory(xyz=simulation.context.getState(getPositions=True).getPositions(asNumpy=True) / nanometer, \
                                        topology=mdt.Topology().from_openmm(simulation.topology))
            self._nextModel += 1
        self._topology.xyz = simulation.context.getState(getPositions=True, enforcePeriodicBox=self._enforcePeriodicBox).getPositions(asNumpy=True) / nanometer
        self._topology.save_xtc('%s.append'%(self._out))
        os.system('cat %s.append >> %s'%(self._out, self._out))
        self._nextModel += 1
        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()
