"""
dcdreporter.py: Outputs simulation trajectories in DCD format

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
__author__ = "Peter Eastman"
__version__ = "1.0"

from openmm.app import DCDFile
from openmm.unit import nanometer

class DCDReporter(object):
    """DCDReporter outputs a series of frames from a Simulation to a DCD file.

    To use it, create a DCDReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, append=False, enforcePeriodicBox=None):
        """Create a DCDReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        append : bool=False
            If True, open an existing DCD file to append to.  If False, create a new file.
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        """
        self._reportInterval = reportInterval
        self._append = append
        self._enforcePeriodicBox = enforcePeriodicBox
        if append:
            mode = 'r+b'
        else:
            mode = 'wb'
        self._out = open(file, mode)
        self._dcd = None

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        dict
            A dictionary describing the required information for the next report
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return {'steps':steps, 'periodic':self._enforcePeriodicBox, 'include':['positions']}

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """

        if self._dcd is None:
            self._dcd = DCDFile(
                self._out, simulation.topology, simulation.integrator.getStepSize(),
                simulation.currentStep, self._reportInterval, self._append
            )
        self._dcd.writeModel(state.getPositions(asNumpy=True), periodicBoxVectors=state.getPeriodicBoxVectors())

    def __del__(self):
        self._out.close()
