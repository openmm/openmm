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
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm.app import Topology
import simtk.openmm as mm
from simtk.openmm.app import DCDFile
from simtk.unit import nanometer
    
class DCDReporter(object):
    """DCDReporter outputs a series of frames from a Simulation to a DCD file.
    
    To use it, create a DCDReporter, then add it to the Simulation's list of reporters.
    """
    
    def __init__(self, file, reportInterval, atomIndices=None):
        """Create a DCDReporter.
    
        Parameters:
         - file (string) The file to write to
         - reportInterval (int) The interval (in time steps) at which to write frames
         - atomIndices (list) The indicies (zero-based) of the atoms you wish to
           write to disk. If not supplied, all of the atoms will be written
           to the file (default)
        """
        self._reportInterval = reportInterval
        self._atomIndices = atomIndices
        self._is_initialized = False
        self._out = open(file, 'wb')
        self._dcd = None

    def _initialize(self, simulation):
        """Delayed initialization

        This is called before the first report is written, once we know the
        simulation we'll be reporting on. It sets up the file-like object
        """
        if self._atomIndices is not None:
            if min(self._atomIndices) < 0:
                raise ValueError('%s is an invalid index. it\'s less than '
                                 'zero!' % min(self._atomIndices))
            if max(self._atomIndices) > max(atom.index for atom in simulation.topology.atoms()):
                raise ValueError('%s is an invalid index. it\'s greater than '
                                 'the index of the last atom.' % max(self._atomIndices))
            topology = simulation.topology.subset(self._atomIndices)

        else:
            topology = simulation.topology

        self._dcd = DCDFile(self._out, topology, simulation.integrator.getStepSize(), 0, self._reportInterval)    
        self._is_initialized = True

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.
        
        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
        Returns: A five element tuple.  The first element is the number of steps until the
        next report.  The remaining elements specify whether that report will require
        positions, velocities, forces, and energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False)
    
    def report(self, simulation, state):
        """Generate a report.
        
        Parameters:
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if not self._is_initialized:
            self._initialize(simulation)
        a, b, c = state.getPeriodicBoxVectors()

        positions = state.getPositions()
        if self._atomIndices is not None:
            positions = [positions[i].value_in_unit(nanometer) for i in self._atomIndices] * nanometer

        self._dcd.writeModel(positions, mm.Vec3(a[0].value_in_unit(nanometer), b[1].value_in_unit(nanometer), c[2].value_in_unit(nanometer))*nanometer)

    def __del__(self):
        self._out.close()
