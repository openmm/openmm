"""
xtcreporter.py: Outputs simulation trajectories in XTC format

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2024 Stanford University and the Authors.
Authors: Raul P. Pelaez
Contributors: Peter Eastman

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
__author__ = "Raul P. Pelaez"
from openmm.app import Topology, XTCFile

class XTCReporter(object):
    """XTCReporter outputs a series of frames from a Simulation to a XTC file.

    To use it, create a XTCReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, append=False, enforcePeriodicBox=None, atomSubset=None):
        """Create a XTCReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        append : bool=False
            If True, open an existing XTC file to append to.  If False, create a new file.
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        atomSubset: list
            Atom indices (zero indexed) of the particles to output.  If None (the default), all particles will be output.
        """
        self._reportInterval = reportInterval
        self._append = append
        self._enforcePeriodicBox = enforcePeriodicBox
        self._atomSubset = atomSubset
        self._fileName = file
        self._xtc = None
        if not append:
            open(file, 'wb').close()

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

        if self._xtc is None:
            if self._atomSubset is None:
                topology = simulation.topology
            else:
                topology = Topology()
                topology.setPeriodicBoxVectors(simulation.topology.getPeriodicBoxVectors())
                atoms = list(simulation.topology.atoms())
                chain = topology.addChain()
                residue = topology.addResidue('', chain)
                for i in self._atomSubset:
                    topology.addAtom(atoms[i].name, atoms[i].element, residue)
            self._xtc = XTCFile(
                self._fileName,
                topology,
                simulation.integrator.getStepSize(),
                self._reportInterval,
                self._reportInterval,
                self._append,
            )
        positions = state.getPositions(asNumpy=True)
        if self._atomSubset is not None:
            positions = [positions[i] for i in self._atomSubset]
        self._xtc.writeModel(positions, periodicBoxVectors=state.getPeriodicBoxVectors())
