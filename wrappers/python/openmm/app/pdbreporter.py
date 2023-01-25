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
__author__ = "Peter Eastman"
__version__ = "1.0"

from openmm.app import PDBFile, PDBxFile, Topology
from openmm.unit import nanometers, Quantity

class PDBReporter(object):
    """PDBReporter outputs a series of frames from a Simulation to a PDB file.

    To use it, create a PDBReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, enforcePeriodicBox=None, atomSubset=None):
        """Create a PDBReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        atomSubset: list
            Atom indices (zero indexed) of the particles to output. if None (the default), all particles will be output.
        """
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._out = open(file, 'w')
        self._topology = None
        self._nextModel = 0
        self._atomSubset = atomSubset

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

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if self._atomSubset is not None:
            if not all(a==int(a) for a in self._atomSubset):
                raise ValueError('all of the indices in atomSubset must be integers')
            if min(self._atomSubset) < 0:
                raise ValueError('The smallest allowed value in atomSubset is zero')
            if max(self._atomSubset) >= simulation.topology.getNumAtoms():
                raise ValueError('The maximum allowed value in atomSubset must be less than the total number of particles')
            if len(set(self._atomSubset)) != len(self._atomSubset):
                raise ValueError('atomSubset must contain unique indices')

            topology = _subsetTopology(simulation.topology, self._atomSubset)
            positions = _subsetPositions(state.getPositions(), self._atomSubset)
        else:
            topology = simulation.topology
            positions = state.getPositions()

        if self._nextModel == 0:
            PDBFile.writeHeader(topology, self._out)
            self._topology = topology
            self._nextModel += 1
        PDBFile.writeModel(topology, positions, self._out, self._nextModel)
        self._nextModel += 1
        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        if self._topology is not None:
            PDBFile.writeFooter(self._topology, self._out)
        self._out.close()

class PDBxReporter(PDBReporter):
    """PDBxReporter outputs a series of frames from a Simulation to a PDBx/mmCIF file.

    To use it, create a PDBxReporter, then add it to the Simulation's list of reporters.
    """

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if self._nextModel == 0:
            PDBxFile.writeHeader(simulation.topology, self._out)
            self._nextModel += 1
        PDBxFile.writeModel(simulation.topology, state.getPositions(), self._out, self._nextModel)
        self._nextModel += 1
        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        self._out.close()

def _subsetPositions(positions, atomSubset):
    """Create a subset of the positions

    Parameters
    ----------
    positions : list
        The positions
    atomSubset : list
        The list of atomic indices in the subset

    Returns
    -------
    subsetPositions : list
        A subset of the input positions that only contains the atoms
        specified in atomSubset.
    """

    return Quantity([positions[i].value_in_unit(nanometers) for i in atomSubset], unit=nanometers)

    
def _subsetTopology(topology, atomSubset):
    """Create a subset of an existing topology.

    Parameters
    ----------
    topology : Topology
        The Topology to create a subset from
    atomSubset : list
        The list of atomic indices in the subset

    Returns
    -------
    subsetTopology : Topology
        A new Topology copied from the input topology that only contains the atoms
        specified in atomSubset.
    """
    subsetTopology = Topology()

    posIndex = 0
    for chain in topology.chains():
        c = subsetTopology.addChain(chain.id)
        residues = list(chain.residues())
        for res in residues:
            r = subsetTopology.addResidue(res.name,c,res.id,res.insertionCode)
            for atom in res.atoms():
                    if posIndex in atomSubset:
                        atom = subsetTopology.addAtom(atom.name, atom.element, r, atom.id)
                    posIndex += 1

    return subsetTopology