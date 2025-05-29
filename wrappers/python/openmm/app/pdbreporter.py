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
from openmm.unit import angstroms

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
        self._subsetTopology = None


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
        if self._atomSubset is not None:
            if self._subsetTopology is None:
                self._createSubsetTopology(simulation.topology)

            topology = self._subsetTopology

            #PDBFile will convert to angstroms so do it here first instead
            positions = state.getPositions(asNumpy=True).value_in_unit(angstroms)
            positions = [positions[i] for i in self._atomSubset]

        else:
            topology = simulation.topology
            positions = state.getPositions(asNumpy=True)

        if self._nextModel == 0:
            PDBFile.writeHeader(topology, self._out)
            self._topology = topology
            self._nextModel += 1
        PDBFile.writeModel(topology, positions, self._out, self._nextModel)
        self._nextModel += 1
        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def _createSubsetTopology(self, topology):
        """Create a subset of an existing topology.

        Parameters
        ----------
        topology : Topology
            The Topology to create a subset from
        """
        # check atomSubset is valid 
        if len(self._atomSubset) == 0:
            self._out.close()
            raise ValueError('atomSubset cannot be an empty list')
        if not all(a == int(a) for a in self._atomSubset):
            self._out.close()
            raise ValueError('all of the indices in atomSubset must be integers')
        if len(set(self._atomSubset)) != len(self._atomSubset):
            self._out.close()
            raise ValueError('atomSubset must contain unique indices')
        if sorted(self._atomSubset) != self._atomSubset:
            self._out.close()
            raise ValueError('atomSubset must be sorted in ascending order')
        if self._atomSubset[0] < 0:
            self._out.close()
            raise ValueError('The smallest allowed value in atomSubset is zero')
        if self._atomSubset[-1] >= topology.getNumAtoms():
            self._out.close()
            raise ValueError('The maximum allowed value in atomSubset must be less than the total number of particles')
        
        self._subsetTopology = Topology()
        
        # convert to set for fast look up
        atomSubsetSet = set(self._atomSubset)

        # store a map from posIndex to Atom object for when we add the bonds
        indexToAtom = {}

        for chain in topology.chains():
            c = self._subsetTopology.addChain(chain.id)
            for res in chain.residues():
                r = self._subsetTopology.addResidue(res.name, c, res.id, res.insertionCode)
                for atom in res.atoms():
                    if atom.index in atomSubsetSet:
                        indexToAtom[atom.index] = self._subsetTopology.addAtom(atom.name, atom.element, r, atom.id)

        self._subsetTopology.setPeriodicBoxVectors(topology.getPeriodicBoxVectors())

        for bond in topology.bonds():
            if bond[0].index in atomSubsetSet and bond[1].index in atomSubsetSet:
                atom1 = indexToAtom[bond[0].index]
                atom2 = indexToAtom[bond[1].index]
                self._subsetTopology.addBond(atom1, atom2, bond.type, bond.order)
        

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

        if self._atomSubset is not None:
            if self._subsetTopology is None:
                self._createSubsetTopology(simulation.topology)

            topology = self._subsetTopology

            #PDBFile will convert to angstroms so do it here first instead
            positions = state.getPositions(asNumpy=True).value_in_unit(angstroms)
            positions = [positions[i] for i in self._atomSubset]

        else:
            topology = simulation.topology
            positions = state.getPositions(asNumpy=True)

        if self._nextModel == 0:
            PDBxFile.writeHeader(topology, self._out)
            self._nextModel += 1
        PDBxFile.writeModel(topology, positions, self._out, self._nextModel)
        self._nextModel += 1
        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        self._out.close()