"""
pdbfile.py: Used for loading PDB files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2014-2015 Stanford University and the Authors.
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

import os
import sys
import math
from simtk.openmm import Vec3
from simtk.openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
from simtk.openmm.app.internal.unitcell import computePeriodicBoxVectors
from simtk.openmm.app import Topology
from simtk.unit import nanometers, angstroms, is_quantity, norm, Quantity
import element as elem
try:
    import numpy
except:
    pass

class PDBxFile(object):
    """PDBxFile parses a PDBx/mmCIF file and constructs a Topology and a set of atom positions from it."""

    def __init__(self, file):
        """Load a PDBx/mmCIF file.

        The atom positions and Topology can be retrieved by calling getPositions() and getTopology().

        Parameters:
         - file (string) the name of the file to load.  Alternatively you can pass an open file object.
        """
        top = Topology()
        ## The Topology read from the PDBx/mmCIF file
        self.topology = top
        self._positions = []

        # Load the file.

        inputFile = file
        if isinstance(file, str):
            inputFile = open(file)
        reader = PdbxReader(inputFile)
        data = []
        reader.read(data)
        block = data[0]

        # Build the topology.

        atomData = block.getObj('atom_site')
        atomNameCol = atomData.getAttributeIndex('label_atom_id')
        atomIdCol = atomData.getAttributeIndex('id')
        resNameCol = atomData.getAttributeIndex('label_comp_id')
        resIdCol = atomData.getAttributeIndex('label_seq_id')
        resNumCol = atomData.getAttributeIndex('auth_seq_id')
        asymIdCol = atomData.getAttributeIndex('label_asym_id')
        chainIdCol = atomData.getAttributeIndex('label_entity_id')
        elementCol = atomData.getAttributeIndex('type_symbol')
        altIdCol = atomData.getAttributeIndex('label_alt_id')
        modelCol = atomData.getAttributeIndex('pdbx_PDB_model_num')
        xCol = atomData.getAttributeIndex('Cartn_x')
        yCol = atomData.getAttributeIndex('Cartn_y')
        zCol = atomData.getAttributeIndex('Cartn_z')
        lastChainId = None
        lastResId = None
        lastAsymId = None
        atomTable = {}
        models = []
        for row in atomData.getRowList():
            atomKey = ((row[resIdCol], row[asymIdCol], row[atomNameCol]))
            model = ('1' if modelCol == -1 else row[modelCol])
            if model not in models:
                models.append(model)
                self._positions.append([])
            modelIndex = models.index(model)
            if row[altIdCol] != '.' and atomKey in atomTable and len(self._positions[modelIndex]) > atomTable[atomKey].index:
                # This row is an alternate position for an existing atom, so ignore it.

                continue
            if modelIndex == 0:
                # This row defines a new atom.

                if lastChainId != row[chainIdCol]:
                    # The start of a new chain.
                    chain = top.addChain(row[chainIdCol])
                    lastChainId = row[chainIdCol]
                    lastResId = None
                    lastAsymId = None
                if lastResId != row[resIdCol] or lastAsymId != row[asymIdCol]:
                    # The start of a new residue.
                    res = top.addResidue(row[resNameCol], chain, None if resNumCol == -1 else row[resNumCol])
                    lastResId = row[resIdCol]
                    if lastResId == '.':
                        lastResId = None
                    lastAsymId = row[asymIdCol]
                element = None
                try:
                    element = elem.get_by_symbol(row[elementCol])
                except KeyError:
                    pass
                atom = top.addAtom(row[atomNameCol], element, res, row[atomIdCol])
                atomTable[atomKey] = atom
            else:
                # This row defines coordinates for an existing atom in one of the later models.

                try:
                    atom = atomTable[atomKey]
                except KeyError:
                    raise ValueError('Unknown atom %s in residue %s %s for model %s' % (row[atomNameCol], row[resNameCol], row[resIdCol], model))
                if atom.index != len(self._positions[modelIndex]):
                    raise ValueError('Atom %s for model %s does not match the order of atoms for model %s' % (row[atomIdCol], model, models[0]))
            self._positions[modelIndex].append(Vec3(float(row[xCol]), float(row[yCol]), float(row[zCol]))*0.1)
        for i in range(len(self._positions)):
            self._positions[i] = self._positions[i]*nanometers
        ## The atom positions read from the PDBx/mmCIF file.  If the file contains multiple frames, these are the positions in the first frame.
        self.positions = self._positions[0]
        self.topology.createStandardBonds()
        self._numpyPositions = None

        # Record unit cell information, if present.

        cell = block.getObj('cell')
        if cell is not None and cell.getRowCount() > 0:
            row = cell.getRow(0)
            (a, b, c) = [float(row[cell.getAttributeIndex(attribute)])*0.1 for attribute in ('length_a', 'length_b', 'length_c')]
            (alpha, beta, gamma) = [float(row[cell.getAttributeIndex(attribute)])*math.pi/180.0 for attribute in ('angle_alpha', 'angle_beta', 'angle_gamma')]
            self.topology.setPeriodicBoxVectors(computePeriodicBoxVectors(a, b, c, alpha, beta, gamma))

        # Add bonds based on struct_conn records.

        connectData = block.getObj('struct_conn')
        if connectData is not None:
            res1Col = connectData.getAttributeIndex('ptnr1_label_seq_id')
            res2Col = connectData.getAttributeIndex('ptnr2_label_seq_id')
            atom1Col = connectData.getAttributeIndex('ptnr1_label_atom_id')
            atom2Col = connectData.getAttributeIndex('ptnr2_label_atom_id')
            asym1Col = connectData.getAttributeIndex('ptnr1_label_asym_id')
            asym2Col = connectData.getAttributeIndex('ptnr2_label_asym_id')
            typeCol = connectData.getAttributeIndex('conn_type_id')
            connectBonds = []
            for row in connectData.getRowList():
                type = row[typeCol][:6]
                if type in ('covale', 'disulf', 'modres'):
                    key1 = (row[res1Col], row[asym1Col], row[atom1Col])
                    key2 = (row[res2Col], row[asym2Col], row[atom2Col])
                    if key1 in atomTable and key2 in atomTable:
                        connectBonds.append((atomTable[key1], atomTable[key2]))
            if len(connectBonds) > 0:
                # Only add bonds that don't already exist.
                existingBonds = set(top.bonds())
                for bond in connectBonds:
                    if bond not in existingBonds and (bond[1], bond[0]) not in existingBonds:
                        top.addBond(bond[0], bond[1])
                        existingBonds.add(bond)

    def getTopology(self):
        """Get the Topology of the model."""
        return self.topology

    def getNumFrames(self):
        """Get the number of frames stored in the file."""
        return len(self._positions)

    def getPositions(self, asNumpy=False, frame=0):
        """Get the atomic positions.

        Parameters:
         - asNumpy (boolean=False) if true, the values are returned as a numpy array instead of a list of Vec3s
         - frame (int=0) the index of the frame for which to get positions
         """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = [None]*len(self._positions)
            if self._numpyPositions[frame] is None:
                self._numpyPositions[frame] = Quantity(numpy.array(self._positions[frame].value_in_unit(nanometers)), nanometers)
            return self._numpyPositions[frame]
        return self._positions[frame]
