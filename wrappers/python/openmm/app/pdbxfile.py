"""
pdbxfile.py: Used for loading PDBx/mmCIF files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2015-2023 Stanford University and the Authors.
Authors: Peter Eastman
Contributors: Jason Swails

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
from __future__ import division, absolute_import, print_function

__author__ = "Peter Eastman"
__version__ = "2.0"

import sys
import math
from openmm import Vec3, Platform
from datetime import date
from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
from openmm.app.internal.unitcell import computePeriodicBoxVectors, computeLengthsAndAngles
from openmm.app import Topology, PDBFile
from openmm.unit import nanometers, angstroms, is_quantity, norm, Quantity
from . import element as elem
try:
    import numpy
except:
    pass

class PDBxFile(object):
    """PDBxFile parses a PDBx/mmCIF file and constructs a Topology and a set of atom positions from it."""

    def __init__(self, file):
        """Load a PDBx/mmCIF file.

        The atom positions and Topology can be retrieved by calling
        getPositions() and getTopology().

        Parameters
        ----------
        file : string
            the name of the file to load.  Alternatively you can pass an open
            file object.
        """
        top = Topology()
        ## The Topology read from the PDBx/mmCIF file
        self.topology = top
        self._positions = []
        PDBFile._loadNameReplacementTables()

        # Load the file.

        inputFile = file
        ownHandle = False
        if isinstance(file, str):
            inputFile = open(file)
            ownHandle = True
        reader = PdbxReader(inputFile)
        data = []
        reader.read(data)
        if ownHandle:
            inputFile.close()
        block = data[0]

        # Build the topology.

        atomData = block.getObj('atom_site')
        atomNameCol = atomData.getAttributeIndex('auth_atom_id')
        if atomNameCol == -1:
            atomNameCol = atomData.getAttributeIndex('label_atom_id')
        atomIdCol = atomData.getAttributeIndex('id')
        resNameCol = atomData.getAttributeIndex('auth_comp_id')
        if resNameCol == -1:
            resNameCol = atomData.getAttributeIndex('label_comp_id')
        resNumCol = atomData.getAttributeIndex('auth_seq_id')
        if resNumCol == -1:
            resNumCol = atomData.getAttributeIndex('label_seq_id')
        resInsertionCol = atomData.getAttributeIndex('pdbx_PDB_ins_code')
        chainIdCol = atomData.getAttributeIndex('auth_asym_id')
        if chainIdCol == -1:
            chainIdCol = atomData.getAttributeIndex('label_asym_id')
            altChainIdCol = -1
        else:
            altChainIdCol = atomData.getAttributeIndex('label_asym_id')
        if altChainIdCol != -1:
            # Figure out which column is best to use for chain IDs.
            
            idSet = set(row[chainIdCol] for row in atomData.getRowList())
            altIdSet = set(row[altChainIdCol] for row in atomData.getRowList())
            if len(altIdSet) > len(idSet):
                chainIdCol, altChainIdCol = (altChainIdCol, chainIdCol)
        elementCol = atomData.getAttributeIndex('type_symbol')
        altIdCol = atomData.getAttributeIndex('label_alt_id')
        modelCol = atomData.getAttributeIndex('pdbx_PDB_model_num')
        xCol = atomData.getAttributeIndex('Cartn_x')
        yCol = atomData.getAttributeIndex('Cartn_y')
        zCol = atomData.getAttributeIndex('Cartn_z')
        lastChainId = None
        lastAltChainId = None
        lastResId = None
        lastInsertionCode = ''
        atomTable = {}
        atomsInResidue = set()
        models = []
        for row in atomData.getRowList():
            atomKey = ((row[resNumCol], row[chainIdCol], row[atomNameCol]))
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

                if resInsertionCol == -1:
                    insertionCode = ''
                else:
                    insertionCode = row[resInsertionCol]
                if insertionCode in ('.', '?'):
                    insertionCode = ''
                if lastChainId != row[chainIdCol] or (altChainIdCol != -1 and lastAltChainId != row[altChainIdCol]):
                    # The start of a new chain.
                    chain = top.addChain(row[chainIdCol])
                    lastChainId = row[chainIdCol]
                    lastResId = None
                    if altChainIdCol != -1:
                        lastAltChainId = row[altChainIdCol]
                if lastResId != row[resNumCol] or lastChainId != row[chainIdCol] or lastInsertionCode != insertionCode or (lastResId == '.' and row[atomNameCol] in atomsInResidue):
                    # The start of a new residue.
                    resId = (None if resNumCol == -1 else row[resNumCol])
                    resIC = insertionCode
                    resName = row[resNameCol]
                    if resName in PDBFile._residueNameReplacements:
                        resName = PDBFile._residueNameReplacements[resName]
                    res = top.addResidue(resName, chain, resId, resIC)
                    if resName in PDBFile._atomNameReplacements:
                        atomReplacements = PDBFile._atomNameReplacements[resName]
                    else:
                        atomReplacements = {}
                    lastResId = row[resNumCol]
                    lastInsertionCode = insertionCode
                    atomsInResidue.clear()
                element = None
                try:
                    element = elem.get_by_symbol(row[elementCol])
                except KeyError:
                    pass
                atomName = row[atomNameCol]
                if atomName in atomReplacements:
                    atomName = atomReplacements[atomName]
                atom = top.addAtom(atomName, element, res, row[atomIdCol])
                atomTable[atomKey] = atom
                atomsInResidue.add(atomName)
            else:
                # This row defines coordinates for an existing atom in one of the later models.

                try:
                    atom = atomTable[atomKey]
                except KeyError:
                    raise ValueError('Unknown atom %s in residue %s %s for model %s' % (row[atomNameCol], row[resNameCol], row[resNumCol], model))
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

        Parameters
        ----------
        asNumpy : bool=False
            if true, the values are returned as a numpy array instead of a list
            of Vec3s
        frame : int=0
            the index of the frame for which to get positions
        """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = [None]*len(self._positions)
            if self._numpyPositions[frame] is None:
                self._numpyPositions[frame] = Quantity(numpy.array(self._positions[frame].value_in_unit(nanometers)), nanometers)
            return self._numpyPositions[frame]
        return self._positions[frame]

    @staticmethod
    def writeFile(topology, positions, file=sys.stdout, keepIds=False,
                  entry=None):
        """Write a PDBx/mmCIF file containing a single model.

        Parameters
        ----------
        topology : Topology
            The Topology defining the model to write
        positions : list
            The list of atomic positions to write
        file : string or file
            the name of the file to write.  Alternatively you can pass an open file object.
        keepIds : bool=False
            If True, keep the residue and chain IDs specified in the Topology
            rather than generating new ones.  Warning: It is up to the caller to
            make sure these are valid IDs that satisfy the requirements of the
            PDBx/mmCIF format.  Otherwise, the output file will be invalid.
        entry : str=None
            The entry ID to assign to the CIF file
        """
        if isinstance(file, str):
            with open(file, 'w') as output:
                PDBxFile.writeFile(topology, positions, output, keepIds, entry)
        else:
            PDBxFile.writeHeader(topology, file, entry, keepIds)
            PDBxFile.writeModel(topology, positions, file, keepIds=keepIds)

    @staticmethod
    def writeHeader(topology, file=sys.stdout, entry=None, keepIds=False):
        """Write out the header for a PDBx/mmCIF file.

        Parameters
        ----------
        topology : Topology
            The Topology defining the molecular system being written
        file : file=stdout
            A file to write the file to
        entry : str=None
            The entry ID to assign to the CIF file
        keepIds : bool=False
            If True, keep the residue and chain IDs specified in the Topology
            rather than generating new ones.  Warning: It is up to the caller to
            make sure these are valid IDs that satisfy the requirements of the
            PDBx/mmCIF format.  Otherwise, the output file will be invalid.
        """
        if entry is not None:
            print('data_%s' % entry, file=file)
        else:
            print('data_cell', file=file)
        print("# Created with OpenMM %s, %s" % (Platform.getOpenMMVersion(), str(date.today())), file=file)
        print('#', file=file)
        vectors = topology.getPeriodicBoxVectors()
        if vectors is not None:
            a, b, c, alpha, beta, gamma = computeLengthsAndAngles(vectors)
            RAD_TO_DEG = 180/math.pi
            print('_cell.length_a     %10.4f' % (a*10), file=file)
            print('_cell.length_b     %10.4f' % (b*10), file=file)
            print('_cell.length_c     %10.4f' % (c*10), file=file)
            print('_cell.angle_alpha  %10.4f' % (alpha*RAD_TO_DEG), file=file)
            print('_cell.angle_beta   %10.4f' % (beta*RAD_TO_DEG), file=file)
            print('_cell.angle_gamma  %10.4f' % (gamma*RAD_TO_DEG), file=file)
            print('#', file=file)

        # Identify bonds that should be listed in the file.

        bonds = []
        for atom1, atom2 in topology.bonds():
            if atom1.residue.name not in PDBFile._standardResidues or atom2.residue.name not in PDBFile._standardResidues:
                bonds.append((atom1, atom2))
            elif atom1.name == 'SG' and atom2.name == 'SG' and atom1.residue.name == 'CYS' and atom2.residue.name == 'CYS':
                bonds.append((atom1, atom2))
        if len(bonds) > 0:

            # Write the bond information.

            print('loop_', file=file)
            print('_struct_conn.id', file=file)
            print('_struct_conn.conn_type_id', file=file)
            print('_struct_conn.ptnr1_label_asym_id', file=file)
            print('_struct_conn.ptnr1_label_comp_id', file=file)
            print('_struct_conn.ptnr1_label_seq_id', file=file)
            print('_struct_conn.ptnr1_label_atom_id', file=file)
            print('_struct_conn.ptnr2_label_asym_id', file=file)
            print('_struct_conn.ptnr2_label_comp_id', file=file)
            print('_struct_conn.ptnr2_label_seq_id', file=file)
            print('_struct_conn.ptnr2_label_atom_id', file=file)
            chainIds = {}
            resIds = {}
            if keepIds:
                for chain in topology.chains():
                    chainIds[chain] = chain.id
                for res in topology.residues():
                    resIds[res] = res.id
            else:
                for (chainIndex, chain) in enumerate(topology.chains()):
                    chainIds[chain] = chr(ord('A')+chainIndex%26)
                    for (resIndex, res) in enumerate(chain.residues()):
                        resIds[res] = resIndex+1
            for i, (atom1, atom2) in enumerate(bonds):
                if atom1.residue.name == 'CYS' and atom2.residue.name == 'CYS':
                    bondType = 'disulf'
                else:
                    bondType = 'covale'
                line = "bond%d %s %s %-4s %5s %-4s %s %-4s %5s %-4s"
                print(line % (i+1, bondType, chainIds[atom1.residue.chain], atom1.residue.name, resIds[atom1.residue], atom1.name,
                              chainIds[atom2.residue.chain], atom2.residue.name, resIds[atom2.residue], atom2.name), file=file)
            print('#', file=file)

        # Write the header for the atom coordinates.

        print('loop_', file=file)
        print('_atom_site.group_PDB', file=file)
        print('_atom_site.id', file=file)
        print('_atom_site.type_symbol', file=file)
        print('_atom_site.label_atom_id', file=file)
        print('_atom_site.label_alt_id', file=file)
        print('_atom_site.label_comp_id', file=file)
        print('_atom_site.label_asym_id', file=file)
        print('_atom_site.label_entity_id', file=file)
        print('_atom_site.label_seq_id', file=file)
        print('_atom_site.pdbx_PDB_ins_code', file=file)
        print('_atom_site.Cartn_x', file=file)
        print('_atom_site.Cartn_y', file=file)
        print('_atom_site.Cartn_z', file=file)
        print('_atom_site.occupancy', file=file)
        print('_atom_site.B_iso_or_equiv', file=file)
        print('_atom_site.Cartn_x_esd', file=file)
        print('_atom_site.Cartn_y_esd', file=file)
        print('_atom_site.Cartn_z_esd', file=file)
        print('_atom_site.occupancy_esd', file=file)
        print('_atom_site.B_iso_or_equiv_esd', file=file)
        print('_atom_site.pdbx_formal_charge', file=file)
        print('_atom_site.auth_seq_id', file=file)
        print('_atom_site.auth_comp_id', file=file)
        print('_atom_site.auth_asym_id', file=file)
        print('_atom_site.auth_atom_id', file=file)
        print('_atom_site.pdbx_PDB_model_num', file=file)

    @staticmethod
    def writeModel(topology, positions, file=sys.stdout, modelIndex=1, keepIds=False):
        """Write out a model to a PDBx/mmCIF file.

        Parameters
        ----------
        topology : Topology
            The Topology defining the model to write
        positions : list
            The list of atomic positions to write
        file : file=stdout
            A file to write the model to
        modelIndex : int=1
            The model number of this frame
        keepIds : bool=False
            If True, keep the residue and chain IDs specified in the Topology
            rather than generating new ones.  Warning: It is up to the caller to
            make sure these are valid IDs that satisfy the requirements of the
            PDBx/mmCIF format.  Otherwise, the output file will be invalid.
        """
        if len(list(topology.atoms())) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
        if is_quantity(positions):
            positions = positions.value_in_unit(angstroms)
        import numpy as np
        positions = np.asarray(positions)
        if np.isnan(positions).any():
            raise ValueError('Particle position is NaN.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan')
        if np.isinf(positions).any():
            raise ValueError('Particle position is infinite.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan')
        nonHeterogens = PDBFile._standardResidues[:]
        nonHeterogens.remove('HOH')
        atomIndex = 1
        posIndex = 0
        for (chainIndex, chain) in enumerate(topology.chains()):
            if keepIds:
                chainName = chain.id
            else:
                chainName = chr(ord('A')+chainIndex%26)
            residues = list(chain.residues())
            for (resIndex, res) in enumerate(residues):
                if keepIds:
                    resId = res.id
                    resIC = (res.insertionCode if res.insertionCode.strip() else '.')
                else:
                    resId = resIndex + 1
                    resIC = '.'
                if res.name in nonHeterogens:
                    recordName = "ATOM"
                else:
                    recordName = "HETATM"
                for atom in res.atoms():
                    coords = positions[posIndex]
                    if atom.element is not None:
                        symbol = atom.element.symbol
                    else:
                        symbol = '?'
                    line = "%s  %5d %-3s %-4s . %-4s %s ? %5s %s %10.4f %10.4f %10.4f  0.0  0.0  ?  ?  ?  ?  ?  .  %5s %4s %s %4s %5d"
                    print(line % (recordName, atomIndex, symbol, atom.name, res.name, chainName, resId, resIC, coords[0], coords[1], coords[2],
                                  resId, res.name, chainName, atom.name, modelIndex), file=file)
                    posIndex += 1
                    atomIndex += 1
