"""
tinkerfiles.py: A loader for TINKER files (.xyz, .prm, .key, .seq)

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2025 Stanford University and the Authors.
Authors: Joao Morado
Contributors: Peter Eastman

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, includ  ing without limitation
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

__author__ = "Joao Morado"

import datetime
import io
import math
import os
import re
import shlex
import xml.etree.ElementTree as etree
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from openmm.app import PDBFile
from openmm.app.internal.unitcell import computePeriodicBoxVectors
from openmm.unit import nanometers
from openmm.vec3 import Vec3

from . import element as elem
from . import forcefield as ff
from . import topology as top


@dataclass
class TinkerAtom:
    """
    A data class to represent Tinker atoms.

    Attributes
    ----------
    symbol : str
        The atom symbol (e.g., 'C', 'N', 'O')
    positions : Vec3
        The atom's position in nanometers
    bonds : List[int]
        List of indices of bonded atoms
    atomType : Optional[str]
        The atom type from the force field
    atomClass : Optional[str]
        The atom class from the force field
    nameShort : Optional[str]
        Short name of the atom
    nameLong : Optional[str]
        Long name of the atom
    element : Optional[str]
        The element symbol
    atomicNumber : Optional[int]
        The atomic number
    mass : Optional[float]
        The atom's mass in amu
    valence : Optional[int]
        The atom's valence
    residueName : Optional[str]
        The name of the residue this atom belongs to
    chain : Optional[Any]
        The chain this atom belongs to
    index : int
        The atom's index in the system
    """

    symbol: str
    positions: Vec3
    bonds: List[int]
    index: int
    atomType: Optional[str] = None
    atomClass: Optional[str] = None
    nameShort: Optional[str] = None
    nameLong: Optional[str] = None
    element: Optional[str] = None
    atomicNumber: Optional[int] = None
    mass: Optional[float] = None
    valence: Optional[int] = None
    residueName: Optional[str] = None
    chain: Optional[Any] = None

    def updateFromAtomType(self, atomTypeData: Dict[str, Any]) -> None:
        """Update atom data from atom type information.

        Parameters
        ----------
        atomTypeData : Dict[str, Any]
            Dictionary containing atom type data
        """
        self.atomClass = atomTypeData.atomClass
        self.nameShort = atomTypeData.nameShort
        self.nameLong = atomTypeData.nameLong
        self.element = atomTypeData.element
        self.atomicNumber = atomTypeData.atomicNumber
        self.mass = atomTypeData.mass
        self.valence = atomTypeData.valence


@dataclass
class TinkerAtomType:
    """
    A data class to represent Tinker atom types.
    atomType : str
        The atom type
    atomClass : str
        The atom class
    nameShort : str
        The short name of the atom
    nameLong : str
        The long name of the atom
    element : str
        The element
    atomicNumber : int
        The atomic number
    mass : float
        The mass
    valence : int
        The valence
    """

    atomType: str
    atomClass: str
    nameShort: str
    nameLong: str
    element: str
    atomicNumber: int
    mass: float
    valence: int


class TinkerFiles:
    """TinkerFiles parses Tinker files (.xyz, .prm, .key), constructs a Topology, and (optionally) an OpenMM System from it."""

    @staticmethod
    def _initialize_class():
        """Initialize class variables RECOGNIZED_FORCES and RECOGNIZED_SCALARS upon class creation."""

        def addMultipole(
            lineIndex: int, allLines: List[List[str]], forces: Dict[str, Any]
        ) -> int:
            """
            Parse and store multipole force data from the key file.

            Parameters
            ----------
            lineIndex : int
                The current line index in the key file.
            allLines : list of list of str
                All lines from the key file, split into fields.
            forces : dict
                The forces dictionary to store the parsed data.

            Returns
            -------
            int
                The updated line index after parsing the multipole force data.
            """
            if "multipole" not in forces:
                forces["multipole"] = []
            fields = allLines[lineIndex]
            multipoles = [fields[-1]]
            axisInfo = fields[1:-1]
            lineIndex += 1
            fields = allLines[lineIndex]
            multipoles.append(fields[0])
            multipoles.append(fields[1])
            multipoles.append(fields[2])
            lineIndex += 1
            fields = allLines[lineIndex]
            multipoles.append(fields[0])
            lineIndex += 1
            fields = allLines[lineIndex]
            multipoles.append(fields[0])
            multipoles.append(fields[1])
            lineIndex += 1
            fields = allLines[lineIndex]
            multipoles.append(fields[0])
            multipoles.append(fields[1])
            multipoles.append(fields[2])
            lineIndex += 1
            multipoleInfo = [axisInfo, multipoles]
            forces["multipole"].append(multipoleInfo)
            return lineIndex

        @staticmethod
        def addTorTor(
            lineIndex: int, allLines: List[List[str]], forces: Dict[str, Any]
        ) -> int:
            """
            Parse and store torsion-torsion force data from the key file.

            Parameters
            ----------
            lineIndex : int
                The current line index in the key file.
            allLines : list of list of str
                All lines from the key file, split into fields.
            forces : dict
                The forces dictionary to store the parsed data.

            Returns
            -------
            int
                The updated line index after parsing the torsion-torsion force data.
            """
            if "tortors" not in forces:
                forces["tortors"] = []
            fields = allLines[lineIndex]
            tortorInfo = fields[1:]
            lastGridLine = lineIndex + int(fields[6]) * int(fields[7])
            grid = []
            while lineIndex < lastGridLine:
                lineIndex += 1
                grid.append(allLines[lineIndex])
            forces["tortors"].append([tortorInfo, grid])
            return lineIndex

        RECOGNIZED_FORCES = {
            "bond": 1,
            "angle": 1,
            "strbnd": 1,
            "ureybrad": 1,
            "opbend": 1,
            "torsion": 1,
            "pitors": 1,
            "strtor": 1,
            "angtor": 1,
            "vdw": 1,
            "multipole": addMultipole,
            "tortors": addTorTor,
        }

        RECOGNIZED_SCALARS = {
            "forcefield": "-2.55",
            "bond-cubic": "-2.55",
            "bond-quartic": "3.793125",
            "angle-cubic": "-0.014",
            "angle-quartic": "0.000056",
            "angle-pentic": "-0.0000007",
            "angle-sextic": "0.000000022",
            "opbendtype": "ALLINGER",
            "opbend-cubic": "-0.014",
            "opbend-quartic": "0.000056",
            "opbend-pentic": "-0.0000007",
            "opbend-sextic": "0.000000022",
            "torsionunit": "0.5",
            "vdwtype": "BUFFERED-14-7",
            "radiusrule": "CUBIC-MEAN",
            "radiustype": "R-MIN",
            "radiussize": "DIAMETER",
            "epsilonrule": "HHG",
            "dielectric": "1.0",
            "polarization": "MUTUAL",
            "vdw-13-scale": "0.0",
            "vdw-14-scale": "1.0",
            "vdw-15-scale": "1.0",
            "mpole-12-scale": "0.0",
            "mpole-13-scale": "0.0",
            "mpole-14-scale": "0.4",
            "mpole-15-scale": "0.8",
            "polar-12-scale": "0.0",
            "polar-13-scale": "0.0",
            "polar-14-scale": "1.0",
            "polar-15-scale": "1.0",
            "polar-14-intra": "0.5",
            "direct-11-scale": "0.0",
            "direct-12-scale": "1.0",
            "direct-13-scale": "1.0",
            "direct-14-scale": "1.0",
            "mutual-11-scale": "1.0",
            "mutual-12-scale": "1.0",
            "mutual-13-scale": "1.0",
            "mutual-14-scale": "1.0",
        }

        return RECOGNIZED_FORCES, RECOGNIZED_SCALARS

    # Call the initialization method
    RECOGNIZED_FORCES, RECOGNIZED_SCALARS = _initialize_class()

    def __init__(
        self,
        xyz: str,
        key: str,
        periodicBoxVectors: Optional[Tuple[Vec3, Vec3, Vec3]] = None,
        unitCellDimensions: Optional[Vec3] = None,
        writeXmlFiles: bool = False,
        writePDB: bool = True,
    ):
        """
        Initialize TinkerFiles.

        Parameters
        ----------
        xyz : str
            Path to the Tinker .xyz file.
        key : str
            Path to the Tinker .key file.
        periodicBoxVectors : Optional[Tuple[Vec3, Vec3, Vec3]]
            The periodic box vectors.
        unitCellDimensions : Optional[Vec3]
            The unit cell dimensions.
        writeXmlFiles : bool
            Whether to write XML files.
        writePDB : bool
            Whether to write PDB files.
        """
        # Position and box information
        self.positions = None
        self.boxVectors = None
        self._numpyPositions = None
        self._numpyBoxVectors = None

        # Store data from files
        self.atoms = None
        self.topology = None
        self._atomTypes = dict()
        self._forces = dict()
        self._scalars = dict()

        # Load xyz file
        self.atoms, self.boxVectors, self.positions = TinkerFiles._loadXyzFile(xyz)
        self.positions = self.positions * nanometers

        # Load the .key or .prm file(s)
        key = key if isinstance(key, list) else [key]
        for keyFile in key:
            atomTypes, forces, scalars = self._loadKeyFile(keyFile)
            self._atomTypes.update(atomTypes)
            self._forces.update(forces)
            self._scalars.update(scalars)

        # Update atoms with atom type information
        for atom in self.atoms:
            atomType = atom.atomType
            if atomType not in self._atomTypes:
                raise ValueError(f"Atom type {atomType} not found in atomTypes")
            atom.updateFromAtomType(self._atomTypes[atomType])

        # Create topology
        self.topology = TinkerFiles._createTopology(self.atoms)

        # Set the periodic box vectors as specified in the xyz file
        if self.boxVectors is not None:
            self.topology.setPeriodicBoxVectors(self.boxVectors)

        # If provided, set the periodic box vectors or unit cell dimensions in the topology
        # This overwrites the box information from the xyz file
        if periodicBoxVectors is not None and unitCellDimensions is not None:
            raise ValueError(
                "Specify either periodicBoxVectors or unitCellDimensions, but not both"
            )

        if periodicBoxVectors is not None:
            self.topology.setPeriodicBoxVectors(periodicBoxVectors)
        elif unitCellDimensions is not None:
            self.topology.setUnitCellDimensions(unitCellDimensions)

        # Write XML files if requested
        if writeXmlFiles:
            TinkerXmlWriter.writeXmlFiles(
                self.topology,
                key,
                self.atomTypes,
                self.forces,
                self.scalars,
                self.atoms,
            )

        # Write PDB file if requested
        if writePDB:
            with open(xyz.replace(".xyz", ".pdb"), "w") as f:
                PDBFile.writeFile(self.topology, self.positions, f)

    def createSystem(
        self,
        nonbondedMethod=ff.NoCutoff,
        nonbondedCutoff=1.0 * nanometers,
        constraints=None,
        rigidWater: bool = False,
        removeCMMotion: bool = True,
        hydrogenMass=None,
        polarization: str = "mutual",
        mutualInducedTargetEpsilon: float = 0.00001,
        implicitSolvent: bool = False,
        *args,
        **kwargs,
    ) -> Any:
        """
        Create an OpenMM System from the parsed Tinker files.

        Parameters
        ----------
        nonbondedMethod : object=PME
            The method to use for nonbonded interactions.
            Allowed values are NoCutoff, and PME.
        nonbondedCutoff : distance=1.0*nanometers
            The cutoff distance to use for nonbonded interactions.
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : bool, optional, default=True
            If true, water molecules will be fully rigid regardless of the value passed for the constraints argument.
            Note that AMOEBA waters are flexible.
        removeCMMotion : bool, optional, default=True
            If True, center of mass motion will be removed.
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.
            Any mass added to a hydrogen is subtracted from the heavy atom to keep their total mass the same.
        polarization : str, optional, default="mutual"
            The method to use for calculating induced dipoles.
            Allowed values are "mutual", "direct", or "extrapolated".
        mutualInducedTargetEpsilon : float, optional, default=0.00001
            The target epsilon for mutual induced dipoles.
            Only used if polarization="mutual".
        implicitSolvent : bool, optional, default=False
            If True, solvent will be modeled implicitly.

        Returns
        -------
        openmm.System
            The created OpenMM System.
        """
        from openmm.app.forcefield import ForceField

        # Re-use SystemData to get bond information
        data = ForceField._SystemData(self.topology)
        rigidResidue = [False] * self.topology.getNumResidues()

        # Set atomtypes in SystemData
        for atom in self.topology.atoms():
            data.atomType[atom] = self.atoms[atom.index].atomType

        # Initialize system
        system = ForceField._initializeSystem(
            self.topology,
            self._atomTypes,
            data,
            hydrogenMass,
            nonbondedMethod,
            rigidResidue,
            constraints,
        )

        # Add AMOEBA bond force
        print(self._forces["bond"])
        for type1, type2, length, k in self._forces["bond"]:
            print(type1, type2, length, k)

        return system

    # ------------------------------------------------------------------------------------------ #
    #                                      TOPOLOGY FUNCTIONS                                    #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _createTopology(atoms: List[TinkerAtom]) -> top.Topology:
        """
        Create a topology from the atom dictionary.

        Parameters
        ----------
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        Topology
            The created topology.
        """
        topology = top.Topology()

        # Find molecules using depth-first search through bond network
        molecules = []
        seenAtoms = set()

        for atomId in range(len(atoms)):
            if atomId in seenAtoms or atoms[atomId] is None:
                continue

            # Start a new molecule from this atom
            molecule = []
            stack = [atomId]
            while stack:
                current = stack.pop()
                if current in seenAtoms or atoms[current] is None:
                    continue
                seenAtoms.add(current)
                molecule.append(current)
                stack.extend(
                    atoms[current].bonds
                )  # Add all connected atoms to traverse

            if molecule:
                molecules.append(sorted(molecule))

        # Process each molecule into a chain
        for molecule in molecules:
            processingFunctions = [
                TinkerFiles._processWaterMolecule,
                TinkerFiles._processIonMolecule,
                TinkerFiles._processPeptideChain,
                TinkerFiles._processNucleicAcidChain,
                TinkerFiles._processGenericMolecule,
            ]

            chain = topology.addChain()
            for processingFunction in processingFunctions:
                processingOutput = processingFunction(molecule, atoms)
                if processingOutput:
                    for residueAtoms, residueLabel in processingOutput:
                        residue = topology.addResidue(residueLabel, chain)
                        for atomId in residueAtoms:
                            atom = atoms[atomId]
                            topology.addAtom(
                                atom.nameShort,
                                elem.Element.getByAtomicNumber(atom.atomicNumber),
                                residue,
                                atomId,
                            )
                    break

        # Add the bonds to the topology
        topology_atoms = list(topology.atoms())
        for atom in atoms:
            for bond in atom.bonds:
                if bond > atom.index:
                    topology.addBond(topology_atoms[atom.index], topology_atoms[bond])
        return topology

    @staticmethod
    def _processGenericMolecule(
        molecule: List[int],
        atoms: List[TinkerAtom],
    ) -> List[Union[List[int], str]]:
        """
        Process a generic molecule (unknown).

        Parameters
        ----------
        molecule : List[int]
            List of atom indices in the molecule.
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        List[Union[List[int], str]]
            A list containing [atomIndices, residueName] for the molecule.
        """
        return [[molecule, "UNK"]]

    @staticmethod
    def _processWaterMolecule(
        molecule: List[int],
        atoms: List[TinkerAtom],
    ) -> List[Union[List[int], str]]:
        """
        Process a water molecule.

        Parameters
        ----------
        molecule : List[int]
            List of atom indices in the molecule.
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        List[Union[List[int], str]] or bool
            If a water molecule is identified, returns a list containing [atomIndices, "HOH"].
            Otherwise, returns False.
        """
        if len(molecule) != 3:
            return False

        atomCounts = {atoms[atomId].atomicNumber: 0 for atomId in molecule}
        for atomId in molecule:
            atomCounts[atoms[atomId].atomicNumber] += 1

        if atomCounts.get(1, 0) != 2 or atomCounts.get(8, 0) != 1:
            return False

        return [[molecule, "HOH"]]

    @staticmethod
    def _processPeptideChain(
        molecule: List[int],
        atoms: List[TinkerAtom],
    ) -> List[Union[List[int], str]]:
        """
        Process a peptide chain, identifying all residues and adding them to the topology.

        Parameters
        ----------
        molecule : List[int]
            List of atom indices in the molecule.
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        List[Union[List[int], str]]
            If peptide chain is found, returns [molecule, residueLabel].
            If no peptide chain is found, returns an empty list.
        """
        moleculeAtoms = [None] * len(atoms)
        for atomId in molecule:
            moleculeAtoms[atomId] = atoms[atomId]
        return TinkerFiles._getPeptideResidueAtoms(moleculeAtoms)

    @staticmethod
    def _processNucleicAcidChain(
        molecule: List[int],
        atoms: List[TinkerAtom],
    ) -> List[Union[List[int], str]]:
        """
        Process a nucleic acid chain, identifying all residues.

        Parameters
        ----------
        molecule : List[int]
            List of atom indices in the molecule.
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        List[Union[List[int], str]]
            If nucleic acid chain is found, returns a list of [residueAtoms, residueLabel] pairs.
            If no nucleic acids are found, returns an empty list.
        """
        # Map molecule atoms into full atom array for residue detection
        moleculeAtoms = [None] * len(atoms)
        for atomId in molecule:
            moleculeAtoms[atomId] = atoms[atomId]

        # Identify nucleic acid residues in this molecule
        residues = TinkerFiles._getNucleicAcidResidueAtoms(moleculeAtoms)

        return residues

    @staticmethod
    def _processIonMolecule(
        molecule: List[int],
        atoms: List[TinkerAtom],
    ) -> List[Union[List[int], str]]:
        """
        Process an ion molecule.

        Parameters
        ----------
        molecule : List[int]
            List of atom indices in the molecule.
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        List[Union[List[int], str]] or bool
            If an ion is identified, returns a list containing [atomIndices, ionResidueName].
            Otherwise, returns False.
        """
        if len(molecule) != 1:
            return False

        atom = atoms[molecule[0]]
        resName = TinkerFiles._getIonResidueName(atom.atomicNumber)
        return [[molecule, resName]]

    @staticmethod
    def _getIonResidueName(atomicNumber: int) -> str:
        """
        Get the residue name for an ion based on its atomic number.

        Parameters
        ----------
        atomicNumber : int
            The atomic number of the ion.

        Returns
        -------
        str
            The residue name for the ion (e.g., "Na+", "Cl-").
            Returns "UNK" if the atomic number does not correspond to a known ion.
        """
        ionMap = {
            3: "Li+",
            11: "Na+",
            19: "K+",
            37: "Rb+",
            55: "Cs+",
            4: "Be2+",
            12: "Mg2+",
            20: "Ca2+",
            30: "Zn2+",
            9: "F-",
            17: "Cl-",
            35: "Br-",
            53: "I-",
        }
        return ionMap.get(atomicNumber, "UNK")

    @staticmethod
    def _findBitorsions(
        atoms: List[TinkerAtom],
    ) -> List[Tuple[int, int, int, int, int]]:
        """
        Find all bitorsions (pairs of adjacent torsional angles) in the molecule.

        Notes
        -----
        Follows the logic of the bitors.f and angles.f in Tinker.

        Parameters
        ----------
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        List[Tuple[int, int, int, int, int]]
            List of bitorsion tuples (ia, ib, ic, id, ie) where:
            - ia, ib, ic, id, ie are atom indices
            - ib-ic-id forms an angle
            - ia is connected to ib
            - ie is connected to id
            - ia, ie are not part of the angle
        """
        # First collect all angles (3-atom sequences)
        angles = set()
        for i in range(len(atoms)):
            if atoms[i] is None:
                continue

            ib = i  # Central atom of the angle
            # For each pair of bonds from the central atom
            for j in range(len(atoms[i].bonds) - 1):
                ic = atoms[i].bonds[j]
                for k in range(j + 1, len(atoms[i].bonds)):
                    id = atoms[i].bonds[k]
                    angles.add((ic, ib, id))

        # Now find bitorsions (5-atom sequences) using angles
        bitorsions = set()
        for ib, ic, id in angles:
            # Find external atom connected to ib
            for ia in atoms[ib].bonds:
                if ia == ic or ia == id:
                    continue

                # Find external atom connected to id
                for ie in atoms[id].bonds:
                    if ie == ic or ie == ib or ie == ia:
                        continue

                        bitorsions.add((ia, ib, ic, id, ie))

        return list(bitorsions)

    @staticmethod
    def _addResidueToTopology(
        topology: top.Topology,
        resAtoms: List[int],
        resName: str,
        chain: top.Chain,
        atoms: List[TinkerAtom],
    ) -> None:
        """
        Add a residue to the topology.

        Parameters
        ----------
        topology : openmm.app.topology.Topology
            The topology to add the residue to.
        resAtoms : List[int]
            List of atom indices in the residue.
        resName : str
            The name of the residue.
        chain : Chain
            The chain to add the residue to.
        atoms : List[TinkerAtom]
            List of atoms in the system.
        """
        residue = topology.addResidue(resName, chain)
        for atomIndex in resAtoms:
            atom = atoms[atomIndex]
            topology.addAtom(
                atom.nameShort or atom.element or str(atom.atomicNumber),
                elem.Element.getByAtomicNumber(atom.atomicNumber),
                residue,
            )

    # ------------------------------------------------------------------------------------------ #
    #                                     PEPTIDE PROCESSING                                     #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _getPeptideResidueAtoms(
        atoms: List[TinkerAtom], atomIndex: Optional[int] = None
    ) -> List[Tuple[List[int], str]]:
        """
        Get all amino acid residues in the molecule.
        If atomIndex is provided, only return residues that contain that atom.

        Parameters
        ----------
        atoms : List[TinkerAtom]
            List of atoms in the system
        atomIndex : Optional[int]
            Optional index of an atom to filter residues by

        Returns
        -------
        List[Tuple[List[int], str]]
            List of (residueAtoms, residueLabel) tuples for each amino acid residue found
        """
        residues = []
        seenAtoms = set()

        # Get all bitorsions in the molecule
        bitorsions = TinkerFiles._findBitorsions(atoms)
        for ia, ib, ic, id, ie in bitorsions:
            if any(atom in seenAtoms for atom in [ib, ic, id]):
                continue

            if atomIndex is not None and atomIndex not in [ib, ic, id]:
                continue

            # Check for N/C-terminal cap groups (ACE, NME)
            capAtoms = TinkerFiles._checkCapGroup(ia, ib, ic, id, ie, atoms)
            if capAtoms is not None:
                # Collect all atoms in the cap group including hydrogens
                allAtoms = set(capAtoms)
                for atom in capAtoms:
                    allAtoms.update(atoms[atom].bonds)
                allAtoms = sorted(list(allAtoms))
                seenAtoms.update(allAtoms)
                # Determine if it's an acetyl (ACE) or N-methyl (NME) cap
                residues.append(
                    (allAtoms, "ACE" if atoms[ia].atomicNumber == 6 else "NME")
                )
                continue

            # Backbone pattern
            if TinkerFiles._checkBackbonePattern(ia, ib, ic, id, ie, atoms):
                residueLabel, sideChain = TinkerFiles._identifyAminoAcid(
                    ic, ib, id, atoms
                )
                if residueLabel:
                    residueAtoms = TinkerFiles._collectResidueAtoms(
                        ic, ib, id, sideChain, atoms
                    )
                    seenAtoms.update(residueAtoms)
                    residues.append((residueAtoms, residueLabel))
                continue

            # Terminal pattern
            if TinkerFiles._checkTerminalPattern(ia, ib, ic, id, ie, atoms):
                residueLabel, sideChain = TinkerFiles._identifyAminoAcid(
                    ic, ib, id, atoms
                )
                if residueLabel:
                    residueAtoms = TinkerFiles._collectResidueAtoms(
                        ic, ib, id, sideChain, atoms
                    )
                    seenAtoms.update(residueAtoms)
                    residues.append((residueAtoms, residueLabel))
                continue

        return sorted(residues, key=lambda x: x[0][0])

    @staticmethod
    def _checkBackbonePattern(
        ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]
    ) -> bool:
        """
        Check if a 5-atom sequence matches a peptide backbone pattern.
        """
        atomA = atoms[ia]
        atomB = atoms[ib]
        atomC = atoms[ic]
        atomD = atoms[id]
        atomE = atoms[ie]

        # Identify peptide backbone pattern
        if (
            len(atomA.bonds) == 3
            and len(atomB.bonds) == 3
            and len(atomC.bonds) == 4
            and len(atomD.bonds) == 3
            and len(atomE.bonds) == 3
        ):
            if (
                atomA.atomicNumber == 6
                and atomB.atomicNumber == 7
                and atomC.atomicNumber == 6
                and atomD.atomicNumber == 6
                and atomE.atomicNumber == 7
            ) or (
                atomA.atomicNumber == 7
                and atomB.atomicNumber == 6
                and atomC.atomicNumber == 6
                and atomD.atomicNumber == 7
                and atomE.atomicNumber == 6
            ):
                # Check for carbonyl oxygen
                if atomA.atomicNumber == 6:
                    for j in range(len(atomA.bonds)):
                        ij = atomA.bonds[j]
                        if len(atoms[ij].bonds) == 1 and atoms[ij].atomicNumber == 8:
                            return True

                if atomB.atomicNumber == 6:
                    for j in range(len(atomB.bonds)):
                        ij = atomB.bonds[j]
                        if len(atoms[ij].bonds) == 1 and atoms[ij].atomicNumber == 8:
                            return True

                if atomD.atomicNumber == 6:
                    for j in range(len(atomD.bonds)):
                        ij = atomD.bonds[j]
                        if len(atoms[ij].bonds) == 1 and atoms[ij].atomicNumber == 8:
                            return True

                if atomE.atomicNumber == 6:
                    for j in range(len(atomE.bonds)):
                        ij = atomE.bonds[j]
                        if len(atoms[ij].bonds) == 1 and atoms[ij].atomicNumber == 8:
                            return True

        return False

    @staticmethod
    def _checkTerminalPattern(
        ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]
    ) -> bool:
        """
        Check if a 5-atom sequence matches a terminal residue pattern.
        """
        atomA = atoms[ia]
        atomB = atoms[ib]
        atomC = atoms[ic]
        atomD = atoms[id]
        atomE = atoms[ie]

        # H-N-C-C-N
        if (
            len(atomA.bonds) == 1
            and len(atomB.bonds) == 4
            and len(atomC.bonds) == 4
            and len(atomD.bonds) == 3
            and len(atomE.bonds) == 3
            and atomA.atomicNumber == 1
            and atomB.atomicNumber == 7
            and atomC.atomicNumber == 6
            and atomD.atomicNumber == 6
            and atomE.atomicNumber == 7
        ):
            return True

        # N-C-C-N-H
        if (
            len(atomA.bonds) == 3
            and len(atomB.bonds) == 3
            and len(atomC.bonds) == 4
            and len(atomD.bonds) == 4
            and len(atomE.bonds) == 1
            and atomA.atomicNumber == 7
            and atomB.atomicNumber == 6
            and atomC.atomicNumber == 6
            and atomD.atomicNumber == 7
            and atomE.atomicNumber == 1
        ):
            return True

        # C-N-C-C-N
        if (
            len(atomA.bonds) == 3
            and len(atomB.bonds) == 3
            and len(atomC.bonds) == 4
            and len(atomD.bonds) == 3
            and len(atomE.bonds) == 1
            and atomA.atomicNumber == 6
            and atomB.atomicNumber == 7
            and atomC.atomicNumber == 6
            and atomD.atomicNumber == 6
            and atomE.atomicNumber == 8
        ):
            return True

        # O-C-C-N-C
        if (
            len(atomA.bonds) == 1
            and len(atomB.bonds) == 3
            and len(atomC.bonds) == 4
            and len(atomD.bonds) == 3
            and len(atomE.bonds) == 3
            and atomA.atomicNumber == 8
            and atomB.atomicNumber == 6
            and atomC.atomicNumber == 6
            and atomD.atomicNumber == 7
            and atomE.atomicNumber == 6
        ):
            return True

        return False

    @staticmethod
    def _checkCapGroup(
        ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]
    ) -> Optional[List[int]]:
        """
        Check if a 5-atom sequence matches a cap group pattern (ACE or NME).

        Parameters
        ----------
        ia, ib, ic, id, ie : int
            Indices of the 5 atoms in the sequence
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        Optional[List[int]]
            List of atom indices in the cap group if found, None otherwise
        """
        atomA = atoms[ia]
        atomB = atoms[ib]
        atomC = atoms[ic]
        atomD = atoms[id]
        atomE = atoms[ie]

        # ACE cap group pattern
        if (
            len(atomA.bonds) == 4
            and len(atomB.bonds) == 3
            and len(atomC.bonds) == 3
            and len(atomD.bonds) == 4
            and len(atomE.bonds) == 3
            and atomA.atomicNumber == 6
            and atomB.atomicNumber == 6
            and atomC.atomicNumber == 7
            and atomD.atomicNumber == 6
            and atomE.atomicNumber == 6
        ):
            # Count hydrogens on first carbon
            nHyd = sum(1 for bond in atomA.bonds if atoms[bond].atomicNumber == 1)
            if nHyd == 3:
                capAtoms = [ia]
                for bond in atomA.bonds:
                    if atoms[bond].atomicNumber == 1:
                        capAtoms.append(bond)
                capAtoms.append(ib)
                for bond in atomB.bonds:
                    if atoms[bond].atomicNumber == 8:
                        capAtoms.append(bond)
                return capAtoms

        # NME cap group pattern
        if (
            len(atomA.bonds) == 3
            and len(atomB.bonds) == 4
            and len(atomC.bonds) == 3
            and len(atomD.bonds) == 4
            and len(atomE.bonds) == 3
            and atomA.atomicNumber == 7
            and atomB.atomicNumber == 6
            and atomC.atomicNumber == 6
            and atomD.atomicNumber == 6
            and atomE.atomicNumber == 7
        ):
            # Count hydrogens on last carbon
            nHyd = sum(1 for bond in atomE.bonds if atoms[bond].atomicNumber == 1)
            if nHyd == 3:
                capAtoms = [id]
                for bond in atomD.bonds:
                    if atoms[bond].atomicNumber == 1:
                        capAtoms.append(bond)
                capAtoms.append(ie)
                for bond in atomE.bonds:
                    if atoms[bond].atomicNumber == 1:
                        capAtoms.append(bond)
                return capAtoms

        return None

    @staticmethod
    def _identifyAminoAcid(
        ic: int, ib: int, id: int, atoms: List[TinkerAtom]
    ) -> Tuple[str, Dict[str, int]]:
        """
        Identify amino acid type based on side chain atoms.

        Notes
        -----
        Follows logic from findpro.f in Tinker.

        Parameters
        ----------
        ic : int
            Alpha carbon index
        ib : int
            Nitrogen index
        id : int
            Carbonyl carbon index
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        Tuple[str, Dict[str, int]]
            Amino acid label and dictionary of side chain atom indices
        """
        # Zero out residue name and possible side chain atoms
        # Format: "i" + element + position + number (e.g., "icb1" = carbon beta 1)
        # Position codes: b=beta, g=gamma, d=delta, e=epsilon, z=zeta
        sideChain = {
            # Carbon atoms along side chain
            "icb1": 0,  # First beta carbon (Cβ1)
            "icb2": 0,  # Second beta carbon (Cβ2)
            "icg1": 0,  # First gamma carbon (Cγ1)
            "icg2": 0,  # Second gamma carbon (Cγ2)
            "icd1": 0,  # First delta carbon (Cδ1)
            "icd2": 0,  # Second delta carbon (Cδ2)
            "ice1": 0,  # First epsilon carbon (Cε1)
            "ice2": 0,  # Second epsilon carbon (Cε2)
            "icz1": 0,  # First zeta carbon (Cζ1)
            "icz2": 0,  # Second zeta carbon (Cζ2)
            # Oxygen atoms
            "iog": 0,  # Gamma oxygen (Oγ)
            "ioh": 0,  # Hydroxyl oxygen (Oη)
            "iod1": 0,  # First delta oxygen (Oδ1)
            "iod2": 0,  # Second delta oxygen (Oδ2)
            "ioe1": 0,  # First epsilon oxygen (Oε1)
            "ioe2": 0,  # Second epsilon oxygen (Oε2)
            # Nitrogen atoms
            "ind": 0,  # Delta nitrogen (Nδ)
            "ine": 0,  # Epsilon nitrogen (Nε)
            "inz": 0,  # Zeta nitrogen (Nζ)
            "inh1": 0,  # First hydrogen on zeta nitrogen (NζH1)
            "inh2": 0,  # Second hydrogen on zeta nitrogen (NζH2)
            # Sulfur atoms
            "isg": 0,  # Gamma sulfur (Sγ)
            "isd": 0,  # Delta sulfur (Sδ)
            # Ring carbon (for aromatic residues)
            "ich": 0,  # Ring carbon (Cη)
        }

        label = "   "
        atomC = atoms[ic]  # Alpha carbon

        # Inspect the beta position of amino acid residue
        nHa = 0
        for bond in atomC.bonds:
            if atoms[bond].atomicNumber == 1:
                nHa += 1
            elif bond != ib and bond != id and atoms[bond].atomicNumber == 6:
                if sideChain["icb1"] != 0:
                    sideChain["icb2"] = bond
                else:
                    sideChain["icb1"] = bond
        if nHa == 2:
            label = "GLY"
            return label, sideChain
        if sideChain["icb2"] != 0:
            label = "AIB"
            return label, sideChain

        # Inspect the gamma position of amino acid residue
        if sideChain["icb1"] != 0:
            nHb = 0
            for bond in atoms[sideChain["icb1"]].bonds:
                if atoms[bond].atomicNumber == 1:
                    nHb += 1
                elif bond != ic and atoms[bond].atomicNumber == 6:
                    if sideChain["icg1"] != 0:
                        sideChain["icg2"] = bond
                    else:
                        sideChain["icg1"] = bond
                elif atoms[bond].atomicNumber == 8:
                    sideChain["iog"] = bond
                elif atoms[bond].atomicNumber == 16:
                    sideChain["isg"] = bond
            if nHb == 3:
                label = "ALA"
                return label, sideChain
            if sideChain["iog"] != 0 and sideChain["icg1"] == 0:
                label = "SER"
                return label, sideChain
            if sideChain["iog"] != 0 and sideChain["icg1"] != 0:
                label = "THR"
                return label, sideChain
            if sideChain["isg"] != 0:
                if len(atoms[sideChain["isg"]].bonds) == 1:
                    label = "CYD"
                else:
                    for bond in atoms[sideChain["isg"]].bonds:
                        if atoms[bond].atomicNumber == 1:
                            label = "CYS"

                        elif atoms[bond].atomicNumber == 16:
                            sideChain["isg"] = 0
                            label = "CYX"
                    return label, sideChain
            if sideChain["icg2"] != 0:
                nhg = 0
                for bond in atoms[sideChain["icg1"]].bonds:
                    if atoms[bond].atomicNumber == 1:
                        nhg += 1
                    elif bond != sideChain["icb1"] and atoms[bond].atomicNumber == 6:
                        sideChain["icd1"] = bond
                for bond in atoms[sideChain["icg2"]].bonds:
                    if atoms[bond].atomicNumber == 1:
                        nhg += 1
                    elif bond != sideChain["icb1"] and atoms[bond].atomicNumber == 6:
                        sideChain["icd1"] = bond
                if nhg == 5:
                    label = "ILE"
                    return label, sideChain
                if nhg == 6:
                    label = "VAL"
                    return label, sideChain

        # Inspect the delta position of amino acid residue
        if sideChain["icg1"] != 0:
            for bond in atoms[sideChain["icg1"]].bonds:
                if bond != sideChain["icb1"]:
                    if atoms[bond].atomicNumber == 6:
                        if sideChain["icd1"] != 0:
                            sideChain["icd2"] = bond
                        else:
                            sideChain["icd1"] = bond
                    elif atoms[bond].atomicNumber == 7:
                        sideChain["ind"] = bond
                    elif atoms[bond].atomicNumber == 8:
                        if sideChain["iod1"] != 0:
                            sideChain["iod2"] = bond
                        else:
                            sideChain["iod1"] = bond
                    elif atoms[bond].atomicNumber == 16:
                        sideChain["isd"] = bond
                        label = "MET"
                        for subbond in atoms[bond].bonds:
                            if (
                                subbond != sideChain["icg1"]
                                and atoms[subbond].atomicNumber == 6
                            ):
                                sideChain["ice1"] = subbond
                        return label, sideChain

            if sideChain["iod1"] != 0 and sideChain["ind"] != 0:
                label = "ASN"
                return label, sideChain
            if sideChain["iod2"] != 0:
                if (
                    len(atoms[sideChain["iod1"]].bonds) == 1
                    and len(atoms[sideChain["iod2"]].bonds) == 1
                ):
                    label = "ASP"
                else:
                    label = "ASH"
                return label, sideChain
            if sideChain["icd2"] != 0:
                nhd = 0
                for bond in atoms[sideChain["icd1"]].bonds:
                    if atoms[bond].atomicNumber == 1:
                        nhd += 1
                for bond in atoms[sideChain["icd2"]].bonds:
                    if atoms[bond].atomicNumber == 1:
                        nhd += 1
                if nhd == 6:
                    label = "LEU"
                    return label, sideChain
            if sideChain["icd1"] != 0:
                for bond in atoms[sideChain["icd1"]].bonds:
                    if bond == ib or bond == id:
                        label = "PRO"
                        return label, sideChain

        # Inspect the epsilon position of amino acid residue
        if sideChain["icd1"] != 0:
            for bond in atoms[sideChain["icd1"]].bonds:
                if bond != sideChain["icg1"]:
                    if atoms[bond].atomicNumber == 6:
                        if sideChain["ice1"] != 0:
                            sideChain["ice2"] = bond
                        else:
                            sideChain["ice1"] = bond
                    elif atoms[bond].atomicNumber == 7:
                        sideChain["ine"] = bond
                    elif atoms[bond].atomicNumber == 8:
                        if sideChain["ioe1"] != 0:
                            sideChain["ioe2"] = bond
                        else:
                            sideChain["ioe1"] = bond

            if sideChain["ine"] != 0:
                if len(atoms[sideChain["ine"]].bonds) == 4:
                    label = "ORN"
                    return label, sideChain
            if sideChain["ioe1"] != 0 and sideChain["ine"] != 0:
                label = "GLN"
                return label, sideChain
            if sideChain["ioe2"] != 0:
                if (
                    len(atoms[sideChain["ioe1"]].bonds) == 1
                    and len(atoms[sideChain["ioe2"]].bonds) == 1
                ):
                    label = "GLU"
                else:
                    label = "GLH"
                return label, sideChain

        if sideChain["icd2"] != 0:
            for bond in atoms[sideChain["icd2"]].bonds:
                if bond != sideChain["icg1"]:
                    if atoms[bond].atomicNumber == 6:
                        if sideChain["ice1"] != 0:
                            sideChain["ice2"] = bond
                        else:
                            sideChain["ice1"] = bond
                    elif atoms[bond].atomicNumber == 7:
                        sideChain["ine"] = bond

        if sideChain["ind"] != 0:
            for bond in atoms[sideChain["ind"]].bonds:
                if bond != sideChain["icg1"] and atoms[bond].atomicNumber == 6:
                    if sideChain["ice1"] != 0:
                        sideChain["ice2"] = bond
                    else:
                        sideChain["ice1"] = bond

        if min(sideChain["ind"], sideChain["ine"]) != 0:
            for bond in atoms[sideChain["ine"]].bonds:
                if bond == sideChain["ice1"]:
                    label = "HIS"
                    if len(atoms[sideChain["ind"]].bonds) == 2:
                        label = "HIE"
                        return label, sideChain
                    if len(atoms[sideChain["ine"]].bonds) == 2:
                        label = "HID"
                        return label, sideChain
                    return label, sideChain

        if sideChain["ine"] != 0:
            for bond in atoms[sideChain["ine"]].bonds:
                if bond == sideChain["ice1"]:
                    label = "TRP"
                    # Store ring atoms for TRP
                    for ebond in atoms[sideChain["ice1"]].bonds:
                        if (
                            ebond != sideChain["icd1"]
                            and atoms[ebond].atomicNumber == 6
                        ):
                            if sideChain["icz1"] != 0:
                                sideChain["icz2"] = ebond
                            else:
                                sideChain["icz1"] = ebond
                    # Look for CH atom in the ring
                    for zbond in atoms[sideChain["icz1"]].bonds:
                        if (
                            zbond != sideChain["ice1"]
                            and atoms[zbond].atomicNumber == 6
                        ):
                            sideChain["ich"] = zbond
                    return label, sideChain
                elif bond == sideChain["ice2"]:
                    label = "TRP"
                    # Store ring atoms for TRP
                    for ebond in atoms[sideChain["ice2"]].bonds:
                        if (
                            ebond != sideChain["icd1"]
                            and atoms[ebond].atomicNumber == 6
                        ):
                            if sideChain["icz1"] != 0:
                                sideChain["icz2"] = ebond
                            else:
                                sideChain["icz1"] = ebond
                    # Look for CH atom in the ring
                    for zbond in atoms[sideChain["icz1"]].bonds:
                        if (
                            zbond != sideChain["ice2"]
                            and atoms[zbond].atomicNumber == 6
                        ):
                            sideChain["ich"] = zbond
                    return label, sideChain

        # Inspect the zeta position of amino acid residue
        if sideChain["ice1"] != 0:
            for bond in atoms[sideChain["ice1"]].bonds:
                if bond != sideChain["icd1"]:
                    if atoms[bond].atomicNumber == 6:
                        if sideChain["icz1"] != 0:
                            sideChain["icz2"] = bond
                        else:
                            sideChain["icz1"] = bond
                    elif atoms[bond].atomicNumber == 7:
                        sideChain["inz"] = bond

            if sideChain["inz"] != 0 and len(atoms[sideChain["ice1"]].bonds) == 4:
                if len(atoms[sideChain["inz"]].bonds) == 3:
                    label = "LYD"
                    return label, sideChain
                elif len(atoms[sideChain["inz"]].bonds) == 4:
                    label = "LYS"
                    return label, sideChain

        if sideChain["ice2"] != 0:
            for bond in atoms[sideChain["ice2"]].bonds:
                if bond != sideChain["icd2"]:
                    if atoms[bond].atomicNumber == 6:
                        if sideChain["icz1"] != 0:
                            sideChain["icz2"] = bond
                        else:
                            sideChain["icz1"] = bond

        if sideChain["icz1"] == sideChain["icz2"]:
            sideChain["icz2"] = 0
            label = "PHE"
            for bond in atoms[sideChain["icz1"]].bonds:
                if atoms[bond].atomicNumber == 8:
                    sideChain["ioh"] = bond
                    if len(atoms[bond].bonds) == 1:
                        label = "TYD"
                    elif len(atoms[bond].bonds) == 2:
                        label = "TYR"
                    return label, sideChain

        if sideChain["ine"] != 0:
            for bond in atoms[sideChain["ine"]].bonds:
                if bond != sideChain["icd1"] and atoms[bond].atomicNumber == 6:
                    if sideChain["icz1"] != 0:
                        sideChain["icz2"] = bond
                    else:
                        sideChain["icz1"] = bond

            if sideChain["icz1"] != 0:
                label = "ARG"
                for bond in atoms[sideChain["icz1"]].bonds:
                    if atoms[bond].atomicNumber != 7:
                        label = "   "
                    elif bond != sideChain["ine"]:
                        if sideChain["inh1"] != 0:
                            sideChain["inh2"] = bond
                        else:
                            sideChain["inh1"] = bond
                return label, sideChain

        return label, sideChain

    @staticmethod
    def _collectResidueAtoms(
        ic: int, ib: int, id: int, sideChain: Dict[str, int], atoms: List[TinkerAtom]
    ) -> List[int]:
        """
        Collect all atoms in the residue.

        Parameters
        ----------
        ic, ib, id : int
            Alpha carbon, nitrogen, and carbonyl carbon indices
        sideChain : Dict[str, int]
            Dictionary of side chain atom indices
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        List[int]
            List of all atom indices in the residue
        """
        residueAtoms = {ic}  # Alpha carbon

        # Add backbone atoms (N-H, C=O)
        residueAtoms.add(ib)
        for bond in atoms[ib].bonds:
            if atoms[bond].atomicNumber == 1:
                residueAtoms.add(bond)

        residueAtoms.add(id)  # Add carbonyl carbon
        for bond in atoms[id].bonds:
            if atoms[bond].atomicNumber == 8:  # Add carbonyl oxygen
                residueAtoms.add(bond)
            elif atoms[bond].atomicNumber == 1:  # Add any hydrogens
                residueAtoms.add(bond)

        # Add side chain atoms
        for key, index in sideChain.items():
            if index != 0:
                residueAtoms.add(index)
                residueAtoms.update(atoms[index].bonds)

        # Add all hydrogen atoms
        for atom in list(residueAtoms):
            for bond in atoms[atom].bonds:
                if atoms[bond].atomicNumber == 1:
                    residueAtoms.add(bond)

        return sorted(list(residueAtoms))

    # ------------------------------------------------------------------------------------------ #
    #                                   NUCLEIC ACID PROCESSING                                  #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _getNucleicAcidResidueAtoms(
        atoms: List[TinkerAtom],
    ) -> List[Tuple[List[int], str]]:
        """
        Get all nucleic acid residues in the molecule.

        Parameters
        ----------
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        List[Tuple[List[int], str]]
            List of (residueAtoms, residueLabel) tuples for each nucleic acid residue found
        """
        residues = []
        seenAtoms = set()

        # Get all bitorsions in the molecule
        bitorsions = TinkerFiles._findBitorsions(atoms)

        for ia, ib, ic, id, ie in bitorsions:
            if any(atom in seenAtoms for atom in [ib, ic, id]):
                continue

            # Check if the sequence matches a nucleotide pattern
            if not TinkerFiles._checkNucleotidePattern(ia, ib, ic, id, ie, atoms):
                continue

            # Identify sugar ring and phosphate atoms
            label, nucleotideAtoms = TinkerFiles._identifyNucleotide(
                ia, ib, ic, id, ie, atoms
            )

            # Extract all atoms in this nucleotide
            residueAtoms = TinkerFiles._collectNucleotideAtoms(nucleotideAtoms, atoms)
            seenAtoms.update(residueAtoms)
            residues.append((sorted(residueAtoms), label))

        return sorted(residues, key=lambda x: x[0][0])

    @staticmethod
    def _checkNucleotidePattern(
        ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]
    ) -> bool:
        """
        Check if a 5-atom sequence matches a nucleotide pattern (O5' -> C5' -> C4' -> C3' -> O3').

        Parameters
        ----------
        ia, ib, ic, id, ie : int
            Indices of the 5 atoms in the sequence
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        bool
            True if the sequence matches a nucleotide pattern, False otherwise
        """
        atomA = atoms[ia]
        atomB = atoms[ib]
        atomC = atoms[ic]
        atomD = atoms[id]
        atomE = atoms[ie]

        # Atomic numbers
        if not (
            atomA.atomicNumber == 8
            and atomB.atomicNumber == 6
            and atomC.atomicNumber == 6
            and atomD.atomicNumber == 6
            and atomE.atomicNumber == 8
        ):
            return False

        # Connectivity
        if not (
            len(atomA.bonds) == 2
            and len(atomB.bonds) == 4
            and len(atomC.bonds) == 4
            and len(atomD.bonds) == 4
            and len(atomE.bonds) == 2
        ):
            return False

        # Count hydrogens and phosphorus atoms connected to O5' and O3'
        nHyd = sum(
            1 for bond in atomA.bonds + atomE.bonds if atoms[bond].atomicNumber == 1
        )
        nPhos = sum(
            1 for bond in atomA.bonds + atomE.bonds if atoms[bond].atomicNumber == 15
        )

        return (nPhos == 1 and nHyd == 1) or (nPhos == 2 and nHyd == 0)

    @staticmethod
    def _identifyNucleotide(
        ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]
    ) -> Tuple[str, Dict[str, int]]:
        """
        Identify a nucleotide and its atoms.

        Parameters
        ----------
        ia, ib, ic, id, ie : int
            Indices of the 5 atoms in the sequence
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        Tuple[str, Dict[str, int]]
            Label of the nucleotide and dictionary of nucleotide atom indices
        """
        # Initialize with empty values
        nucleotideAtoms = {
            "io5s": 0,  # 5' sugar oxygen
            "ic5s": 0,  # 5' carbon
            "ic4s": 0,  # 4' carbon
            "ic3s": 0,  # 3' carbon
            "ic2s": 0,  # 2' carbon
            "io4s": 0,  # Ring oxygen
            "ic1s": 0,  # 1' carbon (connects to base)
            "io3s": 0,  # 3' sugar oxygen
            "io2s": 0,  # 2' hydroxyl oxygen (absent in deoxyribose)
            "ipo": 0,  # Phosphate phosphorus
            "iop1": 0,  # First phosphate oxygen
            "iop2": 0,  # Second phosphate oxygen
            "iop3": 0,  # Third phosphate oxygen
            "in1": 0,  # N1 nitrogen (connects to sugar)
            "in3": 0,  # N3 nitrogen
            "in4": 0,  # N4 nitrogen
            "in5": 0,  # N5 nitrogen
            "in6": 0,  # N6 nitrogen
            "in7": 0,  # N7 nitrogen
            "in8": 0,  # N8 nitrogen
            "in9": 0,  # N9 nitrogen
            "ic2": 0,  # C2 carbon in base
            "ic3": 0,  # C3 carbon in base
            "ic4": 0,  # C4 carbon in base
            "ic5": 0,  # C5 carbon in base
            "ic6": 0,  # C6 carbon in base
            "ic7": 0,  # C7 carbon in base
            "ic9": 0,  # C9 carbon in base
            "icm": 0,  # Methyl carbon (in thymine)
            "io2": 0,  # O2 oxygen (for pyrimidines)
            "io4": 0,  # O4 oxygen
            "io9": 0,  # O9 oxygen
        }

        label = "   "

        # Locate sugar ring and assign corresponding atom names
        # First determine orientation by hydrogen count on C5'
        nHyd = 0
        for bond in atoms[ib].bonds:
            if atoms[bond].atomicNumber == 1:
                nHyd += 1

        if nHyd == 2:
            # 5' to 3' orientation
            nucleotideAtoms["io5s"] = ia
            nucleotideAtoms["ic5s"] = ib
            nucleotideAtoms["ic4s"] = ic
            nucleotideAtoms["ic3s"] = id
            nucleotideAtoms["io3s"] = ie
        else:
            # 3' to 5' orientation
            nucleotideAtoms["io5s"] = ie
            nucleotideAtoms["ic5s"] = id
            nucleotideAtoms["ic4s"] = ic
            nucleotideAtoms["ic3s"] = ib
            nucleotideAtoms["io3s"] = ia

        # Find ring oxygen (O4')
        for bond in atoms[nucleotideAtoms["ic4s"]].bonds:
            if (
                atoms[bond].atomicNumber == 8
                and bond != nucleotideAtoms["ic5s"]
                and bond != nucleotideAtoms["ic3s"]
            ):
                nucleotideAtoms["io4s"] = bond

        # Find C2' carbon
        for bond in atoms[nucleotideAtoms["ic3s"]].bonds:
            if bond != nucleotideAtoms["ic4s"] and atoms[bond].atomicNumber == 6:
                nucleotideAtoms["ic2s"] = bond

        # Check for 2' hydroxyl oxygen (absent in DNA)
        deoxy = True
        if nucleotideAtoms["ic2s"] != 0:
            for bond in atoms[nucleotideAtoms["ic2s"]].bonds:
                if atoms[bond].atomicNumber == 8:
                    nucleotideAtoms["io2s"] = bond
                    deoxy = False

        # Find C1' carbon connected to ring oxygen
        if nucleotideAtoms["io4s"] != 0:
            for bond in atoms[nucleotideAtoms["io4s"]].bonds:
                if bond != nucleotideAtoms["ic4s"] and atoms[bond].atomicNumber == 6:
                    nucleotideAtoms["ic1s"] = bond

        # Find phosphate group attached at 5' sugar oxygen
        for bond in atoms[nucleotideAtoms["io5s"]].bonds:
            if atoms[bond].atomicNumber == 15:
                nucleotideAtoms["ipo"] = bond

        # Find phosphate oxygens
        if nucleotideAtoms["ipo"] != 0:
            for bond in atoms[nucleotideAtoms["ipo"]].bonds:
                if (
                    bond != nucleotideAtoms["io5s"]
                    and atoms[bond].atomicNumber == 8
                    and len(atoms[bond].bonds) == 1
                ):
                    if nucleotideAtoms["iop2"] != 0:
                        nucleotideAtoms["iop3"] = bond
                    elif nucleotideAtoms["iop1"] != 0:
                        nucleotideAtoms["iop2"] = bond
                    else:
                        nucleotideAtoms["iop1"] = bond

        # Find N1 nitrogen connecting sugar to base
        if nucleotideAtoms["ic1s"] != 0:
            for bond in atoms[nucleotideAtoms["ic1s"]].bonds:
                if atoms[bond].atomicNumber == 7:
                    nucleotideAtoms["in1"] = bond

        # Identify if purine or pyrimidine and set base atoms
        baseType = "   "
        if nucleotideAtoms["in1"] != 0:
            for bond in atoms[nucleotideAtoms["in1"]].bonds:
                if bond != nucleotideAtoms["ic1s"]:
                    for subbond in atoms[bond].bonds:
                        if atoms[subbond].atomicNumber == 1:
                            nucleotideAtoms["ic2"] = bond
                        if atoms[subbond].atomicNumber == 8:
                            nucleotideAtoms["io2"] = subbond
                            nucleotideAtoms["ic6"] = bond
                            baseType = "PYR"

                    if nucleotideAtoms["ic6"] == 0:
                        baseType = "PUR"
                        nucleotideAtoms["ic5"] = bond

        # For purines: adenine (A) or guanine (G)
        if baseType == "PUR":
            # Find N3 nitrogen attached to C2
            if nucleotideAtoms["ic2"] != 0:
                for bond in atoms[nucleotideAtoms["ic2"]].bonds:
                    if bond != nucleotideAtoms["in1"] and atoms[bond].atomicNumber == 7:
                        nucleotideAtoms["in3"] = bond

            # Find C4 carbon attached to N3
            if nucleotideAtoms["in3"] != 0:
                for bond in atoms[nucleotideAtoms["in3"]].bonds:
                    if bond != nucleotideAtoms["ic2"] and atoms[bond].atomicNumber == 6:
                        nucleotideAtoms["ic4"] = bond

            # Find N6 nitrogen in the ring
            if nucleotideAtoms["ic5"] != 0:
                for bond in atoms[nucleotideAtoms["ic5"]].bonds:
                    if bond != nucleotideAtoms["in1"] and atoms[bond].atomicNumber == 7:
                        nucleotideAtoms["in6"] = bond

            # Find C7 carbon connected to N6
            if nucleotideAtoms["in6"] != 0:
                for bond in atoms[nucleotideAtoms["in6"]].bonds:
                    if bond != nucleotideAtoms["ic5"] and atoms[bond].atomicNumber == 6:
                        nucleotideAtoms["ic7"] = bond

            # Find N7/N8 based on hydrogen count
            if nucleotideAtoms["ic7"] != 0:
                for bond in atoms[nucleotideAtoms["ic7"]].bonds:
                    if bond != nucleotideAtoms["in6"] and atoms[bond].atomicNumber == 7:
                        # Count hydrogens to determine if it's N7 or N8
                        nHyd = 0
                        for subbond in atoms[bond].bonds:
                            if atoms[subbond].atomicNumber == 1:
                                nHyd += 1

                        if nHyd <= 1:
                            nucleotideAtoms["in8"] = bond
                        elif nHyd == 2:
                            nucleotideAtoms["in7"] = bond

            # Find C9 carbon connected to N8
            if nucleotideAtoms["in8"] != 0:
                for bond in atoms[nucleotideAtoms["in8"]].bonds:
                    if bond != nucleotideAtoms["ic7"] and atoms[bond].atomicNumber == 6:
                        nucleotideAtoms["ic9"] = bond

            # Find N9 or O9 to determine if adenine or guanine
            if nucleotideAtoms["ic9"] != 0:
                for bond in atoms[nucleotideAtoms["ic9"]].bonds:
                    if bond != nucleotideAtoms["in8"]:
                        if atoms[bond].atomicNumber == 7:
                            nucleotideAtoms["in9"] = bond
                        elif atoms[bond].atomicNumber == 8:
                            nucleotideAtoms["io9"] = bond

            # Set final nucleotide label
            if nucleotideAtoms["io9"] != 0:
                label = "  G" if not deoxy else " DG"
            elif nucleotideAtoms["in9"] != 0:
                label = "  A" if not deoxy else " DA"

        # For pyrimidines: uracil (U), cytosine (C), or thymine (T)
        elif baseType == "PYR":
            # Find N5 nitrogen attached to C6
            if nucleotideAtoms["ic6"] != 0:
                for bond in atoms[nucleotideAtoms["ic6"]].bonds:
                    if bond != nucleotideAtoms["in1"] and atoms[bond].atomicNumber == 7:
                        nucleotideAtoms["in5"] = bond

            # Find C4 carbon attached to N5
            if nucleotideAtoms["in5"] != 0:
                for bond in atoms[nucleotideAtoms["in5"]].bonds:
                    if bond != nucleotideAtoms["ic6"] and atoms[bond].atomicNumber == 6:
                        nucleotideAtoms["ic4"] = bond

            # Find O4 oxygen, N4 nitrogen, and C3 carbon from C4
            if nucleotideAtoms["ic4"] != 0:
                for bond in atoms[nucleotideAtoms["ic4"]].bonds:
                    if atoms[bond].atomicNumber == 6:
                        nucleotideAtoms["ic3"] = bond
                    elif (
                        bond != nucleotideAtoms["in5"] and atoms[bond].atomicNumber == 7
                    ):
                        nucleotideAtoms["in4"] = bond
                    elif atoms[bond].atomicNumber == 8:
                        nucleotideAtoms["io4"] = bond

            # Default to uracil
            label = "  U"

            # Check if cytosine (has N-H)
            if (
                nucleotideAtoms["in5"] != 0
                and len(atoms[nucleotideAtoms["in5"]].bonds) == 2
            ):
                # Cytosine identified by N5 with 2 bonds
                label = "  C" if not deoxy else " DC"

            # Check if thymine (has methyl group)
            if nucleotideAtoms["ic3"] != 0:
                for bond in atoms[nucleotideAtoms["ic3"]].bonds:
                    if (
                        bond != nucleotideAtoms["ic2"]
                        and bond != nucleotideAtoms["ic4"]
                        and atoms[bond].atomicNumber == 6
                    ):
                        # Additional verification for methyl group: should have 3 hydrogens
                        methylHydrogens = 0
                        for methylBond in atoms[bond].bonds:
                            if atoms[methylBond].atomicNumber == 1:
                                methylHydrogens += 1

                        # Only identify as thymine if the methyl group has appropriate hydrogens
                        if (
                            methylHydrogens >= 2
                        ):  # Should be 3, but allow some flexibility
                            nucleotideAtoms["icm"] = bond
                            label = "  T" if not deoxy else " DT"
        else:
            raise ValueError(
                f"Could not identify base type for nucleotide with bitorsion {ia}, {ib}, {ic}, {id}, {ie}"
            )

        return label, nucleotideAtoms

    @staticmethod
    def _collectNucleotideAtoms(
        nucleotideAtoms: Dict[str, Any], atoms: List[TinkerAtom]
    ) -> List[int]:
        """
        Collect all atoms in a nucleotide.

        Parameters
        ----------
        nucleotideAtoms : Dict[str, Any]
            Dictionary of nucleotide atoms
        atoms : List[TinkerAtom]
            List of atoms in the system

        Returns
        -------
        List[int]
            List of all atom indices in the nucleotide
        """
        residueAtoms = set()

        # Add all heavy atoms in the nucleotide
        for key, index in nucleotideAtoms.items():
            # Special handling for terminal oxygens (5' and 3') which can have indices of 0
            if key in ("io5s", "io3s"):
                residueAtoms.add(index)
            elif index != 0:
                residueAtoms.add(index)

        # Add all hydrogen atoms
        for atom in list(residueAtoms):
            for bond in atoms[atom].bonds:
                if atoms[bond].atomicNumber == 1:
                    residueAtoms.add(bond)

        return list(residueAtoms)

    # ------------------------------------------------------------------------------------------ #
    #                                     POSITIONS AND BOX                                      #
    # ------------------------------------------------------------------------------------------ #
    def getPositions(self, asNumpy: bool = False) -> Union[List[Vec3], np.ndarray]:
        """
        Get the atomic positions.

        Parameters
        ----------
        asNumpy : bool=False
            If true, the values are returned as a numpy array instead of a list of Vec3s.

        Returns
        -------
        list of Vec3 or np.ndarray
            The atomic positions
        """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = (
                    np.array(self.positions.value_in_unit(nanometers)) * nanometers
                )
            return self._numpyPositions
        return self.positions

    def getBoxVectors(self, asNumpy: bool = False) -> Union[List[Vec3], np.ndarray]:
        """
        Get the periodic box vectors.

        Parameters
        ----------
        asNumpy : bool, optional, default=False
            If true, the values are returned as a numpy array instead of a list of Vec3s.

        Returns
        -------
        list of Vec3 or np.ndarray
            The periodic box vectors.
        """
        if self.boxVectors is None:
            raise AttributeError(f"Box information not found in {self.file}")
        if asNumpy:
            if self._numpyBoxVectors is None:
                self._numpyBoxVectors = []
                self._numpyBoxVectors.append(
                    np.array(self.boxVectors[0].value_in_unit(nanometers)) * nanometers
                )
                self._numpyBoxVectors.append(
                    np.array(self.boxVectors[1].value_in_unit(nanometers)) * nanometers
                )
                self._numpyBoxVectors.append(
                    np.array(self.boxVectors[2].value_in_unit(nanometers)) * nanometers
                )
            return self._numpyBoxVectors
        return self.boxVectors

    # ------------------------------------------------------------------------------------------ #
    #                                   XYZ FILE PARSING                                         #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _loadXyzFile(file: str) -> Tuple[List[TinkerAtom], List[Vec3], List[Vec3]]:
        """
        Load and parse a Tinker .xyz file.

        Parameters
        ----------
        file : str
            Path to the .xyz file to load.

        Returns
        -------
        Tuple[List[TinkerAtom], List[Vec3], List[Vec3]]
            A tuple containing:
            - List of TinkerAtom objects
            - List of atomic positions
            - List of periodic box vectors (if present, otherwise None)
        """
        atoms = []
        boxVectors = None

        try:
            with open(file, "r") as f:
                firstLine = f.readline().strip()
                if not firstLine:
                    raise ValueError("Empty .xyz file")

                try:
                    nAtoms = int(firstLine.split()[0])
                    if nAtoms <= 0:
                        raise ValueError(
                            f"Invalid number of atoms: {nAtoms}. Must be positive."
                        )
                except (ValueError, IndexError):
                    raise ValueError(
                        "First line must contain a positive integer number of atoms"
                    )

                linesLeft = nAtoms
                lineNum = 1

                secondLine = f.readline().strip()
                lineNum += 1
                if not secondLine:
                    raise ValueError("Missing atom coordinates in .xyz file")

                secondLineSplit = secondLine.split()

                if len(secondLineSplit) == 6 and secondLineSplit[0] != "1":
                    try:
                        box = [
                            float(val) * (0.1 if i < 3 else math.pi / 180.0)
                            for i, val in enumerate(secondLineSplit)
                        ]
                        if any(math.isnan(x) or math.isinf(x) for x in box):
                            raise ValueError(
                                "Box vectors contain invalid values (NaN or Inf)"
                            )
                        boxVectors = computePeriodicBoxVectors(*box)
                    except ValueError as e:
                        raise ValueError(
                            f"Line {lineNum}: Error parsing box vectors: {e}"
                        )
                else:
                    TinkerFiles._parseAndStoreXyzLine(secondLine, lineNum, atoms)
                    linesLeft -= 1

                for i in range(linesLeft):
                    lineNum += 1
                    atomLine = f.readline().strip()
                    if not atomLine:
                        raise ValueError(
                            f"Expected {nAtoms} atoms but found only {i + (1 if len(secondLineSplit) != 6 or secondLineSplit[0] == '1' else 0)}"
                        )
                    TinkerFiles._parseAndStoreXyzLine(atomLine, lineNum, atoms)

                if len(atoms) != nAtoms:
                    raise ValueError(f"Expected {nAtoms} atoms but found {len(atoms)}")

            positions = [atom.positions for atom in atoms]
            return atoms, boxVectors, positions

        except FileNotFoundError:
            raise FileNotFoundError(f"XYZ file not found: {file}")
        except IOError as e:
            raise IOError(f"Error reading XYZ file {file}: {e}")
        except Exception as e:
            raise ValueError(f"Error parsing XYZ file {file}: {e}")

    @staticmethod
    def _parseAndStoreXyzLine(line: str, lineNum: int, atoms: List[TinkerAtom]) -> None:
        """
        Parse a single line from a TINKER .xyz file and create an TinkerAtom object.

        Parameters
        ----------
        line : str
            The line to parse from the .xyz file.
        lineNum : int
            The line number in the file (for error reporting).
        atoms : List[TinkerAtom]
            List to store the created TinkerAtom object.
        """
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(
                f"Line {lineNum}: Each line in the TINKER .xyz file must have at least 6 fields"
            )

        try:
            index = int(fields[0]) - 1
            if index < 0:
                raise ValueError(
                    f"Line {lineNum}: Invalid atom index {index}. Must be positive."
                )

            symbol = str(fields[1])
            if not symbol:
                raise ValueError(f"Line {lineNum}: Empty atom symbol")

            try:
                x = float(fields[2]) * 0.1
                y = float(fields[3]) * 0.1
                z = float(fields[4]) * 0.1
                if any(math.isnan(coord) or math.isinf(coord) for coord in (x, y, z)):
                    raise ValueError(
                        "Atom coordinates contain invalid values (NaN or Inf)"
                    )
                position = Vec3(x, y, z)
            except ValueError as e:
                raise ValueError(f"Line {lineNum}: Error parsing atom coordinates: {e}")

            atomType = str(fields[5])
            if not atomType:
                raise ValueError(f"Line {lineNum}: Empty atom type")

            try:
                bonds = [int(bond) - 1 for bond in fields[6:]]
                if any(bond < 0 for bond in bonds):
                    raise ValueError("Invalid bond index (must be >= 0)")
            except ValueError as e:
                raise ValueError(f"Line {lineNum}: Error parsing bond indices: {e}")

            # Create TinkerAtom object
            atom = TinkerAtom(
                symbol=symbol,
                positions=position,
                bonds=bonds,
                index=index,
                atomType=atomType,
            )
            atoms.append(atom)

        except ValueError as e:
            raise ValueError(f"Line {lineNum}: {str(e)}")

    # ------------------------------------------------------------------------------------------ #
    #                                   KEY FILE PARSING                                         #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _addAtomType(
        atomTypeDict: Dict[str, TinkerAtomType],
        atomType: str,
        atomClass: str,
        nameShort: str,
        nameLong: str,
        atomicNumber: int,
        mass: float,
        valence: int,
        element: str,
    ) -> None:
        """
        Add and validate atom type information.

        Parameters
        ----------
        atomTypeDict : Dict[str, TinkerAtomType]
            The dictionary of atom type data.
        atomType : str
            The atom type.
        atomClass : str
            The atom class.
        nameShort : str
            The short name of the atom type.
        nameLong : str
            The long name of the atom type.
        atomicNumber : int
            The atomic number of the atom type.
        mass : float
            The mass of the atom type.
        valence : int
            The valence of the atom type.
        element : str
            The element of the atom type.

        Raises
        ------
        ValueError
            If there is an inconsistency with existing atom type data.
        """
        tinkerType = TinkerAtomType(
            atomType=atomType,
            atomClass=atomClass,
            nameShort=nameShort,
            nameLong=nameLong,
            element=element,
            atomicNumber=atomicNumber,
            mass=mass,
            valence=valence,
        )

        if atomType in atomTypeDict:
            # Validate against existing atom type data
            stored = atomTypeDict[atomType]
            if stored != tinkerType:
                raise ValueError(
                    f"Inconsistent data for atom type '{atomType}': "
                    f"expected '{stored}', got '{tinkerType}'."
                )

        # Add new atom type to the dictionary
        atomTypeDict[atomType] = tinkerType

    @staticmethod
    def _loadKeyFile(
        keyFile: str,
    ) -> Tuple[Dict[str, TinkerAtomType], Dict[str, Dict], Dict[str, str]]:
        """
        Load a TINKER .key or .prm file.

        Parameters
        ----------
        keyFile : str
            The path to the .key or .prm file.

        Returns
        -------
        Tuple[Dict[str, TinkerAtomType], Dict[str, Dict], Dict[str, str]]
            A tuple containing:
            - Atom types dictionary
            - Forces dictionary
            - Scalars dictionary
        """
        atomTypesDict = dict()
        forcesDict = dict()
        scalarsDict = TinkerFiles.RECOGNIZED_SCALARS.copy()

        try:
            with open(keyFile, "r") as file:
                allLines = []
                for line_num, line in enumerate(file, 1):
                    try:
                        if line.count('"') % 2 != 0:
                            # Skip lines with an odd number of quotes to avoid parsing errors
                            # with citations or other non-essential information
                            continue

                        lineStripped = line.lstrip()
                        if lineStripped.startswith("#") or lineStripped == "":
                            continue

                        fields = shlex.split(line)
                        if fields:  # Make sure the line has content after parsing
                            allLines.append(fields)
                    except Exception as e:
                        raise ValueError(f"Line {line_num}: Error parsing line: {e}")

            lineIndex = 0
            while lineIndex < len(allLines):
                fields = allLines[lineIndex]
                if not fields:  # Skip empty lines
                    lineIndex += 1
                    continue

                if fields[0] == "atom":
                    if len(fields) != 8:
                        raise ValueError(
                            f"Invalid atom line: Expected 8 fields, got {len(fields)}. Fields: {fields}"
                        )
                    # Atom type information
                    # atom atomType atomClass nameShort nameLong atomicNumber mass valence
                    (
                        atomType,
                        atomClass,
                        nameShort,
                        nameLong,
                        atomicNumber,
                        mass,
                        valence,
                    ) = fields[1:]

                    try:
                        atomicNumber = int(atomicNumber)
                        mass = float(mass)
                        valence = int(valence)
                    except ValueError:
                        raise ValueError(
                            f"Invalid numeric values in atom line: {fields}"
                        )

                    if atomicNumber > 0:
                        element = elem.Element.getByAtomicNumber(atomicNumber).symbol
                    else:
                        element = None  # For dummy atoms

                    nameLong = re.sub(r"\s+", " ", nameLong.strip())
                    TinkerFiles._addAtomType(
                        atomTypesDict,
                        atomType,
                        atomClass,
                        nameShort,
                        nameLong,
                        atomicNumber,
                        mass,
                        valence,
                        element,
                    )
                    lineIndex += 1
                elif fields[0] == "biotype":
                    # No need to read biotype information because the Tinker
                    # .xyz file directly contains atom type information
                    lineIndex += 1
                elif fields[0] in TinkerFiles.RECOGNIZED_FORCES:
                    if TinkerFiles.RECOGNIZED_FORCES[fields[0]] == 1:
                        if fields[0] not in forcesDict:
                            forcesDict[fields[0]] = []
                        forcesDict[fields[0]].append(fields[1:])
                        lineIndex += 1
                    else:
                        # Call the function to parse the specific force
                        lineIndex = TinkerFiles.RECOGNIZED_FORCES[fields[0]](
                            lineIndex, allLines, forcesDict
                        )
                elif fields[0] in TinkerFiles.RECOGNIZED_SCALARS:
                    scalar, value = fields
                    scalarsDict[scalar] = value
                    lineIndex += 1
                else:
                    # Skip unrecognized fields
                    lineIndex += 1

            return atomTypesDict, forcesDict, scalarsDict

        except FileNotFoundError:
            raise FileNotFoundError(f"File {keyFile} not found")
        except IOError as e:
            raise IOError(f"Error reading file {keyFile}: {e}")
        except Exception as e:
            raise ValueError(f"Error parsing {keyFile}: {e}")


class TinkerXmlWriter:
    """
    Utility class for converting Tinker force field parameters to OpenMM XML format.

    This class provides methods to transform Tinker force field parameters extracted
    from .key/.prm files into OpenMM-compatible XML force field files. It creates
    both standard and implicit solvent (Generalized Kirkwood) parameter files.
    """

    @staticmethod
    def _writeAtomTypes(
        root: etree.Element, atomTypes: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Write the atom types to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        atomTypes : dict
            The atom types dictionary.
        """
        atomTypesElement = etree.SubElement(root, "AtomTypes")
        for atomType in sorted(atomTypes, key=int):
            atomTypeElement = etree.SubElement(atomTypesElement, "Type")
            atomTypeElement.attrib["name"] = atomType
            atomTypeElement.attrib["class"] = atomTypes[atomType]["atomClass"]
            atomTypeElement.attrib["element"] = atomTypes[atomType]["element"]
            atomTypeElement.attrib["mass"] = str(atomTypes[atomType]["mass"])

    @staticmethod
    def _writeResidues(
        root: etree.Element,
        atomDict: Dict[int, TinkerAtom],
        topology: top.Topology,
    ) -> None:
        """
        Write the residues to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        atomDict : Dict[int, TinkerAtom]
            Dictionary mapping atom indices to TinkerAtom objects.
        topology : openmm.app.topology.Topology
            The topology object.
        """
        residuesElement = etree.SubElement(root, "Residues")
        residuesSeen = set()
        for residue in topology.residues():
            if residue.name not in residuesSeen:
                residueElement = etree.SubElement(residuesElement, "Residue")
                residueElement.attrib["name"] = residue.name

                if residue.name == "HOH":
                    # AMOEBA water is flexible
                    residueElement.attrib["rigidWater"] = "false"

                for i, atom in enumerate(residue.atoms()):
                    atomElement = etree.SubElement(residueElement, "Atom")
                    atomElement.attrib["name"] = atomDict[atom.index].nameShort + str(i)
                    atomElement.attrib["type"] = atomDict[atom.index].atomType

                baseIndex = next(atom.index for atom in residue.atoms())
                for bond in residue.bonds():
                    bondElement = etree.SubElement(residueElement, "Bond")
                    bondElement.attrib["from"] = str(bond[0].index - baseIndex)
                    bondElement.attrib["to"] = str(bond[1].index - baseIndex)

                residuesSeen.add(residue.name)

    @staticmethod
    def _writeAmoebaBondForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaBondForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "bond" not in forces:
            return

        cubic = 10.0 * float(scalars["bond-cubic"])
        quartic = 100.0 * float(scalars["bond-quartic"])
        bondForceElement = etree.SubElement(root, "AmoebaBondForce")
        bondForceElement.attrib["bond-cubic"] = str(cubic)
        bondForceElement.attrib["bond-quartic"] = str(quartic)

        for bond in forces["bond"]:
            length = float(bond[3]) * 0.1
            k = float(bond[2]) * 100.0 * 4.184
            bondElement = etree.SubElement(bondForceElement, "Bond")
            bondElement.attrib["class1"] = bond[0]
            bondElement.attrib["class2"] = bond[1]
            bondElement.attrib["length"] = str(length)
            bondElement.attrib["k"] = str(k)

    @staticmethod
    def _writeAmoebaAngleForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaAngleForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "angle" not in forces and "anglep" not in forces:
            return

        cubic = float(scalars["angle-cubic"])
        quartic = float(scalars["angle-quartic"])
        pentic = float(scalars["angle-pentic"])
        sextic = float(scalars["angle-sextic"])
        angleForceElement = etree.SubElement(root, "AmoebaAngleForce")
        angleForceElement.attrib["angle-cubic"] = str(cubic)
        angleForceElement.attrib["angle-quartic"] = str(quartic)
        angleForceElement.attrib["angle-pentic"] = str(pentic)
        angleForceElement.attrib["angle-sextic"] = str(sextic)

        radian = 57.2957795130
        radian2 = 4.184 / (radian * radian)

        for angleSet in ["angle", "anglep"]:
            if angleSet not in forces:
                continue
            for angle in forces[angleSet]:
                k = float(angle[3]) * radian2
                angleElement = etree.SubElement(angleForceElement, "Angle")
                angleElement.attrib["class1"] = angle[0]
                angleElement.attrib["class2"] = angle[1]
                angleElement.attrib["class3"] = angle[2]
                angleElement.attrib["k"] = str(k)
                angleElement.attrib["angle1"] = angle[4]
                if len(angle) > 5:
                    angleElement.attrib["angle2"] = angle[5]
                if len(angle) > 6:
                    angleElement.attrib["angle3"] = angle[6]
                angleElement.attrib["inPlane"] = (
                    "true" if angleSet == "anglep" else "false"
                )

    @staticmethod
    def _writeAmoebaOutOfPlaneBendForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaOutOfPlaneBendForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "opbend" not in forces:
            return

        cubic = float(scalars["opbend-cubic"])
        quartic = float(scalars["opbend-quartic"])
        pentic = float(scalars["opbend-pentic"])
        sextic = float(scalars["opbend-sextic"])
        opbendType = scalars["opbendtype"]

        opbendForceElement = etree.SubElement(root, "AmoebaOutOfPlaneBendForce")
        opbendForceElement.attrib["type"] = opbendType
        opbendForceElement.attrib["opbend-cubic"] = str(cubic)
        opbendForceElement.attrib["opbend-quartic"] = str(quartic)
        opbendForceElement.attrib["opbend-pentic"] = str(pentic)
        opbendForceElement.attrib["opbend-sextic"] = str(sextic)

        radian = 57.2957795130
        radian2 = 4.184 / (radian * radian)

        for opbend in forces["opbend"]:
            k = float(opbend[4]) * radian2
            angleElement = etree.SubElement(opbendForceElement, "Angle")
            angleElement.attrib["class1"] = opbend[0]
            angleElement.attrib["class2"] = opbend[1]
            angleElement.attrib["class3"] = opbend[2]
            angleElement.attrib["class4"] = opbend[3]
            angleElement.attrib["k"] = str(k)

    @staticmethod
    def _writeAmoebaTorsionForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "torsion" not in forces:
            return

        torsionUnit = float(scalars["torsionunit"])
        conversion = 4.184 * torsionUnit
        radian = 57.2957795130

        torsionForceElement = etree.SubElement(root, "PeriodicTorsionForce")
        for torsion in forces["torsion"]:
            torsionElement = etree.SubElement(torsionForceElement, "Proper")
            torsionElement.attrib["class1"] = torsion[0]
            torsionElement.attrib["class2"] = torsion[1]
            torsionElement.attrib["class3"] = torsion[2]
            torsionElement.attrib["class4"] = torsion[3]
            startIndex = 4
            for ii in range(0, (len(torsion) - 4) // 3):
                torsionSuffix = str(ii + 1)
                amplitudeAttributeName = f"k{torsionSuffix}"
                angleAttributeName = f"phase{torsionSuffix}"
                periodicityAttributeName = f"periodicity{torsionSuffix}"
                amplitude = float(torsion[startIndex]) * conversion
                angle = float(torsion[startIndex + 1]) / radian
                periodicity = int(torsion[startIndex + 2])
                torsionElement.attrib[amplitudeAttributeName] = str(amplitude)
                torsionElement.attrib[angleAttributeName] = str(angle)
                torsionElement.attrib[periodicityAttributeName] = str(periodicity)
                startIndex += 3

    @staticmethod
    def _writeAmoebaPiTorsionForce(root: etree.Element, forces: Dict[str, Any]) -> None:
        """
        Write the AmoebaPiTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "pitors" not in forces:
            return

        piTorsionUnit = 1.0
        conversion = 4.184 * piTorsionUnit
        piTorsionForceElement = etree.SubElement(root, "AmoebaPiTorsionForce")
        piTorsionForceElement.attrib["piTorsionUnit"] = str(piTorsionUnit)
        for piTorsion in forces["pitors"]:
            k = float(piTorsion[2]) * conversion
            piTorsionElement = etree.SubElement(piTorsionForceElement, "PiTorsion")
            piTorsionElement.attrib["class1"] = piTorsion[0]
            piTorsionElement.attrib["class2"] = piTorsion[1]
            piTorsionElement.attrib["k"] = str(k)

    @staticmethod
    def _writeAmoebaStretchTorsionForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaStretchTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "strtors" not in forces:
            return

        stretchTorsionForceElement = etree.SubElement(root, "AmoebaStretchTorsionForce")
        for torsion in forces["strtors"]:
            v = [float(x) * 10 * 4.184 for x in torsion[4:]]
            torsionElement = etree.SubElement(stretchTorsionForceElement, "Torsion")
            torsionElement.attrib["class1"] = torsion[0]
            torsionElement.attrib["class2"] = torsion[1]
            torsionElement.attrib["class3"] = torsion[2]
            torsionElement.attrib["class4"] = torsion[3]
            torsionElement.attrib["v11"] = str(v[0])
            torsionElement.attrib["v12"] = str(v[1])
            torsionElement.attrib["v13"] = str(v[2])
            torsionElement.attrib["v21"] = str(v[3])
            torsionElement.attrib["v22"] = str(v[4])
            torsionElement.attrib["v23"] = str(v[5])
            torsionElement.attrib["v31"] = str(v[6])
            torsionElement.attrib["v32"] = str(v[7])
            torsionElement.attrib["v33"] = str(v[8])

    @staticmethod
    def _writeAmoebaAngleTorsionForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaAngleTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "angtors" not in forces:
            return

        angleTorsionForceElement = etree.SubElement(root, "AmoebaAngleTorsionForce")
        for torsion in forces["angtors"]:
            v = [float(x) * 4.184 for x in torsion[4:]]
            torsionElement = etree.SubElement(angleTorsionForceElement, "Torsion")
            torsionElement.attrib["class1"] = torsion[0]
            torsionElement.attrib["class2"] = torsion[1]
            torsionElement.attrib["class3"] = torsion[2]
            torsionElement.attrib["class4"] = torsion[3]
            torsionElement.attrib["v11"] = str(v[0])
            torsionElement.attrib["v12"] = str(v[1])
            torsionElement.attrib["v13"] = str(v[2])
            torsionElement.attrib["v21"] = str(v[3])
            torsionElement.attrib["v22"] = str(v[4])
            torsionElement.attrib["v23"] = str(v[5])

    @staticmethod
    def _writeAmoebaStretchBendForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaStretchBendForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "strbnd" not in forces:
            return

        stretchBendUnit = 1.0
        radian = 57.2957795130
        conversion = 41.84 / radian
        stretchBendForceElement = etree.SubElement(root, "AmoebaStretchBendForce")
        stretchBendForceElement.attrib["stretchBendUnit"] = str(stretchBendUnit)

        for stretchBend in forces["strbnd"]:
            k1 = float(stretchBend[3]) * conversion
            k2 = float(stretchBend[4]) * conversion
            stretchBendElement = etree.SubElement(
                stretchBendForceElement, "StretchBend"
            )
            stretchBendElement.attrib["class1"] = stretchBend[0]
            stretchBendElement.attrib["class2"] = stretchBend[1]
            stretchBendElement.attrib["class3"] = stretchBend[2]
            stretchBendElement.attrib["k1"] = str(k1)
            stretchBendElement.attrib["k2"] = str(k2)

    @staticmethod
    def _writeAmoebaTorsionTorsionForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaTorsionTorsionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "tortors" not in forces:
            return

        torsionTorsionUnit = 1.0

        torsionTorsionForceElement = etree.SubElement(root, "AmoebaTorsionTorsionForce")
        for index, torsionTorsion in enumerate(forces["tortors"]):
            torsionTorsionElement = etree.SubElement(
                torsionTorsionForceElement, "TorsionTorsion"
            )
            torsionTorsionElement.attrib["class1"] = torsionTorsion[0][0]
            torsionTorsionElement.attrib["class2"] = torsionTorsion[0][1]
            torsionTorsionElement.attrib["class3"] = torsionTorsion[0][2]
            torsionTorsionElement.attrib["class4"] = torsionTorsion[0][3]
            torsionTorsionElement.attrib["class5"] = torsionTorsion[0][4]
            torsionTorsionElement.attrib["grid"] = str(index)
            torsionTorsionElement.attrib["nx"] = str(torsionTorsion[0][5])
            torsionTorsionElement.attrib["ny"] = str(torsionTorsion[0][6])

        for index, torsionTorsion in enumerate(forces["tortors"]):
            gridElement = etree.SubElement(
                torsionTorsionForceElement, "TorsionTorsionGrid"
            )
            gridElement.attrib["grid"] = str(index)
            gridElement.attrib["nx"] = str(torsionTorsion[0][5])
            gridElement.attrib["ny"] = str(torsionTorsion[0][6])
            for gridIndex, gridEntry in enumerate(torsionTorsion[1]):
                if len(gridEntry) > 5:
                    f = float(gridEntry[2]) * 4.184
                    fx = float(gridEntry[3]) * 4.184
                    fy = float(gridEntry[4]) * 4.184
                    fxy = float(gridEntry[5]) * 4.184
                    gridEntryElement = etree.SubElement(gridElement, "Grid")
                    gridEntryElement.attrib["angle1"] = gridEntry[0]
                    gridEntryElement.attrib["angle2"] = gridEntry[1]
                    gridEntryElement.attrib["f"] = str(f)
                    gridEntryElement.attrib["fx"] = str(fx)
                    gridEntryElement.attrib["fy"] = str(fy)
                    gridEntryElement.attrib["fxy"] = str(fxy)
                elif len(gridEntry) > 2:
                    f = float(gridEntry[2]) * 4.184
                    gridEntryElement = etree.SubElement(gridElement, "Grid")
                    gridEntryElement.attrib["angle1"] = gridEntry[0]
                    gridEntryElement.attrib["angle2"] = gridEntry[1]
                    gridEntryElement.attrib["f"] = str(f)

    @staticmethod
    def _writeAmoebaVdwForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaVdwForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "vdw" not in forces and "vdwpr" not in forces:
            return

        vdwForceElement = etree.SubElement(root, "AmoebaVdwForce")
        vdwForceElement.attrib["type"] = scalars["vdwtype"]
        vdwForceElement.attrib["radiusrule"] = scalars["radiusrule"]
        vdwForceElement.attrib["radiustype"] = scalars["radiustype"]
        vdwForceElement.attrib["radiussize"] = scalars["radiussize"]
        vdwForceElement.attrib["epsilonrule"] = scalars["epsilonrule"]
        vdwForceElement.attrib["vdw-13-scale"] = scalars["vdw-13-scale"]
        vdwForceElement.attrib["vdw-14-scale"] = scalars["vdw-14-scale"]
        vdwForceElement.attrib["vdw-15-scale"] = scalars["vdw-15-scale"]

        if "vdw" in forces:
            for vdw in forces["vdw"]:
                sigma = float(vdw[1]) * 0.1
                epsilon = float(vdw[2]) * 4.184
                reduction = vdw[3] if len(vdw) > 3 else 1.0
                vdwElement = etree.SubElement(vdwForceElement, "Vdw")
                vdwElement.attrib["class"] = vdw[0]
                vdwElement.attrib["sigma"] = str(sigma)
                vdwElement.attrib["epsilon"] = str(epsilon)
                vdwElement.attrib["reduction"] = str(reduction)

        if "vdwpr" in forces:
            for pair in forces["vdwpr"]:
                sigma = float(pair[2]) * 0.1
                epsilon = float(pair[3]) * 4.184
                pairElement = etree.SubElement(vdwForceElement, "Pair")
                pairElement.attrib["class1"] = pair[0]
                pairElement.attrib["class2"] = pair[1]
                pairElement.attrib["sigma"] = str(sigma)
                pairElement.attrib["epsilon"] = str(epsilon)

    @staticmethod
    def _writeAmoebaMultipoleForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaMultipoleForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "multipole" not in forces and "polarize" not in forces:
            return

        amoebaMultipoleForceElement = etree.SubElement(root, "AmoebaMultipoleForce")
        amoebaMultipoleForceElement.attrib["mpole12Scale"] = scalars["mpole-12-scale"]
        amoebaMultipoleForceElement.attrib["mpole13Scale"] = scalars["mpole-13-scale"]
        amoebaMultipoleForceElement.attrib["mpole14Scale"] = scalars["mpole-14-scale"]
        amoebaMultipoleForceElement.attrib["mpole15Scale"] = scalars["mpole-15-scale"]
        amoebaMultipoleForceElement.attrib["polar12Scale"] = scalars["polar-12-scale"]
        amoebaMultipoleForceElement.attrib["polar13Scale"] = scalars["polar-13-scale"]
        amoebaMultipoleForceElement.attrib["polar14Scale"] = scalars["polar-14-scale"]
        amoebaMultipoleForceElement.attrib["polar15Scale"] = scalars["polar-15-scale"]
        amoebaMultipoleForceElement.attrib["polar14Intra"] = scalars["polar-14-intra"]
        amoebaMultipoleForceElement.attrib["direct11Scale"] = scalars["direct-11-scale"]
        amoebaMultipoleForceElement.attrib["direct12Scale"] = scalars["direct-12-scale"]
        amoebaMultipoleForceElement.attrib["direct13Scale"] = scalars["direct-13-scale"]
        amoebaMultipoleForceElement.attrib["direct14Scale"] = scalars["direct-14-scale"]
        amoebaMultipoleForceElement.attrib["mutual11Scale"] = scalars["mutual-11-scale"]
        amoebaMultipoleForceElement.attrib["mutual12Scale"] = scalars["mutual-12-scale"]
        amoebaMultipoleForceElement.attrib["mutual13Scale"] = scalars["mutual-13-scale"]
        amoebaMultipoleForceElement.attrib["mutual14Scale"] = scalars["mutual-14-scale"]

        bohr = 0.52917720859
        dipoleConversion = 0.1 * bohr
        quadrupoleConversion = 0.01 * bohr * bohr / 3.0
        for multipoleInfo in forces["multipole"]:
            axisInfo = multipoleInfo[0]

            if int(axisInfo[0]) < 0:
                # Some Tinker files include a multipole parameter definition
                # for atom types with negative indices with no multipole parameters (all zeros).
                # We skip these entries as they cause an error in OpenMM.
                continue

            multipoles = multipoleInfo[1]
            multipoleElement = etree.SubElement(
                amoebaMultipoleForceElement, "Multipole"
            )
            multipoleElement.attrib["type"] = axisInfo[0]

            axisInfoLen = len(axisInfo)

            if axisInfoLen > 1:
                multipoleElement.attrib["kz"] = axisInfo[1]
            if axisInfoLen > 2:
                multipoleElement.attrib["kx"] = axisInfo[2]
            if axisInfoLen > 3:
                multipoleElement.attrib["ky"] = axisInfo[3]

            multipoleElement.attrib["c0"] = multipoles[0]
            multipoleElement.attrib["d1"] = str(dipoleConversion * float(multipoles[1]))
            multipoleElement.attrib["d2"] = str(dipoleConversion * float(multipoles[2]))
            multipoleElement.attrib["d3"] = str(dipoleConversion * float(multipoles[3]))
            multipoleElement.attrib["q11"] = str(
                quadrupoleConversion * float(multipoles[4])
            )
            multipoleElement.attrib["q21"] = str(
                quadrupoleConversion * float(multipoles[5])
            )
            multipoleElement.attrib["q22"] = str(
                quadrupoleConversion * float(multipoles[6])
            )
            multipoleElement.attrib["q31"] = str(
                quadrupoleConversion * float(multipoles[7])
            )
            multipoleElement.attrib["q32"] = str(
                quadrupoleConversion * float(multipoles[8])
            )
            multipoleElement.attrib["q33"] = str(
                quadrupoleConversion * float(multipoles[9])
            )

        polarityConversion = 0.001
        for polarize in forces["polarize"]:
            polarizeElement = etree.SubElement(amoebaMultipoleForceElement, "Polarize")
            polarizeElement.attrib["type"] = polarize[0]
            polarizeElement.attrib["polarizability"] = str(
                polarityConversion * float(polarize[1])
            )
            polarizeElement.attrib["thole"] = polarize[2]
            for ii in range(3, len(polarize)):
                polarizeElement.attrib[f"pgrp{ii - 2}"] = polarize[ii]

    @staticmethod
    def _writeAmoebaGeneralizedKirkwoodForce(
        root: etree.Element,
        forces: Dict[str, Any],
        atomTypes: Dict[str, Dict[str, Any]],
    ) -> None:
        """
        Write the AmoebaGeneralizedKirkwoodForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        atomTypes : dict
            The atom types dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "multipole" not in forces:
            return

        solventDielectric = 78.3
        soluteDielectric = 1.0
        includeCavityTerm = 1
        probeRadius = 0.14
        surfaceAreaFactor = -6.0 * 3.1415926535 * 0.0216 * 1000.0 * 0.4184

        gkForceElement = etree.SubElement(root, "AmoebaGeneralizedKirkwoodForce")
        gkForceElement.attrib["solventDielectric"] = str(solventDielectric)
        gkForceElement.attrib["soluteDielectric"] = str(soluteDielectric)
        gkForceElement.attrib["includeCavityTerm"] = str(includeCavityTerm)
        gkForceElement.attrib["probeRadius"] = str(probeRadius)
        gkForceElement.attrib["surfaceAreaFactor"] = str(surfaceAreaFactor)

        for multipoleInfo in forces["multipole"]:
            axisInfo = multipoleInfo[0]
            multipoles = multipoleInfo[1]
            atomType = axisInfo[0]

            shct = 0.8  # Default value for overlap scale factor
            if atomType in atomTypes:
                element = atomTypes[atomType]["element"]
                if element == "H":
                    shct = 0.85
                elif element == "C":
                    shct = 0.72
                elif element == "N":
                    shct = 0.79
                elif element == "O":
                    shct = 0.85
                elif element == "F":
                    shct = 0.88
                elif element == "P":
                    shct = 0.86
                elif element == "S":
                    shct = 0.96
                elif element == "Fe":
                    shct = 0.88
                else:
                    print(
                        f"Warning no overlap scale factor for type={atomType} element={element}"
                    )
            else:
                print(f"Warning no overlap scale factor for type={atomType} ")

            gkElement = etree.SubElement(gkForceElement, "GeneralizedKirkwood")
            gkElement.attrib["type"] = axisInfo[0]
            gkElement.attrib["charge"] = multipoles[0]
            gkElement.attrib["shct"] = str(shct)

    @staticmethod
    def _writeAmoebaWcaDispersionForce(
        root: etree.Element, forces: Dict[str, Any], scalars: Dict[str, str]
    ) -> None:
        """
        Write the AmoebaWcaDispersionForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        """
        if "vdw" not in forces:
            return

        epso = 0.1100
        epsh = 0.0135
        rmino = 1.7025
        rminh = 1.3275
        awater = 0.033428
        slevy = 1.0
        dispoff = 0.26
        shctd = 0.81

        root = etree.SubElement(root, "AmoebaWcaDispersionForce")
        root.attrib["epso"] = str(epso * 4.184)
        root.attrib["epsh"] = str(epsh * 4.184)
        root.attrib["rmino"] = str(rmino * 0.1)
        root.attrib["rminh"] = str(rminh * 0.1)
        root.attrib["awater"] = str(1000.0 * awater)
        root.attrib["slevy"] = str(slevy)
        root.attrib["dispoff"] = str(0.1 * dispoff)
        root.attrib["shctd"] = str(shctd)
        convert = 0.1

        if scalars["radiustype"] == "SIGMA":
            convert *= 1.122462048309372

        if scalars["radiussize"] == "DIAMETER":
            convert *= 0.5

        for vdw in forces["vdw"]:
            sigma = float(vdw[1])
            sigma *= convert
            epsilon = float(vdw[2]) * 4.184
            vdwElement = etree.SubElement(root, "WcaDispersion")
            vdwElement.attrib["class"] = vdw[0]
            vdwElement.attrib["radius"] = str(sigma)
            vdwElement.attrib["epsilon"] = str(epsilon)

    @staticmethod
    def _writeAmoebaUreyBradleyForce(
        root: etree.Element, forces: Dict[str, Any]
    ) -> None:
        """
        Write the AmoebaUreyBradleyForce to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        forces : dict
            The forces dictionary.
        """
        if "ureybrad" not in forces:
            return

        cubic = 0.0
        quartic = 0.0
        ureyBradleyForceElement = etree.SubElement(root, "AmoebaUreyBradleyForce")
        ureyBradleyForceElement.attrib["cubic"] = str(cubic)
        ureyBradleyForceElement.attrib["quartic"] = str(quartic)

        for ub in forces["ureybrad"]:
            k = float(ub[3]) * 4.184 * 100.0
            d = float(ub[4]) * 0.1
            ureyBradleyElement = etree.SubElement(
                ureyBradleyForceElement, "UreyBradley"
            )
            ureyBradleyElement.attrib["class1"] = ub[0]
            ureyBradleyElement.attrib["class2"] = ub[1]
            ureyBradleyElement.attrib["class3"] = ub[2]
            ureyBradleyElement.attrib["k"] = str(k)
            ureyBradleyElement.attrib["d"] = str(d)

    @staticmethod
    def writeXmlFiles(topology, filename, atomTypes, forces, scalars, atomDict):
        """
        Create XML force field files from Tinker parameters.

        Generates two XML files:
        1. Standard force field file (filename.xml)
        2. Implicit solvent version with Generalized Kirkwood model (filename_gk.xml)

        Parameters
        ----------
        topology : openmm.app.Topology
            The molecular topology.
        filename : str
            The Tinker .key file path used to derive output filenames.
        atomTypes : dict
            Dictionary of atom types and their properties.
        forces : dict
            Dictionary of force parameters from the Tinker file.
        scalars : dict
            Dictionary of scalar parameters from the Tinker file.
        atomDict : dict
            Dictionary mapping atom indices to TinkerAtom objects.

        Returns
        -------
        tuple(file, file)
            Tuple containing file handles for the standard and implicit solvent XML files.
        """
        root = etree.Element("ForceField")

        # Write the header
        info = etree.SubElement(root, "Info")
        etree.SubElement(info, "Source").text = os.path.basename(filename)
        etree.SubElement(info, "DateGenerated").text = datetime.date.today().isoformat()
        etree.SubElement(info, "Reference")

        # Copy the root element for the Generalized Kirkwood force
        gkRoot = root.__copy__()

        # Write the AtomTypes
        TinkerXmlWriter._writeAtomTypes(root, atomTypes)

        # Write the Residues
        TinkerXmlWriter._writeResidues(root, atomDict, topology)

        # Write the Forces
        # Bonded forces
        TinkerXmlWriter._writeAmoebaBondForce(root, forces, scalars)
        TinkerXmlWriter._writeAmoebaAngleForce(root, forces, scalars)
        TinkerXmlWriter._writeAmoebaOutOfPlaneBendForce(root, forces, scalars)
        TinkerXmlWriter._writeAmoebaTorsionForce(root, forces, scalars)
        TinkerXmlWriter._writeAmoebaPiTorsionForce(root, forces)
        TinkerXmlWriter._writeAmoebaStretchTorsionForce(root, forces)
        TinkerXmlWriter._writeAmoebaAngleTorsionForce(root, forces)
        TinkerXmlWriter._writeAmoebaStretchBendForce(root, forces)
        TinkerXmlWriter._writeAmoebaTorsionTorsionForce(root, forces)

        # Non-bonded forces
        TinkerXmlWriter._writeAmoebaVdwForce(root, forces, scalars)
        TinkerXmlWriter._writeAmoebaMultipoleForce(root, forces, scalars)
        TinkerXmlWriter._writeAmoebaUreyBradleyForce(root, forces)

        # Write the Generalized Kirkwood Force
        TinkerXmlWriter._writeAmoebaGeneralizedKirkwoodForce(gkRoot, forces, atomTypes)
        TinkerXmlWriter._writeAmoebaWcaDispersionForce(gkRoot, forces, scalars)

        # Write footer
        tree = etree.ElementTree(root)
        etree.indent(root, "  ")
        gkTree = etree.ElementTree(gkRoot)
        etree.indent(gkRoot, "  ")

        # Store the XML file in memory
        tinkerXmlFile = io.StringIO()
        tree.write(
            tinkerXmlFile, xml_declaration=True, encoding="unicode", method="xml"
        )

        gkTinkerXmlFile = io.StringIO()
        gkTree.write(
            gkTinkerXmlFile, xml_declaration=True, encoding="unicode", method="xml"
        )

        if scalars["forcefield"] == "-2.55":
            scalars["forcefield"] = filename.rsplit(".", 1)[0]

        tinkerXmlFileName = scalars["forcefield"] + ".xml"
        with open(tinkerXmlFileName, "w") as f:
            f.write(tinkerXmlFile.getvalue())

        gkTinkerXmlFileName = scalars["forcefield"] + "_gk.xml"
        with open(gkTinkerXmlFileName, "w") as f:
            f.write(gkTinkerXmlFile.getvalue())

        return tinkerXmlFile, gkTinkerXmlFile
