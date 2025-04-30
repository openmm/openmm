"""
tinkerfiles.py: Used for loading TINKER key/prm and xyz files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2025 Stanford University and the Authors.
Authors: Joao Morado
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

from __future__ import absolute_import

__author__ = "Joao Morado"

import datetime
import io
import math
import os
import re
import shlex
import xml.etree.ElementTree as etree
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np

from openmm.app.internal.unitcell import computePeriodicBoxVectors
from openmm.unit import nanometers
from openmm.vec3 import Vec3

from . import element as elem
from . import forcefield as ff
from . import topology as top


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
            "polarization": "mutual",
            "polar-eps": "0.00001",
            "ewald": "yes",
            "ewald-cutoff": "7.0",
            "mpole-cutoff": "9.0",
            "vdw-cutoff": "9.0",
            "radiusrule": "ARITHMETIC",
            "radiustype": "SIGMA",
            "radiussize": "DIAMETER",
            "epsilonrule": "GEOMETRIC",
            "vdwtype": "BUFFERED-14-7",
            "a-expterm": "3.43",
            "b-expterm": "2.0",
            "c-expterm": "0.0",
            "gamma": "0.45",
            "dielectric": "1.0",
        }

        return RECOGNIZED_FORCES, RECOGNIZED_SCALARS

    # Call the initialization method
    RECOGNIZED_FORCES, RECOGNIZED_SCALARS = _initialize_class()

    _AMINO_ACID_LIST = [
        "GLY",
        "ALA",
        "VAL",
        "LEU",
        "ILE",
        "SER",
        "THR",
        "CYS",
        "CYX",
        "CYD",
        "PRO",
        "PHE",
        "TYR",
        "TYD",
        "TRP",
        "HIS",
        "HID",
        "HIE",
        "ASP",
        "ASH",
        "ASN",
        "GLU",
        "GLH",
        "GLN",
        "MET",
        "LYS",
        "LYD",
        "ARG",
        "ORN",
        "AIB",
        "PCA",
        "H2N",
        "FOR",
        "ACE",
        "COH",
        "NH2",
        "NME",
        "UNK",
    ]

    # Nucleotide Codes
    _NUCLEOTIDE_LIST = {
        "  A": "A",
        "  G": "G",
        "  C": "C",
        "  U": "U",
        " DA": "D",
        " DG": "B",
        " DC": "I",
        " DT": "T",
        " MP": "1",
        " DP": "2",
        " TP": "3",
        "UNK": "X",
    }

    def __init__(
        self,
        xyz: str,
        key: str,
        seq: str = None,
        periodicBoxVectors: Tuple[Vec3, Vec3, Vec3] = None,
        unitCellDimensions: Vec3 = None,
        writeXmlFiles: bool = False,
    ):
        """
        Load exactly one .xyz file and one or more .key or .prm files.

        Parameters
        ----------
        xyz : str
            The path to the xyz file to load.
        key : str or list of str
            The path(s) to the key/prm file(s) to load.
        seq : str, optional
            The path to the .seq file to load.
        periodicBoxVectors : tuple of Vec3
            The vectors defining the periodic box.
            If provided, this overwrites the box information from the xyz file.
        unitCellDimensions : Vec3, optional
            The dimensions of the crystallographic unit cell.
            For non-rectangular unit cells, specify periodicBoxVectors instead.
            If provided, this overwrites the box information from the xyz file.
        writeXmlFiles : bool, optional, default=False
            If True, the residue and force field XML files are written to disk.
        """
        # ----------------------- INTERNAL VARIABLES -----------------------
        # Populate parser functions for the recognized forces
        self.RECOGNIZED_FORCES["tortors"] = TinkerFiles.__addTorTor
        self.RECOGNIZED_FORCES["multipole"] = TinkerFiles.__addMultipole

        # Store the input parameters
        self._writeXmlFiles = writeXmlFiles
        self._XmlFilesList = None

        # Internal variables to store the data from the files
        self._xyzDict = None
        self._seqDict = None
        self._atomDict = None
        self._scalars, self._forces, self._atomTypes, self._bioTypes = [], [], [], []

        # Topology
        self.topology = None

        # Position and box information
        self.positions = None
        self.boxVectors = None
        self._numpyPositions = None
        self._numpyBoxVectors = None

        # ----------------------- LOAD FILES -----------------------
        # Load the .xyz file
        self._xyzDict, self.boxVectors, self.positions = self._loadXyzFile(xyz)
        self.positions = self.positions * nanometers

        # Load the .key or .prm file(s)
        key = key if isinstance(key, list) else [key]
        for keyFile in key:
            atomTypes, bioTypes, forces, scalars = self._loadKeyFile(keyFile)
            self._atomTypes.append(atomTypes)
            self._bioTypes.append(bioTypes)
            self._forces.append(forces)
            self._scalars.append(scalars)

        # Load the .seq file
        if seq is not None:
            self._seqDict = self._loadSeqFile(seq)
        else:
            self._seqDict = None

        # ----------------------- COMBINE DATA -----------------------
        # Combine the data from the .xyz and .key files
        self._atomDict = TinkerFiles._combineXyzAndKeyData(
            self._xyzDict, self._atomTypes, self._bioTypes
        )

        # ----------------------- CREATE TOPOLOGY -----------------------
        # Create the topology
        self.topology = self._createTopology(self._atomDict, self._seqDict)

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

        # Write to PDB
        from openmm.app import PDBFile

        with open("output.pdb", "w") as f:
            PDBFile.writeFile(self.topology, self.positions, f)

        exit()

        # ----------------------- CREATE XML FILES -----------------------
        self.forceFields = []
        self.implictForceFields = []
        for keyFile, atomTypes, bioTypes, forces, scalars in zip(
            key, self._atomTypes, self._bioTypes, self._forces, self._scalars
        ):
            xmlFile, implicitXmlFile = self._createXmlFile(
                keyFile, atomTypes, bioTypes, forces, scalars, self._atomData
            )

            # Reset the file pointers
            xmlFile.seek(0)
            implicitXmlFile.seek(0)

            # Store the XML files
            self.forceFields.append(xmlFile)
            self.implictForceFields.append(implicitXmlFile)

    # ------------------------------------------------------------------------------------------ #
    #                                      HELPER FUNCTIONS                                      #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _findNeighbours(
        atomData: Dict,
        index: Optional[int] = None,
        distance: Optional[int] = None,
        exact: bool = False,
        exclusionList: Optional[List[int]] = None,
    ) -> Union[List[int], List[List[int]]]:
        """
        Find atoms within or exactly at a specified distance from a given atom,
        or all connected groups (molecules) if no index is specified.

        Parameters
        ----------
        atomData : dict
            The atom data dictionary.
        index : int, optional
            The index of the atom to find neighbours for. If None, returns multiple molecules.
        distance : int, optional
            The distance to search for neighbours in number of bonds.
            If None, all neighbours are returned, including the atom itself.
        exact : bool, default=False
            If True, only atoms exactly at the specified distance are returned.
            If False, atoms within the distance (inclusive) are returned.
        exclusionList : list[int], optional
            A list of atom indices to exclude from the search.

        Returns
        -------
        list[int] or list[list[int]]
            If `index` is provided, returns a flat list of neighbors.
            If `index` is None, returns a list of lists (multiple groups of neighbors).
        """
        exclusionList = exclusionList or []
        visited = set()

        def bfs(start: int) -> List[int]:
            local_neighbours = set()
            queue = [(start, 0)]

            while queue:
                atom, dist = queue.pop(0)

                if atom in exclusionList or atom in visited:
                    continue

                visited.add(atom)

                if (
                    distance is None
                    or (exact and dist == distance)
                    or (not exact and dist <= distance)
                ):
                    local_neighbours.add(atom)

                if distance is None or dist < distance:
                    for neighbor in atomData[atom]["bonds"]:
                        if neighbor not in visited and neighbor not in exclusionList:
                            queue.append((neighbor, dist + 1))

            return sorted(local_neighbours)

        if index is not None:
            return bfs(index)

        molecules = []
        for atom in atomData:
            if atom not in visited and atom not in exclusionList:
                molecule = bfs(atom)
                if molecule:
                    molecules.append(molecule)

        return molecules

    @staticmethod
    def _createTopology(atomDict: Dict, seqDict: Optional[Dict] = None) -> top.Topology:
        """Create the topology from the molecules and sequence data.

        Parameters
        ----------
        atomDict : dict
            The atom data dictionary.
        seqDict : dict, optional
            The sequence data dictionary.
            If None, generic molecules will be used for all chains.

        Returns
        -------
        topology : Topology
            The created topology
        """
        topology = top.Topology()
        molecules = TinkerFiles._findNeighbours(atomDict)

        if seqDict is None:
            seqDict = {
                str(chainId): {"chainType": "GENERIC"}
                for chainId in range(len(molecules))
            }
        else:
            # Ensure seqDict has entries for all molecules
            for i in range(len(molecules) - len(seqDict)):
                seqDict[str(len(seqDict))] = {"chainType": "GENERIC"}

        for molecule, chainId in zip(molecules, seqDict):
            chain = topology.addChain(id=chainId)
            if seqDict[chainId]["chainType"] == "GENERIC":
                # Either a water molecule or a generic molecule
                if len(molecule) == 3:
                    atomCounts = {
                        atomDict[atomId]["atomicNumber"]: 0 for atomId in molecule
                    }
                    for atomId in molecule:
                        atomCounts[atomDict[atomId]["atomicNumber"]] += 1

                    if atomCounts[1] == 2 and atomCounts[8] == 1:
                        resName = "HOH"
                    else:
                        resName = "UNK"
                elif len(molecule) == 1:
                    # Ions
                    atom = atomDict[molecule[0]]
                    if atom["atomicNumber"] == 3:
                        resName = "Li+"
                    elif atom["atomicNumber"] == 11:
                        resName = "Na+"
                    elif atom["atomicNumber"] == 19:
                        resName = "K+"
                    elif atom["atomicNumber"] == 37:
                        resName = "Rb+"
                    elif atom["atomicNumber"] == 55:
                        resName = "Cs+"
                    elif atom["atomicNumber"] == 4:
                        resName = "Be2+"
                    elif atom["atomicNumber"] == 12:
                        resName = "Mg2+"
                    elif atom["atomicNumber"] == 20:
                        resName = "Ca2+"
                    elif atom["atomicNumber"] == 30:
                        resName = "Zn2+"
                    elif atom["atomicNumber"] == 9:
                        resName = "F-"
                    elif atom["atomicNumber"] == 17:
                        resName = "Cl-"
                    elif atom["atomicNumber"] == 35:
                        resName = "Br-"
                    elif atom["atomicNumber"] == 53:
                        resName = "I-"
                    else:
                        resName = "UNK"
                else:
                    resName = "UNK"
                residue = topology.addResidue(resName, chain)
                for atomId in molecule:
                    atom = atomDict[atomId]
                    topology.addAtom(
                        atom["nameShort"],
                        elem.Element.getByAtomicNumber(atom["atomicNumber"]),
                        residue,
                    )
            elif seqDict[chainId]["chainType"] == "NUCLEIC":
                raise NotImplementedError("Nucleic acids not implemented")
            elif seqDict[chainId]["chainType"] == "PEPTIDE":
                seenAtomIds = set()
                resId = 0
                for atomId in molecule:
                    if atomId in seenAtomIds:
                        continue

                    resName = seqDict[chainId]["residues"][resId]
                    atom = atomDict[atomId]

                    if resName == "H2N":
                        raise NotImplementedError("H2N Residue not implemented")
                    elif resName == "FOR":
                        # Formyl (FOR) - N-terminus
                        # O=CH-- // Bonded to N of next residue
                        obone = False
                        nbone = False
                        if atom["atomicNumber"] == 6:
                            neighbours = TinkerFiles._findNeighbours(
                                atomDict, atomId, 1, exact=True
                            )
                            for neighbour in neighbours:
                                if atomDict[neighbour]["atomicNumber"] == 7:
                                    nbone = True
                                elif atomDict[neighbour]["atomicNumber"] == 8:
                                    obone = True

                            if nbone and obone:
                                cai = atomId
                                ci = atomId + 1
                                oi = atomId + 2

                                otherAtoms = []
                                # Get the atoms attached to the carbonyl carbon atom (carbonyl oxygen + hydrogen)
                                ciNeighbours = TinkerFiles._findNeighbours(
                                    atomDict, ci, 1, exact=True, exclusionList=[cai]
                                )
                                ciNeighbours = [
                                    neighbour
                                    for neighbour in ciNeighbours
                                    if atomDict[neighbour]["atomicNumber"] != 7
                                ]
                                otherAtoms.extend(ciNeighbours)

                                # Get the atoms attached to the carbonyl oxygen atom
                                # Atoms of the residue
                                resAtoms = [ci, oi] + otherAtoms
                    elif resName == "ACE":
                        # Acetyl (ACE) - C-terminus
                        # O=C-CαH3
                        #   |
                        #  // Bonded to N of next residue
                        obone = False
                        nbone = False
                        if atom["atomicNumber"] == 6:
                            neighbours = TinkerFiles._findNeighbours(
                                atomDict, atomId, 2, exact=True
                            )
                            for neighbour in neighbours:
                                if atomDict[neighbour]["atomicNumber"] == 7:
                                    nbone = True
                                elif atomDict[neighbour]["atomicNumber"] == 8:
                                    obone = True
                        if obone and nbone:
                            cai = atomId
                            ci = atomId + 1
                            oi = atomId + 2

                            otherAtoms = []
                            # Get the hydrogens attached to the alpha carbon atom
                            caiNeighbours = TinkerFiles._findNeighbours(
                                atomDict, cai, 1, exact=True, exclusionList=[ci]
                            )
                            otherAtoms.extend(caiNeighbours)

                            resAtoms = [cai, ci, oi] + otherAtoms
                    elif resName == "COH":
                        raise NotImplementedError("COH Residue not implemented")
                    elif resName == "NH2":
                        raise NotImplementedError("NH2 Residue not implemented")
                    elif resName == "NME":
                        raise NotImplementedError("NME Residue not implemented")
                    else:
                        if atom["atomicNumber"] == 7:
                            # Amino acid backbone atoms
                            #              O                  O
                            #             ||                 ||
                            # H – N – Cα – C – N – Cα – C – N – Cα – C – OH
                            #     |     |     |     |     |     |     ||
                            #     H    Cβ     H    Cβ     H    Cβ     O
                            #           |           |           |
                            #          R1          R2          R3
                            neighbours = TinkerFiles._findNeighbours(
                                atomDict, atomId, 3, exact=True
                            )
                            for neighbour in neighbours:
                                if atomDict[neighbour]["atomicNumber"] == 8:
                                    # We have found the backbone atoms for an amino acid.
                                    # Tinker uses a strict numbering scheme for the atoms in a peptide bond.
                                    ni = atomId  # Amide nitrogen atom
                                    cai = atomId + 1  # Alpha carbon atom
                                    ci = atomId + 2  # Carbonyl carbon atom
                                    oi = atomId + 3  # Carbonyl oxygen atom

                                    # We need to find the atoms attached to the ni, cai, and ci and add them to the side chain atoms list
                                    otherAtoms = []

                                    # Get the hydrogen atom attached to the amide nitrogen atom
                                    niNeighbours = TinkerFiles._findNeighbours(
                                        atomDict, ni, 1, exact=True, exclusionList=[cai]
                                    )
                                    hydrogenAtoms = [
                                        neighbour
                                        for neighbour in niNeighbours
                                        if atomDict[neighbour]["atomicNumber"] == 1
                                    ]
                                    otherAtoms.extend(hydrogenAtoms)

                                    # Get the atoms attached to the alpha carbon atom (H + side chain atoms)
                                    caiNeighbours = TinkerFiles._findNeighbours(
                                        atomDict,
                                        cai,
                                        None,
                                        exact=True,
                                        exclusionList=[ni, ci],
                                    )
                                    caiNeighbours = [
                                        neighbour
                                        for neighbour in caiNeighbours
                                        if neighbour != cai
                                    ]
                                    otherAtoms.extend(caiNeighbours)

                                    # Get the atoms attached to the carbonyl carbon atom (Usually just the carbonyl oxygen atom, but sometimes 2 oxygens)
                                    ciNeighbours = TinkerFiles._findNeighbours(
                                        atomDict,
                                        ci,
                                        1,
                                        exact=True,
                                        exclusionList=[oi, cai],
                                    )

                                    # Remove oi from the list and any other amide nitrogen atoms
                                    ciNeighbours = [
                                        neighbour
                                        for neighbour in ciNeighbours
                                        if neighbour != oi
                                        and atomDict[neighbour]["atomicNumber"] != 7
                                    ]
                                    otherAtoms.extend(ciNeighbours)

                                    # Atoms of the residue
                                    resAtoms = [ni, cai, ci, oi] + otherAtoms

                    # Assign the correct residue name to the atoms
                    residue = topology.addResidue(resName, chain)
                    sortedAtoms = sorted(resAtoms)
                    for atom in sortedAtoms:
                        atomDict[atom]["residueName"] = resName
                        atomDict[atom]["chain"] = chain

                        # Add the atom to the residue
                        topology.addAtom(
                            atomDict[atom]["nameShort"],
                            elem.Element.getByAtomicNumber(
                                atomDict[atom]["atomicNumber"]
                            ),
                            residue,
                        )
                    seenAtomIds.update(sortedAtoms)
                    resId += 1

        # Add the bonds to the topology
        topology_atoms = list(topology.atoms())
        for atomId, atomIdDict in atomDict.items():
            for bond in atomIdDict["bonds"]:
                topology.addBond(topology_atoms[atomId], topology_atoms[bond])

        return topology

    def createSystem(
        self,
        nonbondedMethod=ff.PME,
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
        xmlFilesList = self.forceFields
        if implicitSolvent:
            xmlFilesList += self.implictForceFields

        # Reset the file pointers
        for f in xmlFilesList:
            f.seek(0)

        forcefield = ff.ForceField(*xmlFilesList)

        system = forcefield.createSystem(
            self.topology,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            rigidWater=rigidWater,
            removeCMMotion=removeCMMotion,
            hydrogenMass=hydrogenMass,
            polarization=polarization,
            mutualInducedTargetEpsilon=mutualInducedTargetEpsilon,
            *args,
            **kwargs,
        )

        return system

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
    #                                    SEQ FILE PARSING                                        #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _loadSeqFile(file: str) -> Dict[str, Dict[str, Union[List[str], str]]]:
        """
        Load the biopolymer sequence from a TINKER .seq file.

        Parameters
        ----------
        file : str
            The path to the .seq file.

        Returns
        -------
        Dict[str, Dict[str, Union[List[str], str]]]
            A dictionary where keys are chain IDs and values are dictionaries containing
            'residues': list of three-letter residue codes
            'chainType': string indicating the type ('PEPTIDE', 'NUCLEIC', or 'GENERIC')

        Notes
        -----
        This function is heavily based on the TINKER readseq.f routine.

        Tinker .seq files list residues belonging to one or more chains.
        Each line typically starts with a residue number, optionally preceded
        by a single-character chain ID. Lines starting with residue number 1
        indicate the beginning of a new chain. If no chain ID is provided
        for the first residue, a default ID (A, B, C, ...) is assigned.

        With Chain ID:
            A     1  CYS PHE GLU PRO PRO PRO ALA THR THR THR GLN THR GLY PHE ARG
            A    16  LEU ...

        Without Chain ID (Chain 'A' assumed for first block, 'B' for second, etc.):
            1  CYS PHE GLU PRO PRO PRO ALA THR THR THR GLN THR GLY PHE ARG
           16  LEU ...
            1  MET VAL ...  (Starts a new chain, 'B')
        """
        from collections import OrderedDict

        seqDict = OrderedDict()
        defaultChainIndex = 0
        currentChainId = None

        try:
            with open(file, "r") as f:
                for line_num, line in enumerate(f, 1):
                    lineSplit = line.split()
                    if not lineSplit:
                        continue

                    parsedChainId = None
                    startResidue = None
                    residuesStartIndex = -1

                    # Determine line format: Chain ID + Res Num or just Res Num
                    if (
                        lineSplit[0].isalpha()
                        and len(lineSplit[0]) == 1
                        and lineSplit[1].isdigit()
                    ):
                        parsedChainId = lineSplit[0]
                        startResidue = int(lineSplit[1])
                        residuesStartIndex = 2
                    elif lineSplit[0].isdigit():
                        startResidue = int(lineSplit[0])
                        residuesStartIndex = 1
                    else:
                        raise ValueError(
                            f"Line {line_num}: Does not have the expected format "
                            f"(ChainID ResNum ... or ResNum ...): {line.strip()}"
                        )

                    # Force 3-letter residues, otherwise add whitespaces to the left
                    residues = [
                        residue.rjust(3) for residue in lineSplit[residuesStartIndex:]
                    ]

                    if not residues:
                        raise ValueError(
                            f"Line {line_num}: No residues found after parsing!"
                        )

                    if startResidue == 1:
                        if parsedChainId:
                            if parsedChainId in seqDict:
                                raise ValueError(
                                    f"Line {line_num}: Duplicate start definition (residue 1) "
                                    f"for chain ID '{parsedChainId}'"
                                )
                            currentChainId = parsedChainId
                        else:
                            currentChainId = chr(65 + (defaultChainIndex % 26))
                            defaultChainIndex += 1

                        seqDict[currentChainId] = {
                            "residues": [],
                            "chainType": "GENERIC",
                        }

                    else:
                        if currentChainId is None:
                            raise ValueError(
                                f"Line {line_num}: Continuation line encountered before "
                                f"any chain was started (residue 1)."
                            )

                        if (
                            parsedChainId is not None
                            and parsedChainId != currentChainId
                        ):
                            raise ValueError(
                                f"Line {line_num}: Chain ID '{parsedChainId}' provided for residue "
                                f"{startResidue}, but expected continuation of chain '{currentChainId}'"
                            )

                        expectedResidueNum = (
                            len(seqDict[currentChainId]["residues"]) + 1
                        )
                        if startResidue != expectedResidueNum:
                            raise ValueError(
                                f"Line {line_num}: Residue numbering inconsistency for chain "
                                f"'{currentChainId}'. Expected residue {expectedResidueNum}, "
                                f"but line starts with {startResidue}."
                            )

                    # Add residues to the current chain
                    seqDict[currentChainId]["residues"].extend(residues)

        except FileNotFoundError:
            raise FileNotFoundError(f"File {file} not found")
        except IOError as e:
            raise IOError(f"Error reading file {file}: {e}")
        except Exception as e:
            raise ValueError(f"Error parsing {file}: {e}")

        # Determine chain types
        for chainId, chainData in seqDict.items():
            residues = chainData["residues"]
            aminoAcidIntersect = set(residues) & set(TinkerFiles._AMINO_ACID_LIST)
            if aminoAcidIntersect:
                chainData["chainType"] = "PEPTIDE"
                continue

            nucleotideIntersect = set(residues) & set(TinkerFiles._NUCLEOTIDE_LIST)
            if nucleotideIntersect:
                chainData["chainType"] = "NUCLEIC"
                continue

            chainData["chainType"] = "GENERIC"

        return seqDict

    # ------------------------------------------------------------------------------------------ #
    #                                   XYZ FILE PARSING                                         #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _parseAndStoreXyzLine(
        line: str, lineNum: int, xyzDict: Dict[int, Dict]
    ) -> None:
        """
        Parse a line of an .xyz file and store the data in a dictionary.

        Parameters
        ----------
        line : str
            The line containing atom data.
        lineNum : int
            The line number in the file, for error reporting.
        xyzDict : Dict[int, Dict]
            The dictionary to store parsed atom data.

        Raises
        ------
        ValueError
            If the line format is invalid.
        """
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(
                f"Line {lineNum}: Each line in the TINKER .xyz file must have at least 6 fields"
            )

        try:
            index = int(fields[0]) - 1
            symbol = str(fields[1])
            # Convert from Angstroms to nanometers
            x = float(fields[2]) * 0.1
            y = float(fields[3]) * 0.1
            z = float(fields[4]) * 0.1
            position = Vec3(x, y, z)
            atomType = str(fields[5])
            bonds = [int(bond) - 1 for bond in fields[6:]]

            xyzDict[index] = {
                "symbol": symbol,
                "positions": position,
                "atomType": atomType,
                "bonds": bonds,
            }
        except ValueError as e:
            raise ValueError(f"Line {lineNum}: Error parsing atom data: {e}")

    @staticmethod
    def _loadXyzFile(file: str) -> Tuple[Dict[int, Dict], List[Vec3], List[Vec3]]:
        """
        Load a TINKER .xyz file.

        Parameters
        ----------
        file : str
            The path to the .xyz file to load.

        Returns
        -------
        Tuple[Dict[int, Dict], List[Vec3], List[Vec3]]
            A tuple containing:
            - Dictionary with atom data (key: atom index, value: atom properties)
            - Box vectors (if periodic, otherwise None)
            - List of atomic positions

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        ValueError
            If there is an error parsing the file.
        """
        xyzDict = {}
        boxVectors = None

        try:
            with open(file, "r") as f:
                # Read number of atoms
                firstLine = f.readline().strip()
                if not firstLine:
                    raise ValueError("Empty .xyz file")

                try:
                    nAtoms = int(firstLine.split()[0])
                except (ValueError, IndexError):
                    raise ValueError("First line must contain the number of atoms")

                linesLeft = nAtoms
                lineNum = 1

                # Read the second line
                secondLine = f.readline().strip()
                lineNum += 1
                if not secondLine:
                    raise ValueError("Missing atom coordinates in .xyz file")

                secondLineSplit = secondLine.split()

                # Check for box information (second line contains box dimensions if it has 6 fields and doesn't start with "1")
                if len(secondLineSplit) == 6 and secondLineSplit[0] != "1":
                    try:
                        # Read box unit vectors and angles from the second line
                        # Convert lengths from Angstroms to nm and angles from degrees to radians
                        box = [
                            float(val) * (0.1 if i < 3 else math.pi / 180.0)
                            for i, val in enumerate(secondLineSplit)
                        ]
                        boxVectors = computePeriodicBoxVectors(*box)
                    except ValueError as e:
                        raise ValueError(
                            f"Line {lineNum}: Error parsing box vectors: {e}"
                        )
                else:
                    # No box information, so treat the second line as atom positions
                    TinkerFiles._parseAndStoreXyzLine(secondLine, lineNum, xyzDict)
                    linesLeft -= 1

                # Process the remaining atom lines
                for i in range(linesLeft):
                    lineNum += 1
                    atomLine = f.readline().strip()
                    if not atomLine:
                        raise ValueError(
                            f"Expected {nAtoms} atoms but found only {i + (1 if len(secondLineSplit) != 6 or secondLineSplit[0] == '1' else 0)}"
                        )
                    TinkerFiles._parseAndStoreXyzLine(atomLine, lineNum, xyzDict)

                # Check if we have the expected number of atoms
                if len(xyzDict) != nAtoms:
                    raise ValueError(
                        f"Expected {nAtoms} atoms but found {len(xyzDict)}"
                    )

            # Store the positions
            positions = [xyzDict[i]["positions"] for i in range(nAtoms)]
            return xyzDict, boxVectors, positions

        except FileNotFoundError:
            raise FileNotFoundError(f"File {file} not found")
        except IOError as e:
            raise IOError(f"Error reading file {file}: {e}")
        except Exception as e:
            raise ValueError(f"Error parsing {file}: {e}")

    # ------------------------------------------------------------------------------------------ #
    #                                   KEY FILE PARSING                                         #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _addAtomType(
        atomTypeDict: Dict[str, Dict],
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
        atomTypeDict : Dict[str, Dict]
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
        if atomType in atomTypeDict:
            # Validate against existing atom type data
            stored = atomTypeDict[atomType]
            mismatches = {
                "atomClass": atomClass,
                "nameShort": nameShort,
                "nameLong": nameLong,
                "element": element,
                "atomicNumber": atomicNumber,
                "mass": mass,
                "valence": valence,
            }

            for key, new_value in mismatches.items():
                if stored[key] != new_value:
                    raise ValueError(
                        f"Inconsistent {key} for atom type '{atomType}': "
                        f"expected '{stored[key]}', got '{new_value}'."
                    )

        # Add new atom type to the dictionary
        atomTypeDict[atomType] = {
            "atomClass": atomClass,
            "nameShort": nameShort,
            "nameLong": nameLong,
            "element": element,
            "atomicNumber": atomicNumber,
            "mass": mass,
            "valence": valence,
        }

    @staticmethod
    def _addBioType(
        bioTypeDict: Dict[str, Dict],
        bioType: str,
        nameShort: str,
        nameLong: str,
        atomType: str,
        element: str,
    ) -> None:
        """
        Add and validate biotype information.

        Parameters
        ----------
        bioTypeDict : Dict[str, Dict]
            The dictionary of biotype data.
        bioType : str
            The bio type.
        nameShort : str
            The short name of the biotype.
        nameLong : str
            The long name of the biotype.
        atomType : str
            The atom type counterpart of the biotype.
        element : str
            The element of the biotype.

        Raises
        ------
        ValueError
            If there is an inconsistency with existing biotype data.
        """
        if bioType in bioTypeDict:
            # Validate against existing atom type data
            stored = bioTypeDict[bioType]
            mismatches = {
                "nameShort": nameShort,
                "nameLong": nameLong,
                "atomType": atomType,
                "element": element,
            }

            for key, new_value in mismatches.items():
                if stored[key] != new_value:
                    raise ValueError(
                        f"Inconsistent {key} for biotype '{bioType}': "
                        f"expected '{stored[key]}', got '{new_value}'."
                    )

        # Add new biotype to the dictionary
        bioTypeDict[bioType] = {
            "nameShort": nameShort,
            "nameLong": nameLong,
            "atomType": atomType,
            "element": element,
        }

    @staticmethod
    def _loadKeyFile(
        keyFile: str,
    ) -> Tuple[Dict[str, Dict], Dict[str, Dict], Dict[str, List], Dict[str, str]]:
        """
        Load a TINKER .key or .prm file.

        Parameters
        ----------
        keyFile : str
            The path to the .key or .prm file.

        Returns
        -------
        Tuple[Dict[str, Dict], Dict[str, Dict], Dict[str, List], Dict[str, str]]
            A tuple containing:
            - Atom types dictionary
            - Bio types dictionary
            - Forces dictionary
            - Scalars dictionary

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        ValueError
            If there is an error parsing the file.
        """
        atomTypesDict = dict()
        bioTypesDict = dict()
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
                    # Biotype information
                    if len(fields) != 5:
                        raise ValueError(
                            f"Invalid biotype line: Expected 5 fields, got {len(fields)}. Fields: {fields}"
                        )
                    bioType, nameShort, nameLong, atomType = fields[1:]

                    # Look up element from atom type definition
                    if atomType in atomTypesDict:
                        element = atomTypesDict[atomType]["element"]
                    else:
                        # If atom type not found, assume it will be defined later
                        element = None

                    nameLong = re.sub(r"\s+", " ", nameLong.strip())
                    TinkerFiles._addBioType(
                        bioTypesDict,
                        bioType,
                        nameShort,
                        nameLong,
                        atomType,
                        element,
                    )
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

            # Update biotypes with elements from atom types if they weren't available initially
            for bioType, bioData in bioTypesDict.items():
                if bioData["element"] is None and bioData["atomType"] in atomTypesDict:
                    bioData["element"] = atomTypesDict[bioData["atomType"]]["element"]

            return atomTypesDict, bioTypesDict, forcesDict, scalarsDict

        except FileNotFoundError:
            raise FileNotFoundError(f"File {keyFile} not found")
        except IOError as e:
            raise IOError(f"Error reading file {keyFile}: {e}")
        except Exception as e:
            raise ValueError(f"Error parsing {keyFile}: {e}")

    @staticmethod
    def __addMultipole(
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
    def __addTorTor(
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

    @staticmethod
    def _combineXyzAndKeyData(xyzDict: Dict, atomTypes: List, bioTypes: List) -> Dict:
        """
        Combine the data from the .xyz and .key files into one dictionary.

        Parameters
        ----------
        xyzDict : dict
            The data from the .xyz file.
        atomTypes : dict
            The atom types from the .key file.
        bioTypes : dict
            The bio types from the .key file.

        Returns
        -------
        atomDataDict : dict
            The combined data.

        Notes
        -----
        atomDataDict[index]                          index = 0, 1, 2, ...
        atomDataDict[index]['symbol']                atom symbol
        atomDataDict[index]['positions']             atom position
        atomDataDict[index]['bonds']                 list of bonded atom indices
        atomDataDict[index]['residue']               residue name

        # if atomType is present in the .key file
        atomDataDict[index]['atomType']              atom type
        atomDataDict[index]['atomClass']             atom class
        atomDataDict[index]['nameShort']             short name
        atomDataDict[index]['nameLong']              long name
        atomDataDict[index]['element']               element
        atomDataDict[index]['atomicNumber']          atomic number
        atomDataDict[index]['mass']                  mass
        atomDataDict[index]['valence']               valence

        # if bioType is present in the .key file
        atomDataDict[index]['bioType']               bio type
        atomDataDict[index]['nameShort']             short name
        atomDataDict[index]['nameLong']              long name
        atomDataDict[index]['atomType']              atom type
        atomDataDict[index]['element']               element
        """
        atomDataDict = dict()
        for atomIndex in xyzDict:
            atomDataDict[atomIndex] = xyzDict[atomIndex]

            if "atomType" in atomDataDict[atomIndex]:
                # Add all the atom type data to the atom data
                for atomTypeDict in atomTypes:
                    if atomDataDict[atomIndex]["atomType"] in atomTypeDict:
                        atomDataDict[atomIndex].update(
                            atomTypeDict[atomDataDict[atomIndex]["atomType"]]
                        )
                        break
            elif "bioType" in atomDataDict[atomIndex]:
                # Add all the biotype data to the atom data
                for bioTypeDict in bioTypes:
                    if atomDataDict[atomIndex]["bioType"] in bioTypeDict:
                        atomDataDict[atomIndex].update(
                            bioTypeDict[atomDataDict[atomIndex]["bioType"]]
                        )
                        break

        return atomDataDict

    # ---------------------------------------------------------------------------------------- #
    #                                      WRITE XML FILE                                      #
    # ---------------------------------------------------------------------------------------- #
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
            atomTypeElement.attrib["mass"] = atomTypes[atomType]["mass"]

    @staticmethod
    def _writeResidues(
        root: etree.Element,
        atomData: Dict[int, Dict[str, Any]],
        topology: top.Topology,
        residuesSet: Set[str],
    ) -> None:
        """
        Write the residues to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        atomData : dict
            The atom data dictionary.
        topology : openmm.app.topology.Topology
            The topology object.
        """
        residuesElement = etree.SubElement(root, "Residues")
        residuesSeen = set()
        for residue in topology.residues():
            if residue.name in residuesSet and residue.name not in residuesSeen:
                residueElement = etree.SubElement(residuesElement, "Residue")
                residueElement.attrib["name"] = residue.name

                if residue.name == "HOH":
                    # AMOEBA water is flexible
                    residueElement.attrib["rigidWater"] = "false"

                for i, atom in enumerate(residue.atoms()):
                    atomElement = etree.SubElement(residueElement, "Atom")
                    atomElement.attrib["name"] = atomData[atom.index][
                        "nameShort"
                    ] + str(i)
                    atomElement.attrib["type"] = atomData[atom.index]["atomType"]

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
            for ii in range(0, 3):
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
        bioTypes : dict
            The bio types dictionary.
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

    def _createXmlFile(self, filename, atomTypes, bioTypes, forces, scalars, atomData):
        """
        Create the XML file.

        Parameters
        ----------
        filename : str
            The name of .key file.
        atomTypes : dict
            The atom types dictionary.
        bioTypes : dict
            The bio types dictionary.
        forces : dict
            The forces dictionary.
        scalars : dict
            The scalars dictionary.
        atomData : dict
            The atom data dictionary.

        Returns
        -------
        tinkXmlFile, gkXmlFile : tuple(file, file)
            The XML files.
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
        TinkerFiles._writeAtomTypes(root, atomTypes)

        # Extract residues present in this file
        residuesSet = {
            atomInfo["residue"]
            for atomInfo in atomData.values()
            if atomInfo["atomType"] in atomTypes
        }

        # Write the Residues
        TinkerFiles._writeResidues(root, atomData, self.topology, residuesSet)

        # Write the Forces
        # Bonded forces
        TinkerFiles._writeAmoebaBondForce(root, forces, scalars)
        TinkerFiles._writeAmoebaAngleForce(root, forces, scalars)
        TinkerFiles._writeAmoebaOutOfPlaneBendForce(root, forces, scalars)
        TinkerFiles._writeAmoebaTorsionForce(root, forces, scalars)
        TinkerFiles._writeAmoebaPiTorsionForce(root, forces)
        TinkerFiles._writeAmoebaStretchTorsionForce(root, forces)
        TinkerFiles._writeAmoebaAngleTorsionForce(root, forces)
        TinkerFiles._writeAmoebaStretchBendForce(root, forces)
        TinkerFiles._writeAmoebaTorsionTorsionForce(root, forces)

        # Non-bonded forces
        TinkerFiles._writeAmoebaVdwForce(root, forces, scalars)
        TinkerFiles._writeAmoebaMultipoleForce(root, forces, scalars)
        TinkerFiles._writeAmoebaUreyBradleyForce(root, forces)

        # Write the Generalized Kirkwood Force
        TinkerFiles._writeAmoebaGeneralizedKirkwoodForce(gkRoot, forces, atomTypes)
        TinkerFiles._writeAmoebaWcaDispersionForce(gkRoot, forces, scalars)

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

        if self._writeXmlFiles:
            if scalars["forcefield"] == "-2.55":
                scalars["forcefield"] = filename.rsplit(".", 1)[0]

            tinkerXmlFileName = scalars["forcefield"] + ".xml"
            with open(tinkerXmlFileName, "w") as f:
                f.write(tinkerXmlFile.getvalue())

            gkTinkerXmlFileName = scalars["forcefield"] + "_gk.xml"
            with open(gkTinkerXmlFileName, "w") as f:
                f.write(gkTinkerXmlFile.getvalue())

        return tinkerXmlFile, gkTinkerXmlFile
