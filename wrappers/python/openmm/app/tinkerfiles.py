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

__author__ = "Joao Morado"

import datetime
import io
import math
import os
import re
import shlex
import xml.etree.ElementTree as etree
from collections import OrderedDict, deque
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np

from openmm.app.internal.unitcell import computePeriodicBoxVectors
from openmm.unit import nanometers
from openmm.vec3 import Vec3

from . import element as elem
from . import forcefield as ff
from . import topology as top


class TinkerFiles:
    """TinkerFiles parses Tinker files (.xyz, .prm, .key, .seq), constructs a Topology, and (optionally) an OpenMM System from it."""

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

    _AMINO_ACID_LIST = [
        "GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "CYS", "CYX", "CYD", "PRO", 
        "PHE", "TYR", "TYD", "TRP", "HIS", "HID", "HIE", "ASP", "ASH", "ASN", "GLU", 
        "GLH", "GLN", "MET", "LYS", "LYD", "ARG", "ORN", "AIB", "PCA", "H2N", "FOR", 
        "ACE", "COH", "NH2", "NME"
        ] 

    _NUCLEOTIDE_LIST = [
        "  A", "  G", "  C", "  U", " DA", " DG", " DC", " DT", " MP", " DP", " TP"
    ] 
    
    def __init__(
        self,
        xyz: str,
        key: str,
        seq: Optional[str] = None,
        periodicBoxVectors: Optional[Tuple[Vec3, Vec3, Vec3]] = None,
        unitCellDimensions: Optional[Vec3] = None,
        writeXmlFiles: bool = False,
        writePDB: bool = False,
    ):
        """
        Load and parse Tinker files to create a Topology.

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
        writePDB : bool, optional, default=False
            If True, the PDB file is written to disk.
        """
        # Internal variables to store the data from the files
        self._xyzDict = None
        self._seqDict = None
        self._atomDict = None
        self.topology = None
        self._scalars, self._forces, self._atomTypes, self._bioTypes = [], [], [], []

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
            print(
                "WARNING: No sequence (.seq) file provided. Generic molecules will be used for all chains."
            )
            self._seqDict = None

        # Combine the data from the .xyz and .key files
        self._atomDict = TinkerFiles._combineXyzAndKeyData(
            self._xyzDict, self._atomTypes, self._bioTypes
        )

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

        # Write to PDB # TODO: Remove this after testing
        if writePDB:
            from openmm.app import PDBFile
            with open("output.pdb", "w") as f:
                PDBFile.writeFile(self.topology, self.positions, f)

        if writeXmlFiles:
            self.writeXmlFiles(key)

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
        raise NotImplementedError("createSystem not implemented")

    def writeXmlFiles(self, key):
        """
        Write the Tinker XML files.
        """
        for keyFile, atomTypes, bioTypes, forces, scalars in zip(
            key, self._atomTypes, self._bioTypes, self._forces, self._scalars
        ):
            xmlFile, implicitXmlFile = TinkerXmlWriter.writeXmlFiles(
                self.topology,
                keyFile,
                atomTypes,
                bioTypes,
                forces,
                scalars,
                self._atomDict,
            )

            # Reset the file pointers
            xmlFile.seek(0)
            implicitXmlFile.seek(0)

    # ------------------------------------------------------------------------------------------ #
    #                                      TOPOLOGY FUNCTIONS                                    #
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
        Find neighboring atoms using breadth-first search.

        Parameters
        ----------
        atomData : Dict
            Dictionary containing atom data including bond information.
        index : Optional[int], default=None
            Starting atom index. If None, returns neighbors for all atoms.
        distance : Optional[int], default=None
            Maximum distance to search (in number of bonds). If None, searches entire connected graph.
        exact : bool, default=False
            If True, only returns atoms at exactly the specified distance.
            If False, returns atoms up to and including the specified distance.
        exclusionList : Optional[List[int]], default=None
            List of atom indices to exclude from the search.

        Returns
        -------
        Union[List[int], List[List[int]]]
            If index is provided, returns list of neighbor indices.
            If index is None, returns list of lists containing neighbor indices for each atom.
        """
        exclusionList = set(exclusionList or [])
        visited = set()

        def bfs(start: int) -> List[int]:
            if start in exclusionList:
                return []

            local_neighbours = set()
            queue = deque([(start, 0)])
            visited.add(start)

            while queue:
                atom, dist = queue.popleft()

                if (
                    distance is None
                    or (exact and dist == distance)
                    or (not exact and dist <= distance)
                ):
                    local_neighbours.add(atom)

                if distance is None or dist < distance:
                    for neighbor in atomData[atom]["bonds"]:
                        if neighbor not in visited and neighbor not in exclusionList:
                            visited.add(neighbor)
                            queue.append((neighbor, dist + 1))

            return sorted(local_neighbours)

        if index is not None:
            return bfs(index)

        return [
            bfs(atom)
            for atom in atomData
            if atom not in visited and atom not in exclusionList
        ]

    @staticmethod
    def _createTopology(atomDict: Dict, seqDict: Optional[Dict] = None) -> top.Topology:
        """
        Create the topology from the molecules and sequence data.

        Parameters
        ----------
        atomDict : dict
            The atom data dictionary.
        seqDict : dict, optional
            The sequence data dictionary.
            If None, generic molecules will be used for all chains.

        Notes
        -----
        Much of the logic here attempts to reproduce the Fortran code in xyzpdb.f
        which can be found in the Tinker source code.

        Returns
        -------
        topology : Topology
            The created topology
        """
        topology = top.Topology()
        molecules = TinkerFiles._findNeighbours(atomDict)

        # Initialize sequence dictionary if not provided
        if seqDict is None:
            seqDict = {
                str(chainId): {"chainType": "GENERIC"}
                for chainId in range(len(molecules))
            }
        else:
            # Ensure seqDict has entries for all molecules
            for i in range(len(molecules) - len(seqDict)):
                seqDict[str(len(seqDict))] = {"chainType": "GENERIC"}

        # Process each molecule
        for molecule, chainId in zip(molecules, seqDict):
            chain = topology.addChain(id=chainId)
            chainType = seqDict[chainId]["chainType"]

            if chainType == "GENERIC":
                TinkerFiles._processGenericMolecule(topology, molecule, atomDict, chain)
            elif chainType == "NUCLEIC":
                TinkerFiles._processNucleicAcidChain(
                    topology, molecule, atomDict, chain, seqDict[chainId]
                )
            elif chainType == "PEPTIDE":
                TinkerFiles._processPeptideChain(
                    topology, molecule, atomDict, chain, seqDict[chainId]
                )

        # Add the bonds to the topology
        topology_atoms = list(topology.atoms())
        for atomId, atomIdDict in atomDict.items():
            for bond in atomIdDict["bonds"]:
                topology.addBond(topology_atoms[atomId], topology_atoms[bond])

        return topology

    @staticmethod
    def _processGenericMolecule(
        topology: top.Topology, molecule: List[int], atomDict: Dict, chain: top.Chain
    ) -> None:
        """
        Process a generic molecule (water, ion, or unknown).

        Parameters
        ----------
        topology : Topology
            The topology to add the molecule to.
        molecule : List[int]
            List of atom indices in the molecule.
        atomDict : Dict
            Dictionary containing atom data.
        chain : Chain
            The chain to add the molecule to.
        """
        if len(molecule) == 3:
            # Check if it's water
            atomCounts = {atomDict[atomId]["atomicNumber"]: 0 for atomId in molecule}
            for atomId in molecule:
                atomCounts[atomDict[atomId]["atomicNumber"]] += 1
            resName = "HOH" if atomCounts[1] == 2 and atomCounts[8] == 1 else "UNK"
        elif len(molecule) == 1:
            # Check if it's an ion
            atom = atomDict[molecule[0]]
            resName = TinkerFiles._getIonResidueName(atom["atomicNumber"])
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
            The residue name for the ion.
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
    def _processNucleicAcidChain(
        topology: top.Topology,
        molecule: List[int],
        atomDict: Dict,
        chain: top.Chain,
        chainData: Dict,
    ) -> None:
        """
        Process a nucleic acid chain.

        Parameters
        ----------
        topology : Topology
            The topology to add the chain to.
        molecule : List[int]
            List of atom indices in the molecule.
        atomDict : Dict
            Dictionary containing atom data.
        chain : Chain
            The chain to add the residues to.
        chainData : Dict
            Dictionary containing chain data including residues.
        """
        seenAtomIds = set()
        resId = 0

        for atomId in molecule:
            if atomId in seenAtomIds:
                continue
            resName = chainData["residues"][resId]

            # Validate residue
            assert resName in TinkerFiles._NUCLEOTIDE_LIST, (
                f"Residue {resName} is not a valid or supported nucleic acid"
            )

            # Get atoms for the residue
            resAtoms = TinkerFiles._getNucleicAcidResidueAtoms(
                atomId, resName, atomDict
            )
            if not resAtoms:
                continue

            # Add residue to topology
            TinkerFiles._addResidueToTopology(
                topology, resAtoms, resName, chain, atomDict, seenAtomIds
            )
            resId += 1

    @staticmethod
    def _getNucleicAcidResidueAtoms(
        atomId: int, resName: str, atomDict: Dict
    ) -> List[int]:
        """
        Get the atoms for a nucleic acid residue.

        Parameters
        ----------
        atomId : int
            The index of the first atom in the residue.
        resName : str
            The name of the residue.
        atomDict : Dict
            Dictionary containing atom data.

        Returns
        -------
        List[int]
            List of atom indices in the residue.
        """
        # Regular nucleotide check: Atom must be carbon (C) with 4 bonds
        atom = atomDict.get(atomId)
        if not atom or atom["atomicNumber"] != 6 or len(atom["bonds"]) != 4:
            return []

        # Ensure C1' atom is correct
        if not TinkerFiles._isC1Atom(atomId, atomDict):
            return []

        # Get sugar backbone atoms
        sugarBackboneAtoms = TinkerFiles._getNucleicSugarBackboneAtoms(atomId)
        o5, c5, c4, o4, c1, c3, c2 = sugarBackboneAtoms

        # Find neighbors for C2', C3', and C4'
        c2Neighbours = TinkerFiles._findNeighbours(
            atomDict, c2, None, exact=True, exclusionList=[c1, c3]
        )[1:]  # Exclude c2 itself
        c3Neighbours = TinkerFiles._findNeighbours(
            atomDict, c3, 1, exact=True, exclusionList=[c2, c4]
        )
        c4Neighbours = TinkerFiles._findNeighbours(
            atomDict, c4, 1, exact=True, exclusionList=[c3, o4, c5]
        )

        # Get hydrogen atoms
        c1NeighboursHydrogen = TinkerFiles._getHydrogenAtoms([c1], atomDict)
        c5NeighboursHydrogens = TinkerFiles._getHydrogenAtoms([c5], atomDict)
        o5NeighboursHydrogens = TinkerFiles._getHydrogenAtoms([o5], atomDict)

        # Determine phosphate atoms if H5T is not present
        phosphateAtoms = []
        if not o5NeighboursHydrogens:
            phosphateAtoms = [
                atomId - 7,
                atomId - 6,
                atomId - 5,
            ]  # No H5T, so include phosphate group

        # Get terminal hydrogen (H3T) if present
        o3 = [atom for atom in c3Neighbours if atomDict[atom]["atomicNumber"] == 8]
        o3NeighboursHydrogens = TinkerFiles._getHydrogenAtoms(o3, atomDict)

        # Get organic base atoms
        baseAtoms = TinkerFiles._findNeighbours(
            atomDict,
            c1,
            None,
            exact=False,
            exclusionList=[*c1NeighboursHydrogen, o4, c2],
        )[1:]  # Exclude c1 itself

        # Combine all atoms into one list
        nucleicAtoms = (
            sugarBackboneAtoms
            + c2Neighbours
            + c3Neighbours
            + c4Neighbours
            + c5NeighboursHydrogens
            + c1NeighboursHydrogen
            + baseAtoms
            + o5NeighboursHydrogens
            + phosphateAtoms
            + o3NeighboursHydrogens
        )

        return nucleicAtoms

    @staticmethod
    def _isC1Atom(atomId: int, atomDict: Dict) -> bool:
        """
        Check if an atom is a C1' atom in a nucleic acid.

        Parameters
        ----------
        atomId : int
            The index of the atom to check.
        atomDict : Dict
            Dictionary containing atom data.

        Returns
        -------
        bool
            True if the atom is a C1' atom, False otherwise.
        """
        atom = atomDict[atomId]
        hasCarbonBond = False
        hasNitrogenBond = False
        hasOxygenBond = False

        if atom["atomicNumber"] != 6:
            return False

        for bond in atom["bonds"]:
            bondAtom = atomDict[bond]
            if bondAtom["atomicNumber"] == 6:
                hasCarbonBond = True
            elif bondAtom["atomicNumber"] == 7 and len(bondAtom["bonds"]) == 3:
                hasNitrogenBond = True
            elif bondAtom["atomicNumber"] == 8 and len(bondAtom["bonds"]) == 2:
                hasOxygenBond = True

        return hasCarbonBond and hasNitrogenBond and hasOxygenBond

    @staticmethod
    def _getNucleicSugarBackboneAtoms(atomId: int) -> List[int]:
        """
        Get the sugar backbone atoms for a nucleic acid residue.

        Parameters
        ----------
        atomId : int
            The index of the C1' atom.

        Returns
        -------
        List[int]
            List of sugar backbone atom indices.
        """

        # Atomic numbers of the backbone sugar atoms
        return [
            atomId - 4,  # O5'
            atomId - 3,  # C5'
            atomId - 2,  # C4'
            atomId - 1,  # O4'
            atomId,       # C1'
            atomId + 1,  # C3'
            atomId + 2,  # C2'
        ]

    @staticmethod
    def _getHydrogenAtoms(atoms: List[int], atomDict: Dict) -> List[int]:
        """
        Get all hydrogen atoms bonded to the given atoms.

        Parameters
        ----------
        atoms : List[int]
            List of atom indices to find hydrogens for.
        atomDict : Dict
            Dictionary containing atom data.

        Returns
        -------
        List[int]
            List of hydrogen atom indices.
        """
        hydrogenAtoms = []
        for atom in atoms:
            neighbours = TinkerFiles._findNeighbours(atomDict, atom, 1, exact=True)
            hydrogenAtoms.extend(
                neighbour
                for neighbour in neighbours
                if atomDict[neighbour]["atomicNumber"] == 1
            )
        return hydrogenAtoms

    @staticmethod
    def _processPeptideChain(
        topology: top.Topology,
        molecule: List[int],
        atomDict: Dict,
        chain: top.Chain,
        chainData: Dict,
    ) -> None:
        """
        Process a peptide chain.

        Parameters
        ----------
        topology : Topology
            The topology to add the chain to.
        molecule : List[int]
            List of atom indices in the molecule.
        atomDict : Dict
            Dictionary containing atom data.
        chain : Chain
            The chain to add the residues to.
        chainData : Dict
            Dictionary containing chain data including residues.
        """
        seenAtomIds = set()
        resId = 0

        for atomId in molecule:
            if atomId in seenAtomIds:
                continue

            resName = chainData["residues"][resId]

            # Validate residue
            assert resName in TinkerFiles._AMINO_ACID_LIST, (
                f"Residue {resName} is not a valid or supported amino acid"
            )

            # Get atoms for the residue
            resAtoms = TinkerFiles._getPeptideResidueAtoms(atomId, resName, atomDict)
            if not resAtoms:
                continue

            # Add residue to topology
            TinkerFiles._addResidueToTopology(
                topology, resAtoms, resName, chain, atomDict, seenAtomIds
            )
            resId += 1

    @staticmethod
    def _getPeptideResidueAtoms(atomId: int, resName: str, atomDict: Dict) -> List[int]:
        """
        Get the atoms for a peptide residue.

        Parameters
        ----------
        atomId : int
            The index of the first atom in the residue.
        resName : str
            The name of the residue.
        atomDict : Dict
            Dictionary containing atom data.

        Returns
        -------
        List[int]
            List of atom indices in the residue.
        """
        if resName == "H2N":
            # N-terminal deprotonated residue
            return []

        if resName == "FOR":
            return TinkerFiles._getFormylResidueAtoms(atomId, atomDict)
        elif resName == "ACE":
            return TinkerFiles._getAcetylResidueAtoms(atomId, atomDict)
        elif resName == "COH":
            return TinkerFiles._getCarboxylResidueAtoms(atomId, atomDict)
        elif resName == "NH2":
            return TinkerFiles._getAmideResidueAtoms(atomId, atomDict)
        elif resName == "NME":
            return TinkerFiles._getMethylAmideResidueAtoms(atomId, atomDict)
        else:
            return TinkerFiles._getStandardAminoAcidResidueAtoms(
                atomId, resName, atomDict
            )

    @staticmethod
    def _getFormylResidueAtoms(atomId: int, atomDict: Dict) -> List[int]:
        """Get atoms for a formyl residue."""
        assert atomDict[atomId]["atomicNumber"] == 6, (
            "FOR residue must start with a carbon atom"
        )
        cai = atomId
        oi = atomId + 1
        hydrogenAtoms = TinkerFiles._getHydrogenAtoms([cai], atomDict)
        return [cai, oi] + hydrogenAtoms

    @staticmethod
    def _getAcetylResidueAtoms(atomId: int, atomDict: Dict) -> List[int]:
        """Get atoms for an acetyl residue."""
        assert atomDict[atomId]["atomicNumber"] == 6, (
            "ACE residue must start with a carbon atom"
        )
        cai = atomId
        ci = atomId + 1
        oi = atomId + 2
        hydrogenAtoms = TinkerFiles._getHydrogenAtoms([cai], atomDict)
        return [cai, ci, oi] + hydrogenAtoms

    @staticmethod
    def _getCarboxylResidueAtoms(atomId: int, atomDict: Dict) -> List[int]:
        """Get atoms for a carboxyl residue."""
        assert atomDict[atomId]["atomicNumber"] == 8, (
            "COH residue must start with an oxygen atom"
        )
        oi = atomId
        hydrogenAtoms = TinkerFiles._getHydrogenAtoms([oi], atomDict)
        return [oi] + hydrogenAtoms

    @staticmethod
    def _getAmideResidueAtoms(atomId: int, atomDict: Dict) -> List[int]:
        """Get atoms for an amide residue."""
        assert atomDict[atomId]["atomicNumber"] == 7, (
            "NH2 residue must start with a nitrogen atom"
        )
        ni = atomId
        hydrogenAtoms = TinkerFiles._getHydrogenAtoms([ni], atomDict)
        return [ni] + hydrogenAtoms

    @staticmethod
    def _getMethylAmideResidueAtoms(atomId: int, atomDict: Dict) -> List[int]:
        """Get atoms for a methyl amide residue."""
        assert atomDict[atomId]["atomicNumber"] == 7, (
            "NME residue must start with a nitrogen atom"
        )
        ni = atomId
        cai = atomId + 1
        niHydrogenAtoms = TinkerFiles._getHydrogenAtoms([ni], atomDict)
        caiHydrogenAtoms = TinkerFiles._getHydrogenAtoms([cai], atomDict)
        return [ni, cai] + niHydrogenAtoms + caiHydrogenAtoms

    @staticmethod
    def _getStandardAminoAcidResidueAtoms(
        atomId: int, resName: str, atomDict: Dict
    ) -> List[int]:
        """Get atoms for a standard amino acid residue."""
        assert atomDict[atomId]["atomicNumber"] == 7, (
            "Amino acid residue must start with a nitrogen atom"
        )

        # Backbone atoms
        ni = atomId  # Amide nitrogen
        cai = atomId + 1  # Alpha carbon
        ci = atomId + 2  # Carbonyl carbon
        oi = atomId + 3  # Carbonyl oxygen

        # Get hydrogen atoms for amide nitrogen
        niHydrogenAtoms = TinkerFiles._getHydrogenAtoms([ni], atomDict)

        # Get side chain atoms
        # For CYX residues, limit side chain connectivity to 2 bonds
        # to prevent the side chain from growing past the sulfur atom
        nBondsSideChain = 2 if resName == "CYX" else None
        caiNeighbours = TinkerFiles._findNeighbours(
            atomDict, cai, nBondsSideChain, exact=False, exclusionList=[ni, ci]
        )
        caiNeighbours = [n for n in caiNeighbours if n != cai]

        # Get atoms attached to carbonyl carbon
        ciNeighbours = TinkerFiles._findNeighbours(
            atomDict, ci, 1, exact=True, exclusionList=[oi, cai]
        )
        ciNeighbours = [n for n in ciNeighbours if atomDict[n]["atomicNumber"] != 7]

        return [ni, cai, ci, oi] + niHydrogenAtoms + caiNeighbours + ciNeighbours

    @staticmethod
    def _addResidueToTopology(
        topology: top.Topology,
        resAtoms: List[int],
        resName: str,
        chain: top.Chain,
        atomDict: Dict,
        seenAtomIds: Set[int],
    ) -> None:
        """
        Add a residue to the topology.

        Parameters
        ----------
        topology : Topology
            The topology to add the residue to.
        resAtoms : List[int]
            List of atom indices in the residue.
        resName : str
            The name of the residue.
        chain : Chain
            The chain to add the residue to.
        atomDict : Dict
            Dictionary containing atom data.
        seenAtomIds : Set[int]
            Set of atom indices that have been processed.
        """
        residue = topology.addResidue(resName, chain)
        sortedAtoms = sorted(resAtoms)

        for atom in sortedAtoms:
            # Add residue info to atom dictionary
            atomDict[atom]["residueName"] = resName
            atomDict[atom]["chain"] = chain

            # Add atom to topology
            topology.addAtom(
                atomDict[atom]["nameShort"],
                elem.Element.getByAtomicNumber(atomDict[atom]["atomicNumber"]),
                residue,
            )

        seenAtomIds.update(sortedAtoms)

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
        """Load and parse a Tinker .seq file containing sequence information.

        Parameters
        ----------
        file : str
            Path to the .seq file to load.

        Returns
        -------
        Dict[str, Dict[str, Union[List[str], str]]]
            Dictionary mapping chain IDs to dictionaries containing:
            - 'chainType': Type of chain (str)
            - 'residues': List of 3-letter residue codes (List[str])

        Raises
        ------
        ValueError
            If the file format is invalid or contains duplicate chain IDs.
        """
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
                            f"Line {line_num}: Invalid format. Expected 'ChainID ResNum ...' or 'ResNum ...'. Got: {line.strip()}"
                        )

                    if startResidue < 1:
                        raise ValueError(
                            f"Line {line_num}: Invalid residue number {startResidue}. Must be positive."
                        )

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
                                    f"Line {line_num}: Duplicate chain ID '{parsedChainId}' found. Each chain must have a unique ID."
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
                                f"Line {line_num}: Continuation line encountered before any chain was started (residue 1)."
                            )

                        if parsedChainId is not None and parsedChainId != currentChainId:
                            raise ValueError(
                                f"Line {line_num}: Chain ID '{parsedChainId}' provided for residue {startResidue}, but expected continuation of chain '{currentChainId}'"
                            )

                        expectedResidueNum = len(seqDict[currentChainId]["residues"]) + 1
                        if startResidue != expectedResidueNum:
                            raise ValueError(
                                f"Line {line_num}: Residue numbering inconsistency for chain '{currentChainId}'. Expected residue {expectedResidueNum}, but line starts with {startResidue}."
                            )

                    seqDict[currentChainId]["residues"].extend(residues)

        except FileNotFoundError:
            raise FileNotFoundError(f"Sequence file not found: {file}")
        except IOError as e:
            raise IOError(f"Error reading sequence file {file}: {e}")
        except Exception as e:
            raise ValueError(f"Error parsing sequence file {file}: {e}")

        if not seqDict:
            raise ValueError(f"No valid chains found in sequence file {file}")

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
    def _loadXyzFile(file: str) -> Tuple[Dict[int, Dict], List[Vec3], List[Vec3]]:
        """
        Load and parse a Tinker .xyz file.

        Parameters
        ----------
        file : str
            Path to the .xyz file to load.

        Returns
        -------
        Tuple[Dict[int, Dict], List[Vec3], List[Vec3]]
            A tuple containing:
            - Dictionary mapping atom indices to their properties
            - List of atomic positions
            - List of periodic box vectors (if present, otherwise None)
        """
        xyzDict = {}
        boxVectors = None

        try:
            with open(file, "r") as f:
                firstLine = f.readline().strip()
                if not firstLine:
                    raise ValueError("Empty .xyz file")

                try:
                    nAtoms = int(firstLine.split()[0])
                    if nAtoms <= 0:
                        raise ValueError(f"Invalid number of atoms: {nAtoms}. Must be positive.")
                except (ValueError, IndexError):
                    raise ValueError("First line must contain a positive integer number of atoms")

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
                            raise ValueError("Box vectors contain invalid values (NaN or Inf)")
                        boxVectors = computePeriodicBoxVectors(*box)
                    except ValueError as e:
                        raise ValueError(f"Line {lineNum}: Error parsing box vectors: {e}")
                else:
                    TinkerFiles._parseAndStoreXyzLine(secondLine, lineNum, xyzDict)
                    linesLeft -= 1

                for i in range(linesLeft):
                    lineNum += 1
                    atomLine = f.readline().strip()
                    if not atomLine:
                        raise ValueError(
                            f"Expected {nAtoms} atoms but found only {i + (1 if len(secondLineSplit) != 6 or secondLineSplit[0] == '1' else 0)}"
                        )
                    TinkerFiles._parseAndStoreXyzLine(atomLine, lineNum, xyzDict)

                if len(xyzDict) != nAtoms:
                    raise ValueError(
                        f"Expected {nAtoms} atoms but found {len(xyzDict)}"
                    )

            positions = [xyzDict[i]["positions"] for i in range(nAtoms)]
            return xyzDict, boxVectors, positions

        except FileNotFoundError:
            raise FileNotFoundError(f"XYZ file not found: {file}")
        except IOError as e:
            raise IOError(f"Error reading XYZ file {file}: {e}")
        except Exception as e:
            raise ValueError(f"Error parsing XYZ file {file}: {e}")

    @staticmethod
    def _parseAndStoreXyzLine(line: str, lineNum: int, xyzDict: Dict[int, Dict]) -> None:
        """
        Parse a single line from a TINKER .xyz file and store the data in the xyzDict.

        Parameters
        ----------
        line : str
            The line to parse from the .xyz file.
        lineNum : int
            The line number in the file (for error reporting).
        xyzDict : Dict[int, Dict]
            Dictionary to store the parsed atom data, keyed by atom index.
        """
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(
                f"Line {lineNum}: Each line in the TINKER .xyz file must have at least 6 fields"
            )

        try:
            index = int(fields[0]) - 1
            if index < 0:
                raise ValueError(f"Line {lineNum}: Invalid atom index {index + 1}. Must be positive.")

            symbol = str(fields[1])
            if not symbol:
                raise ValueError(f"Line {lineNum}: Empty atom symbol")

            try:
                x = float(fields[2]) * 0.1
                y = float(fields[3]) * 0.1
                z = float(fields[4]) * 0.1
                if any(math.isnan(coord) or math.isinf(coord) for coord in (x, y, z)):
                    raise ValueError("Atom coordinates contain invalid values (NaN or Inf)")
                position = Vec3(x, y, z)
            except ValueError as e:
                raise ValueError(f"Line {lineNum}: Error parsing atom coordinates: {e}")

            atomType = str(fields[5])
            if not atomType:
                raise ValueError(f"Line {lineNum}: Empty atom type")

            try:
                bonds = [int(bond) - 1 for bond in fields[6:]]
                if any(bond < 0 for bond in bonds):
                    raise ValueError("Invalid bond index (must be positive)")
            except ValueError as e:
                raise ValueError(f"Line {lineNum}: Error parsing bond indices: {e}")

            if index in xyzDict:
                raise ValueError(f"Line {lineNum}: Duplicate atom index {index + 1}")

            xyzDict[index] = {
                "symbol": symbol,
                "positions": position,
                "atomType": atomType,
                "bonds": bonds,
            }
        except ValueError as e:
            raise ValueError(f"Line {lineNum}: {str(e)}")

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


class TinkerXmlWriter:
    """
    Write the Tinker XML file.
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
        atomDict: Dict[int, Dict[str, Any]],
        topology: top.Topology,
    ) -> None:
        """
        Write the residues to the XML file.

        Parameters
        ----------
        root : xml.etree.ElementTree.Element
            The root element of the XML tree.
        atomDict : dict
            The atom data dictionary.
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
                    atomElement.attrib["name"] = atomDict[atom.index][
                        "nameShort"
                    ] + str(i)
                    atomElement.attrib["type"] = atomDict[atom.index]["atomType"]

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

    @staticmethod
    def writeXmlFiles(
        topology, filename, atomTypes, bioTypes, forces, scalars, atomDict
    ):
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
        atomDict : dict
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
