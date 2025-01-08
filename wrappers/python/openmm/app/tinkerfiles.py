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
import math
import os
import re
import shlex
import xml.etree.ElementTree as etree
from functools import wraps
from typing import List, Union

from openmm.app.internal.unitcell import computePeriodicBoxVectors
from openmm.unit import Quantity, nanometers
from openmm.vec3 import Vec3

from . import element as elem
from . import forcefield as ff
from . import topology as top

try:
    import numpy as np
except ImportError:
    np = None


def numpy_protector(func):
    """
    Decorator to emit useful error messages if users try to request numpy
    processing if numpy is not available. Raises ImportError if numpy could not
    be found
    """

    @wraps(func)
    def wrapper(self, asNumpy=False):
        if asNumpy and np is None:
            raise ImportError("Could not import numpy. Cannot set asNumpy=True")
        return func(self, asNumpy=asNumpy)

    return wrapper


class TinkerFiles:
    """TinkerFiles parses Tinker files (.xyz, .prm, .key), constructs a Topology, and (optionally) an OpenMM System from it."""

    def __init__(
        self,
        xyz: str,
        key: str,
        periodicBoxVectors=None,
        unitCellDimensions=None,
    ):
        """
        Load exactly one .xyz file and one or more .key or .prm files.

        Parameters
        ----------
        xyz : str
            The path to the xyz file to load.
        key : str or list of str
            The path(s) to the key/prm file(s) to load.
        periodicBoxVectors : tuple of Vec3
            The vectors defining the periodic box.
            If provided, this overwrites the box information from the xyz file.
        unitCellDimensions : Vec3, optional
            The dimensions of the crystallographic unit cell.
            For non-rectangular unit cells, specify periodicBoxVectors instead.
            If provided, this overwrites the box information from the xyz file.
        """
        # Position and box information
        self.positions = None
        self.boxVectors = None
        self._numpyPositions = None
        self._numpyBoxVectors = None

        # Internal variables to store the data from the xyz file
        self._symbols = []
        self._names = []
        self._bonds = []
        self._atomTypes = []

        # Internal variables to store the data from the key file
        self._atomTypesDict = {}
        self._residueAtomTypes = None
        self._residues = None
        self._residueNames = None

        # Load the key/prm file(s)
        key = key if isinstance(key, list) else [key]
        for keyFile in key:
            self._loadKeyFile(keyFile)

        # Load the xyz file
        self._loadXyzFile(xyz)
        self.positions = Quantity(self.positions, nanometers)

        # Create the topology
        self.topology = top.Topology()
        self._createTopology()

        # Set the periodic box vectors as specified in the xyz file
        if self.boxVectors is not None:
            self.topology.setPeriodicBoxVectors(self.boxVectors)

        # If provided, set the periodic box vectors or unit cell dimensions in the topology
        # This overwrites the box information from the xyz file
        if periodicBoxVectors is not None:
            if unitCellDimensions is not None:
                raise ValueError(
                    "Specify either periodicBoxVectors or unitCellDimensions, but not both"
                )
            self.topology.setPeriodicBoxVectors(periodicBoxVectors)
        else:
            self.topology.setUnitCellDimensions(unitCellDimensions)

        # Create the XML file for the residues and the residue dictionary
        residueXmlFile = self._createResidueXmlFiles()
        residueDict = TinkerXmlFileWriter.createResidueDict(
            residueXmlFile, self._residueAtomTypes
        )

        # Create the XML files
        xmlFilesList = []
        for keyFile in key:
            atomTypes, bioTypes, forces, scalars = TinkerXmlFileWriter.processKeyFile(
                keyFile, residueDict
            )
            xmlFiles = TinkerXmlFileWriter.createXmlFiles(
                keyFile, atomTypes, bioTypes, forces, residueDict, scalars
            )
            xmlFilesList.extend(xmlFiles)

    def createSystem(
        self,
        nonbondedMethod=ff.NoCutoff,
        nonbondedCutoff=1.0 * nanometers,
        constraints=None,
        rigidWater=True,
        removeCMMotion=True,
        hydrogenMass=None,
    ):
        forcefield = ff.ForceField(
            "/home/joaomorado/tests/amoeba-poltype/SolvE0.1_V1/AMOEBA-BIO-2018.xml",
            "/home/joaomorado/tests/amoeba-poltype/SolvE0.1_V1/final.xml",
        )

        system = forcefield.createSystem(
            self.topology,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            rigidWater=rigidWater,
            removeCMMotion=removeCMMotion,
            hydrogenMass=hydrogenMass,
        )

        return system

    # ------------------------------------------------------------------------------------------ #
    #                                     POSITIONS AND BOX                                      #
    # ------------------------------------------------------------------------------------------ #
    @numpy_protector
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
                self._numpyPositions = Quantity(
                    np.array(self.positions.value_in_unit(nanometers)), nanometers
                )
            return self._numpyPositions
        return self.positions

    @numpy_protector
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
                    Quantity(
                        np.array(self.boxVectors[0].value_in_unit(nanometers)),
                        nanometers,
                    )
                )
                self._numpyBoxVectors.append(
                    Quantity(
                        np.array(self.boxVectors[1].value_in_unit(nanometers)),
                        nanometers,
                    )
                )
                self._numpyBoxVectors.append(
                    Quantity(
                        np.array(self.boxVectors[2].value_in_unit(nanometers)),
                        nanometers,
                    )
                )
            return self._numpyBoxVectors
        return self.boxVectors

    # ------------------------------------------------------------------------------------------ #
    #                                   XYZ FILE PARSING                                         #
    # ------------------------------------------------------------------------------------------ #
    def _parseXyzLine(self, line: str) -> None:
        """
        Parse a line from a TINKER .xyz file.

        Parameters
        ----------
        line : str
            The line to parse
        """
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(
                "Each line in the TINKER .xyz file must have at least 6 fields"
            )

        # Extract info
        index = int(fields[0])
        symbol = str(fields[1])
        x = float(fields[2]) * 0.1
        y = float(fields[3]) * 0.1
        z = float(fields[4]) * 0.1

        self.positions.append(Vec3(x, y, z))
        self._symbols.append(symbol)
        self._atomTypes.append(fields[5])
        self._bonds.append([int(i) - 1 for i in fields[6:]])

    def _loadXyzFile(self, file: str) -> None:
        """
        Load a TINKER .xyz file.

        Parameters
        ----------
        file : str
            The name of the .xyz file to load.
        """
        try:
            with open(file, "r") as f:
                self.positions = []
                nAtoms = int(f.readline())
                secondLine = f.readline().strip()
                secondLineSplit = secondLine.split()

                if len(secondLineSplit) == 6 and secondLineSplit[0] != "1":
                    # Read box unit vectors and angles from the second line
                    box = [
                        float(val) * (0.1 if i < 3 else math.pi / 180.0)
                        for i, val in enumerate(secondLineSplit)
                    ]
                    self.boxVectors = computePeriodicBoxVectors(*box)
                else:
                    # No box information, so treat the second line as atom positions
                    self._parseXyzLine(secondLine)
                    nAtoms -= 1

                # Read the remaining atom positions
                for _ in range(nAtoms):
                    self._parseXyzLine(f.readline().strip())
        except FileNotFoundError:
            raise IOError(f"Could not find file {file}")
        except Exception as e:
            raise ValueError(f"Error parsing {file}: {e}")

    # ------------------------------------------------------------------------------------------ #
    #                                   KEY FILE PARSING                                         #
    # ------------------------------------------------------------------------------------------ #
    def _parseAtomLine(self, fields: List[str]) -> None:
        """
        Parse an atom type line from a TINKER .key file.

        Parameters
        ----------
        fields : list of str
            The fields of the line to parse. Must contain at least 8 fields in the following order:
            [atomIndex, atomType, atomClass, nameShort, nameLong, atomicNumber, mass, valence].

        Raises
        ------
        ValueError
            If the line has fewer than 8 fields or if the atom information is inconsistent.
        """
        if len(fields) < 8:
            raise ValueError(
                f"Invalid atom line: Expected at least 8 fields, got {len(fields)}. Fields: {fields}"
            )

        # Extract atom information
        _, atomType, atomClass, nameShort, nameLong, atomicNumber, mass, valence = (
            fields
        )

        nameLong = nameLong.replace('"', "")

        if atomType in self._atomTypesDict:
            # Validate against existing atom type data
            stored = self._atomTypesDict[atomType]
            mismatches = {
                "atomClass": atomClass,
                "nameShort": nameShort,
                "nameLong": nameLong,
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
        else:
            # Add new atom type to the dictionary
            self._atomTypesDict[atomType] = {
                "atomClass": str(atomClass),
                "nameShort": str(nameShort),
                "nameLong": str(nameLong),
                "atomicNumber": int(atomicNumber),
                "mass": float(mass),
                "valence": int(valence),
            }

    def _loadKeyFile(self, file: str):
        """
        Load a TINKER .key or .prm file to extract atom type information.

        Parameters
        ----------
        file : str
            The name of the .key or .prm file to load.
        """
        try:
            with open(file, "r") as f:
                for line in f:
                    fields = re.findall(r"\"[^\"]*\"|\S+", line)
                    if len(fields) < 2:
                        continue
                    if fields[0] == "atom":
                        self._parseAtomLine(fields)
        except FileNotFoundError:
            raise IOError(f"Could not find file {file}")
        except Exception as e:
            raise ValueError(f"Error parsing {file}: {e}")

    # ------------------------------------------------------------------------------------------ #
    #                                        XML FILES                                           #
    # ------------------------------------------------------------------------------------------ #
    def _createResidueXmlFiles(self):
        """Create the residue XML files."""
        root = etree.Element("Residues")
        seenResidues = set()
        self._residueAtomTypes = dict()
        for residue in self.topology.residues():
            if residue.name in seenResidues:
                continue

            if "water" in residue.name.lower():
                resType = "water"
            elif "ion" in residue.name.lower():
                resType = "ion"
            else:
                resType = "smallMolecule"

            # Create a new residue template
            res = etree.SubElement(root, "Residue")
            res.attrib["abbreviation"] = residue.name[:5].upper()
            res.attrib["loc"] = "free"
            res.attrib["type"] = resType
            res.attrib["tinkerLookupName"] = residue.name
            res.attrib["fullName"] = residue.name

            # Add this residue to the residue dictionary
            self._residueAtomTypes[residue.name] = dict()

            # Add atoms to the residue
            for atomId, atom in enumerate(residue.atoms()):
                atom.name = f"{atom.name}{atomId}"
                at = etree.SubElement(res, "Atom")
                at.attrib["name"] = atom.name
                at.attrib["tinkerLookupName"] = atom.name
                at.attrib["bonds"] = str(len(self._bonds[atom.index]))

                # Add the atom type to the residue dictionary
                self._residueAtomTypes[residue.name][atom.name] = self._atomTypes[
                    atom.index
                ]

            # Add bonds to the residue
            for bond in residue.bonds():
                bo = etree.SubElement(res, "Bond")
                bo.attrib["from"] = bond.atom1.name
                bo.attrib["to"] = bond.atom2.name

            seenResidues.add(residue.name)

        tree = etree.ElementTree(root)
        etree.indent(root, "  ")
        filename = "residues_trial.xml"
        tree.write(filename, xml_declaration=True, encoding="utf-8", method="xml")

        return filename

    # ------------------------------------------------------------------------------------------ #
    #                                        TOPOLOGY                                            #
    # ------------------------------------------------------------------------------------------ #
    def _createTopology(self) -> top.Topology:
        """
        Build the topology from the data parsed from the .xyz and .key files.

        Returns
        -------
        openmm.app.topology.Topology
            The topology object.
        """
        # Infer the residues from the bonds
        self._getResiduesFromBonds()

        # Add chain to the topology
        for residueName, residueAtoms in zip(self._residueNames, self._residues):
            # Add chain to the topology
            chain = self.topology.addChain()
            # Add residues to the topology
            residue = self.topology.addResidue(residueName, chain)
            for atomIndex in residueAtoms:
                # Add atoms to the topology
                atomType = self._atomTypes[atomIndex]
                atomTypeData = self._atomTypesDict[atomType]
                element = elem.Element.getByAtomicNumber(atomTypeData["atomicNumber"])
                self.topology.addAtom(self._symbols[atomIndex], element, residue)

        # Add bonds to the topology
        atoms = list(self.topology.atoms())
        seenBonds = set()
        for atomIndex, bondedAtoms in enumerate(self._bonds):
            for bondedAtomIndex in bondedAtoms:
                if (bondedAtomIndex, atomIndex) not in seenBonds:
                    self.topology.addBond(atoms[atomIndex], atoms[bondedAtomIndex])
                    seenBonds.add((bondedAtomIndex, atomIndex))

        return self.topology

    def _getResiduesFromBonds(self):
        """
        Form residues by recursively traversing bonded atoms.

        Notes
        -----
        Residues are defined as connected atoms in the graph formed by the bonds between atoms.
        For example, connected amino acids in a protein are considered to be part of the same residue.

        Returns
        -------
        list of list of int
            A list of residues, where each residue is a list of atom indices.
        """

        def _getResidueName(atomIndex):
            atomType = self._atomTypes[atomIndex]
            residueName = self._atomTypesDict[atomType]["nameLong"].replace('"', "")
            nameParts = residueName.split()
            if len(nameParts) == 2:
                residueName = nameParts[0]
            else:
                residueName = " ".join(nameParts[:2])
            return residueName

        self._residues = []
        self._residueNames = []
        seen = set()

        lastResidue = None
        for atom1 in range(len(self._bonds)):
            if atom1 not in seen:
                # Start a new residue
                residue = []

                # Get the residue name
                residueName = _getResidueName(atom1)

                # Iterative way of traversing a graph-like structure
                # Recursive approach led to stack overflow
                stack = [atom1]
                while stack:
                    atom = stack.pop()
                    if atom not in seen:
                        seen.add(atom)
                        residue.append(atom)
                        stack.extend(self._bonds[atom])

                self._residues.append(sorted(residue))
                self._residueNames.append(residueName)

        return self._residues, self._residueNames


class TinkerXmlFileWriter:
    """
    TinkerXmlFileWriter is a class that provides methods for parsing Tinker .key and .prm files
    and converting them into XML files that are compatible with OpenMM.

    Notes
    -----
    The methods in this class are largely inspired by the functions in `devtools/forcefield-scripts/processTinkerForceField.py`.
    """

    RECOGNIZED_FORCES = {
        "bond": 1,
        "angle": 1,
        "anglep": 1,
        "strbnd": 1,
        "ureybrad": 1,
        "opbend": 1,
        "torsion": 1,
        "pitors": 1,
        "strtors": 1,
        "angtors": 1,
        "vdw": 1,
        "vdwpr": 1,
        "polarize": 1,
        "tortors": None,
        "multipole": None,
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

    @staticmethod
    def __addMultipole(lineIndex, allLines, forces):
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
    def __addTorTor(lineIndex, allLines, forces):
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

    RECOGNIZED_FORCES["tortors"] = __addTorTor
    RECOGNIZED_FORCES["multipole"] = __addMultipole

    # ------------------------------------------------------------------------------------------ #
    #                                  PRIVATE STATIC METHODS                                    #
    # ------------------------------------------------------------------------------------------ #
    @staticmethod
    def _getDefaultAtom():
        """
        Get the default atom dictionary.

        Returns
        -------
        dict
            The default atom dictionary.
        """
        atom = dict()
        atom["tinkerLookupName"] = "XXX"
        atom["type"] = -1
        atom["bonds"] = dict()
        return atom

    @staticmethod
    def _addBond(atomDict, atom1, atom2):
        """
        Add a bond between two atoms to the atom dictionary.

        Parameters
        ----------
        atomDict : dict
            The atom dictionary.
        atom1 : str
            The name of the first atom.
        atom2 : str
            The name of the second atom.
        """
        if atom1 not in atomDict:
            atomDict[atom1] = TinkerXmlFileWriter._getDefaultAtom()
        if atom2 not in atomDict:
            atomDict[atom2] = TinkerXmlFileWriter._getDefaultAtom()
        atomDict[atom2]["bonds"][atom1] = 1
        atomDict[atom1]["bonds"][atom2] = 1

    @staticmethod
    def _getXmlAtoms(atoms):
        """
        Get the atom information from the XML file.

        Parameters
        ----------
        atoms : list
            The list of atoms.

        Returns
        -------
        dict
            The atom information.
        """
        atomInfo = dict()
        for atom in atoms:
            name = atom.attrib["name"]
            atomInfo[name] = TinkerXmlFileWriter._getDefaultAtom()
            atomInfo[name]["tinkerLookupName"] = atom.attrib["tinkerLookupName"]
        return atomInfo

    @staticmethod
    def _getXmlBonds(bonds):
        """
        Get the bond information from the XML file.

        Parameters
        ----------
        bonds : list
            The list of bonds.

        Returns
        -------
        dict
            The bond information.
        """
        bondInfo = dict()
        for bond in bonds:
            atom1 = bond.attrib["from"]
            atom2 = bond.attrib["to"]
            if atom1 not in bondInfo:
                bondInfo[atom1] = dict()
            if atom2 not in bondInfo:
                bondInfo[atom2] = dict()
            bondInfo[atom1][atom2] = 1
            bondInfo[atom2][atom1] = 1
        return bondInfo

    @staticmethod
    def _copyBonds(bonds):
        """
        Copy the bond information to a new dictionary.

        Parameters
        ----------
        bonds : dict
            The bond information.

        Returns
        -------
        dict
            The new dictionary with the bond information.
        """
        bondCopy = dict()
        for key in bonds.keys():
            bondCopy[key] = bonds[key]
        return bondCopy

    @staticmethod
    def _copyAtom(atom):
        """
        Copy the atom information to a new dictionary.

        Parameters
        ----------
        atom : dict
            The atom information.

        Returns
        -------
        dict
            The new dictionary with the atom information.
        """
        atomCopy = dict()
        for key in atom.keys():
            if key != "bonds":
                atomCopy[key] = atom[key]
            else:
                atomCopy["bonds"] = TinkerXmlFileWriter._copyBonds(atom[key])
        return atomCopy

    @staticmethod
    def _copyProteinResidue(residue):
        """
        Copy the protein residue information to a new dictionary.

        Parameters
        ----------
        residue : dict
            The protein residue information.

        Returns
        -------
        dict
            The new dictionary with the protein residue information.
        """
        residueCopy = dict()
        residueCopy["atoms"] = dict()
        residueCopy["type"] = residue["type"]
        residueCopy["loc"] = residue["loc"]
        residueCopy["tinkerLookupName"] = residue["tinkerLookupName"]
        residueCopy["residueName"] = residue["residueName"]
        residueCopy["include"] = residue["include"]
        for atom in residue["atoms"]:
            residueCopy["atoms"][atom] = TinkerXmlFileWriter._copyAtom(
                residue["atoms"][atom]
            )
        return residueCopy

    @staticmethod
    def _buildProteinResidue(
        residueDict,
        atoms,
        bondInfo,
        abbr,
        loc,
        tinkerLookupName,
        include,
        residueName,
        resType,
    ):
        """
        Build a protein residue.

        Parameters
        ----------
        residueDict : dict
            The residue dictionary.
        atoms : dict
            The atom information.
        bondInfo : dict
            The bond information.
        abbr : str
            The abbreviation of the residue.
        loc : str
            The location of the residue.
        tinkerLookupName : str
            The Tinker lookup name of the residue.
        include : int
            The include flag.
        residueName : str
            The name of the residue.
        resType : str
            The type of the residue.
        """
        residueDict[abbr] = dict()
        residueDict[abbr]["atoms"] = atoms
        residueDict[abbr]["type"] = resType
        residueDict[abbr]["loc"] = loc
        residueDict[abbr]["tinkerLookupName"] = tinkerLookupName
        residueDict[abbr]["residueName"] = residueName
        residueDict[abbr]["include"] = include

        for atom in bondInfo:
            if atom in residueDict[abbr]["atoms"]:
                if "bonds" not in residueDict[abbr]["atoms"][atom]:
                    residueDict[abbr]["atoms"][atom]["bonds"] = dict()
                for bondedAtom in bondInfo[atom]:
                    if bondedAtom in residueDict[abbr]["atoms"]:
                        if "bonds" not in residueDict[abbr]["atoms"][bondedAtom]:
                            residueDict[abbr]["atoms"][bondedAtom]["bonds"] = dict()
                        residueDict[abbr]["atoms"][bondedAtom]["bonds"][atom] = 1
                        residueDict[abbr]["atoms"][atom]["bonds"][bondedAtom] = 1
                    else:
                        print(f"Error: bonded atom={atom} not in residue={abbr}")
            else:
                print(f"Error: bonded atom={atom} not in residue={abbr}")

        return residueDict

    @classmethod
    def createResidueDict(cls, residueXmlFileName, residueAtomTypes):
        """
        Create the residue dictionary from defined residues in the XML file.

        Parameters
        ----------
        residueXmlFileName : str
            The name of the XML file containing the residue information.
        residueAtomTypes : dict
            Dictionary containing the atom types for each residue.

        Returns
        -------
        dict
            The residue dictionary.
        """
        residueTree = etree.parse(residueXmlFileName)
        print(f"Read {residueXmlFileName}")
        root = residueTree.getroot()
        residueDict = dict()
        for residue in root.findall("Residue"):
            abbr = residue.attrib["abbreviation"]
            loc = residue.attrib["loc"]
            resType = residue.attrib["type"]
            tinkerName = residue.attrib["tinkerLookupName"]
            residueName = residue.attrib["fullName"]

            if resType not in (
                "smallMolecule",
                "protein",
                "water",
                "dna",
                "rna",
                "ion",
            ):
                raise ValueError(f"Unknown residue type: {resType}")

            atoms = TinkerXmlFileWriter._getXmlAtoms(residue.findall("Atom"))

            # Add atom types to the atoms
            for atom in atoms:
                try:
                    atoms[atom]["type"] = residueAtomTypes[residueName][atom]
                except KeyError:
                    raise ValueError(
                        f"Error: atom={atom} not in residueAtomTypes for residue={residueName}"
                    )

            bondInfo = TinkerXmlFileWriter._getXmlBonds(residue.findall("Bond"))
            if resType == "water":
                TinkerXmlFileWriter._buildProteinResidue(
                    residueDict,
                    atoms,
                    bondInfo,
                    abbr,
                    "x",
                    tinkerName,
                    1,
                    "HOH",
                    "water",
                )
            elif resType == "protein":
                TinkerXmlFileWriter._buildProteinResidue(
                    residueDict,
                    atoms,
                    bondInfo,
                    abbr,
                    "m",
                    tinkerName,
                    1,
                    residueName,
                    "protein",
                )
                cResidueName = "C" + abbr
                residueDict[cResidueName] = TinkerXmlFileWriter._copyProteinResidue(
                    residueDict[abbr]
                )
                residueDict[cResidueName]["loc"] = "c"
                if residueDict[abbr]["tinkerLookupName"].find("(") > -1:
                    begin = residueDict[abbr]["tinkerLookupName"].find("(")
                    end = residueDict[abbr]["tinkerLookupName"].find(")") + 1
                    sub = residueDict[abbr]["tinkerLookupName"][begin:end]
                    if sub == "(HD)" or sub == "(HE)":
                        residueDict[cResidueName]["tinkerLookupName"] = (
                            "C-Terminal " + "HIS " + sub
                        )
                    else:
                        residueDict[cResidueName]["tinkerLookupName"] = (
                            "C-Terminal " + abbr + " " + sub
                        )
                    print(
                        f"tinkerLookupName {abbr} {residueDict[cResidueName]['tinkerLookupName']}"
                    )
                else:
                    residueDict[cResidueName]["tinkerLookupName"] = "C-Terminal " + abbr
                residueDict[cResidueName]["atoms"]["OXT"] = (
                    TinkerXmlFileWriter._copyAtom(residueDict[abbr]["atoms"]["O"])
                )
                residueDict[cResidueName]["atoms"]["OXT"]["tinkerLookupName"] = "OXT"
                residueDict[cResidueName]["atoms"]["O"]["tinkerLookupName"] = "OXT"
                residueDict[cResidueName]["parent"] = residueDict[abbr]
                nResidueName = "N" + abbr
                residueDict[nResidueName] = TinkerXmlFileWriter._copyProteinResidue(
                    residueDict[abbr]
                )
                residueDict[nResidueName]["loc"] = "n"
                residueDict[nResidueName]["tinkerLookupName"] = "N-Terminal " + abbr
                residueDict[nResidueName]["parent"] = residueDict[abbr]
                if abbr == "PRO":
                    residueDict[nResidueName]["atoms"][
                        "H2"
                    ] = TinkerXmlFileWriter._getDefaultAtom()
                    residueDict[nResidueName]["atoms"][
                        "H3"
                    ] = TinkerXmlFileWriter._getDefaultAtom()
                    residueDict[nResidueName]["atoms"]["H2"]["tinkerLookupName"] = "HN"
                    residueDict[nResidueName]["atoms"]["H3"]["tinkerLookupName"] = "HN"
                    TinkerXmlFileWriter._addBond(
                        residueDict[nResidueName]["atoms"], "H2", "N"
                    )
                    TinkerXmlFileWriter._addBond(
                        residueDict[nResidueName]["atoms"], "H3", "N"
                    )
                else:
                    residueDict[nResidueName]["atoms"]["H2"] = (
                        TinkerXmlFileWriter._copyAtom(residueDict[abbr]["atoms"]["H"])
                    )
                    residueDict[nResidueName]["atoms"]["H3"] = (
                        TinkerXmlFileWriter._copyAtom(residueDict[abbr]["atoms"]["H"])
                    )
            elif resType == "dna" or resType == "rna":
                TinkerXmlFileWriter._buildProteinResidue(
                    residueDict,
                    atoms,
                    bondInfo,
                    abbr,
                    loc,
                    tinkerName,
                    1,
                    residueName,
                    resType,
                )
            else:
                TinkerXmlFileWriter._buildProteinResidue(
                    residueDict,
                    atoms,
                    bondInfo,
                    abbr,
                    loc,
                    tinkerName,
                    1,
                    residueName,
                    resType,
                )

        print(residueDict)
        return residueDict

    @staticmethod
    def _setBioTypes(bioTypes, residueDict):
        """
        Set the biotypes for all atoms.

        Parameters
        ----------
        bioTypes : dict
            The biotype information.
        residueDict : dict
            The residue dictionary.

        Returns
        -------
        int
            The return code.
        """
        for resname, res in residueDict.items():
            for atom in res["atoms"]:
                atomLookup = res["atoms"][atom]["tinkerLookupName"]
                resLookup = []
                if res["type"] == "dna":
                    if res["loc"] in ("5", "N"):
                        resLookup.append("5'-Hydroxyl DNA")
                    if res["loc"] in ("3", "N"):
                        resLookup.append("3'-Hydroxyl DNA")
                    resLookup.append("Phosphodiester DNA")
                if res["type"] == "rna":
                    if res["loc"] in ("5", "N"):
                        resLookup.append("5'-Hydroxyl RNA")
                    if res["loc"] in ("3", "N"):
                        resLookup.append("3'-Hydroxyl RNA")
                    resLookup.append("Phosphodiester RNA")
                resLookup.append(res["tinkerLookupName"])
                for suffix in resLookup:
                    lookupName = atomLookup + "_" + suffix
                    if lookupName in bioTypes:
                        break
                if lookupName in bioTypes:
                    res["atoms"][atom]["type"] = bioTypes[lookupName][3]
                else:
                    print(f"For {atom} lookupName={lookupName} not in biotype")
                    if "parent" in res:
                        lookupName = f"{res['atoms'][atom]['tinkerLookupName']}_{res['parent']['tinkerLookupName']}"
                        if lookupName in bioTypes:
                            res["atoms"][atom]["type"] = bioTypes[lookupName][3]
                        else:
                            print(f"Missing lookupName={lookupName} from biotype")
        return 0

    @classmethod
    def processKeyFile(cls, keyFile, residueDict):
        """
        Process the .key or .prm file.

        Parameters
        ----------
        keyFile : str
            The name of the file.
        """
        # Get all interesting lines from the file
        allLines = []
        for line in open(keyFile):
            try:
                fields = shlex.split(line)
            except:
                continue
            if len(fields) == 0:
                continue
            if fields[0][0] == "#":
                continue
            allLines.append(fields)

        # Load lines in lists/scalar values
        atomTypes = dict()
        bioTypes = dict()
        forces = dict()
        scalars = cls.RECOGNIZED_SCALARS.copy()

        lineIndex = 0
        while lineIndex < len(allLines):
            fields = allLines[lineIndex]
            if fields[0] == "atom":
                element = elem.Element.getByAtomicNumber(int(fields[5])).symbol
                print("Element", element)
                atomTypes[int(fields[1])] = (fields[2], element, fields[6])
                lineIndex += 1

            elif fields[0] == "biotype":
                lookUp = f"{fields[2]}_{fields[3]}"
                if lookUp in bioTypes:
                    # Workaround for Tinker using the same name but different types for H2', H2'', and for H5', H5''
                    lookUp = f"{fields[2]}*_{fields[3]}"
                bioTypes[lookUp] = fields[1:]
                lineIndex += 1

            elif fields[0] in cls.RECOGNIZED_FORCES:
                if cls.RECOGNIZED_FORCES[fields[0]] == 1:
                    if fields[0] not in forces:
                        forces[fields[0]] = []
                    forces[fields[0]].append(fields[1:])
                    lineIndex += 1
                else:
                    lineIndex = cls.RECOGNIZED_FORCES[fields[0]](
                        lineIndex, allLines, forces
                    )
            elif fields[0] in cls.RECOGNIZED_SCALARS:
                scalars[fields[0]] = fields[1]
                lineIndex += 1
            else:
                print(f"Field {fields[0]} not recognized: line=<{allLines[lineIndex]}>")
                lineIndex += 1

        # Set biotypes for all atoms
        cls._setBioTypes(bioTypes, residueDict)

        return atomTypes, bioTypes, forces, scalars

    # ---------------------------------------------------------------------------------------- #
    #                                      WRITE XML FILE                                      #
    # ---------------------------------------------------------------------------------------- #
    @staticmethod
    def _writeAtomTypes(tinkerXmlFile, atomTypes):
        """
        Write the atom types to the XML file.

        Parameters
        ----------
        tinkerXmlFile : file
            The file to write to.
        atomTypes : dict
            The atom types dictionary.
        """
        tinkerXmlFile.write(" <AtomTypes>\n")
        for atomType in sorted(atomTypes):
            outputString = f'  <Type name="{atomType}" class="{atomTypes[atomType][0]}" element="{atomTypes[atomType][1]}" mass="{atomTypes[atomType][2]}" />'
            tinkerXmlFile.write(f"{outputString}\n")
        tinkerXmlFile.write(" </AtomTypes>\n")

    @staticmethod
    def _writeResidues(tinkerXmlFile, residueDict):
        tinkerXmlFile.write(" <Residues>\n")
        for resname, res in sorted(residueDict.items()):
            if res["include"]:
                if resname == "HOH":
                    # AMBER water is flexible
                    tinkerXmlFile.write(
                        f'  <Residue name="{resname}" rigidWater="false">\n'
                    )
                else:
                    tinkerXmlFile.write(f'  <Residue name="{resname}">\n')

                # Write atom fields
                atomIndex = dict()
                for atomCount, atom in enumerate(res["atoms"].keys()):
                    atomType = res["atoms"][atom]["type"]
                    if int(atomType) < 0:
                        print(
                            f"Error: type={atomType} for atom={atom} of residue={resname}"
                        )
                    tag = f'   <Atom name="{atom}" type="{atomType}" />'
                    atomIndex[atom] = atomCount
                    tinkerXmlFile.write(f"{tag}\n")

                # Write bond fields
                includedBonds = dict()
                for atomName in res["atoms"].keys():
                    bondsInfo = res["atoms"][atomName]["bonds"]
                    for bondedAtom in bondsInfo:
                        if (
                            bondedAtom not in includedBonds
                            or atomName not in includedBonds[bondedAtom]
                        ):
                            outputString = f'   <Bond from="{atomIndex[atomName]}" to="{atomIndex[bondedAtom]}" />'
                            if atomName not in includedBonds:
                                includedBonds[atomName] = dict()
                            if bondedAtom not in includedBonds:
                                includedBonds[bondedAtom] = dict()
                            includedBonds[atomName][bondedAtom] = 1
                            includedBonds[bondedAtom][atomName] = 1
                            tinkerXmlFile.write(f"{outputString}\n")

                tinkerXmlFile.write("  </Residue>\n")
        tinkerXmlFile.write(" </Residues>\n")

    @staticmethod
    def _writeAmoebaBondForce(tinkerXmlFile, forces, scalars):
        if "bond" not in forces:
            return

        cubic = 10.0 * float(scalars["bond-cubic"])
        quartic = 100.0 * float(scalars["bond-quartic"])
        outputString = (
            f""" <AmoebaBondForce bond-cubic="{cubic}" bond-quartic="{quartic}">"""
        )
        tinkerXmlFile.write(f"{outputString}\n")
        bonds = forces["bond"]
        for bond in bonds:
            length = float(bond[3]) * 0.1
            k = float(bond[2]) * 100.0 * 4.184
            outputString = f"""  <Bond class1="{bond[0]}" class2="{bond[1]}" length="{length}" k="{k}"/>"""
            tinkerXmlFile.write(f"{outputString}\n")
        tinkerXmlFile.write(" </AmoebaBondForce>\n")

    @staticmethod
    def _writeAmoebaAngleForce(tinkerXmlFile, forces, scalars):
        angleForces = [
            angleSet for angleSet in ["angle", "anglep"] if angleSet in forces
        ]
        if len(angleForces) == 0:
            return

        cubic = float(scalars["angle-cubic"])
        quartic = float(scalars["angle-quartic"])
        pentic = float(scalars["angle-pentic"])
        sextic = float(scalars["angle-sextic"])
        outputString = f""" <AmoebaAngleForce angle-cubic="{cubic}" angle-quartic="{quartic}" angle-pentic="{pentic}" angle-sextic="{sextic}">"""
        tinkerXmlFile.write(f"{outputString}\n")
        radian = 57.2957795130
        radian2 = 4.184 / (radian * radian)
        for angleSet in angleForces:
            for angle in forces[angleSet]:
                k = float(angle[3]) * radian2
                outputString = f'  <Angle class1="{angle[0]}" class2="{angle[1]}" class3="{angle[2]}" k="{k}" angle1="{angle[4]}"'
                if len(angle) > 5:
                    outputString += f' angle2="{angle[5]}"'
                if len(angle) > 6:
                    outputString += f' angle3="{angle[6]}"'
                outputString += f' inPlane="{angleSet == "anglep"}"/>'
                tinkerXmlFile.write(f"{outputString}\n")
        tinkerXmlFile.write(" </AmoebaAngleForce>\n")

    @staticmethod
    def _writeAmoebaOutOfPlaneBendForce(tinkerXmlFile, forces, scalars, radian):
        if "opbend" not in forces:
            return

        cubic = float(scalars["opbend-cubic"])
        quartic = float(scalars["opbend-quartic"])
        pentic = float(scalars["opbend-pentic"])
        sextic = float(scalars["opbend-sextic"])
        type = scalars["opbendtype"]
        outputString = f""" <AmoebaOutOfPlaneBendForce type="{type}" opbend-cubic="{cubic}" opbend-quartic="{quartic}" opbend-pentic="{pentic}" opbend-sextic="{sextic}">"""
        tinkerXmlFile.write(f"{outputString}\n")
        opbends = forces["opbend"]
        radian2 = 4.184 / (radian * radian)
        for opbend in opbends:
            k = float(opbend[4]) * radian2
            for i in range(4):
                if opbend[i] == "0":
                    opbend[i] = ""
            outputString = f"""  <Angle class1="{opbend[0]}" class2="{opbend[1]}" class3="{opbend[2]}" class4="{opbend[3]}" k="{k}"/>"""
            tinkerXmlFile.write(f"{outputString}\n")
        tinkerXmlFile.write(" </AmoebaOutOfPlaneBendForce>\n")

    @staticmethod
    def _writeAmoebaTorsionForce(tinkerXmlFile, forces, scalars, radian):
        if "torsion" not in forces:
            return

        torsionUnit = float(scalars["torsionunit"])
        outputString = """ <PeriodicTorsionForce>"""
        tinkerXmlFile.write(f"{outputString}\n")
        torsions = forces["torsion"]
        conversion = 4.184 * torsionUnit
        for torsion in torsions:
            outputString = f"""  <Proper class1="{torsion[0]}" class2="{torsion[1]}" class3="{torsion[2]}" class4="{torsion[3]}" """
            startIndex = 4
            for ii in range(0, 3):
                torsionSuffix = str(ii + 1)
                amplitudeAttributeName = f"k{torsionSuffix}"
                angleAttributeName = f"phase{torsionSuffix}"
                periodicityAttributeName = f"periodicity{torsionSuffix}"
                amplitude = float(torsion[startIndex]) * conversion
                angle = float(torsion[startIndex + 1]) / radian
                periodicity = int(torsion[startIndex + 2])
                outputString += f"""  {amplitudeAttributeName}="{amplitude}" {angleAttributeName}="{angle}" {periodicityAttributeName}="{periodicity}" """
                startIndex += 3
            outputString += "/>"
            tinkerXmlFile.write(f"{outputString}\n")
        tinkerXmlFile.write(" </PeriodicTorsionForce>\n")

    @staticmethod
    def _writeAmoebaPiTorsionForce(tinkerXmlFile, forces):
        if "pitors" not in forces:
            return

        piTorsionUnit = 1.0
        outputString = f""" <AmoebaPiTorsionForce piTorsionUnit="{piTorsionUnit}">"""
        tinkerXmlFile.write(f"{outputString}\n")
        piTorsions = forces["pitors"]
        conversion = 4.184 * piTorsionUnit
        for piTorsion in piTorsions:
            k = float(piTorsion[2]) * conversion
            outputString = f"""  <PiTorsion class1="{piTorsion[0]}" class2="{piTorsion[1]}" k="{k}" />"""
            tinkerXmlFile.write(f"{outputString}\n")
        tinkerXmlFile.write(" </AmoebaPiTorsionForce>\n")

    @staticmethod
    def _writeAmoebaStretchTorsionForce(tinkerXmlFile, forces):
        if "strtors" not in forces:
            return

        tinkerXmlFile.write(" <AmoebaStretchTorsionForce>\n")
        for torsion in forces["strtors"]:
            v = [float(x) * 10 * 4.184 for x in torsion[4:]]
            tinkerXmlFile.write(
                f'  <Torsion class1="{torsion[0]}" class2="{torsion[1]}" class3="{torsion[2]}" class4="{torsion[3]}" v11="{v[0]}" v12="{v[1]}" v13="{v[2]}" v21="{v[3]}" v22="{v[4]}" v23="{v[5]}" v31="{v[6]}" v32="{v[7]}" v33="{v[8]}"/>\n'
            )
        tinkerXmlFile.write(" </AmoebaStretchTorsionForce>\n")

    @staticmethod
    def _writeAmoebaAngleTorsionForce(tinkerXmlFile, forces):
        if "angtors" not in forces:
            return

        tinkerXmlFile.write(" <AmoebaAngleTorsionForce>\n")
        for torsion in forces["angtors"]:
            v = [float(x) * 4.184 for x in torsion[4:]]
            tinkerXmlFile.write(
                f'  <Torsion class1="{torsion[0]}" class2="{torsion[1]}" class3="{torsion[2]}" class4="{torsion[3]}" v11="{v[0]}" v12="{v[1]}" v13="{v[2]}" v21="{v[3]}" v22="{v[4]}" v23="{v[5]}"/>\n'
            )
        tinkerXmlFile.write(" </AmoebaAngleTorsionForce>\n")

    @staticmethod
    def _writeAmoebaStretchBendForce(tinkerXmlFile, forces, radian):
        if "strbnd" not in forces:
            return

        stretchBendUnit = 1.0
        outputString = (
            f""" <AmoebaStretchBendForce stretchBendUnit="{stretchBendUnit}">"""
        )
        tinkerXmlFile.write(f"{outputString}\n")
        conversion = 41.84 / radian
        stretchBends = forces["strbnd"]
        for stretchBend in stretchBends:
            k1 = float(stretchBend[3]) * conversion
            k2 = float(stretchBend[4]) * conversion
            outputString = f"""  <StretchBend class1="{stretchBend[0]}" class2="{stretchBend[1]}" class3="{stretchBend[2]}" k1="{k1}" k2="{k2}" />"""
            tinkerXmlFile.write(f"{outputString}\n")
        tinkerXmlFile.write("</AmoebaStretchBendForce>\n")

    @staticmethod
    def _writeAmoebaTorsionTorsionForce(tinkerXmlFile, forces):
        if "tortors" not in forces:
            return

        torsionTorsionUnit = 1.0
        tinkerXmlFile.write(""" <AmoebaTorsionTorsionForce >\n""")
        torsionTorsions = forces["tortors"]
        for index, torsionTorsion in enumerate(torsionTorsions):
            torInfo = torsionTorsion[0]
            grid = torsionTorsion[1]
            outputString = f"""  <TorsionTorsion class1="{torInfo[0]}" class2="{torInfo[1]}" class3="{torInfo[2]}" class4="{torInfo[3]}" class5="{torInfo[4]}" grid="{index}" nx="{torInfo[5]}" ny="{torInfo[6]}" />"""
            tinkerXmlFile.write(f"{outputString}\n")

        for index, torsionTorsion in enumerate(torsionTorsions):
            torInfo = torsionTorsion[0]
            grid = torsionTorsion[1]
            outputString = f"""  <TorsionTorsionGrid grid="{index}" nx="{torInfo[5]}" ny="{torInfo[6]}" >"""
            tinkerXmlFile.write(f"{outputString}\n")
            for gridIndex, gridEntry in enumerate(grid):
                print(f"Gxx {gridIndex}  {gridEntry}")
                if len(gridEntry) > 5:
                    f = float(gridEntry[2]) * 4.184
                    fx = float(gridEntry[3]) * 4.184
                    fy = float(gridEntry[4]) * 4.184
                    fxy = float(gridEntry[5]) * 4.184
                    outputString = f"""  <Grid angle1="{gridEntry[0]}" angle2="{gridEntry[1]}" f="{f}" fx="{fx}" fy="{fy}" fxy="{fxy}" />"""
                    tinkerXmlFile.write(f"  {outputString}\n")
                elif len(gridEntry) > 2:
                    f = float(gridEntry[2]) * 4.184
                    outputString = f"""  <Grid angle1="{gridEntry[0]}" angle2="{gridEntry[1]}" f="{f}" />"""
                    tinkerXmlFile.write(f"  {outputString}\n")
            tinkerXmlFile.write("</TorsionTorsionGrid >\n")

        tinkerXmlFile.write("</AmoebaTorsionTorsionForce>\n")

    @staticmethod
    def _writeAmoebaVdwForce(tinkerXmlFile, forces, scalars):
        if "vdw" not in forces and "vdwpr" not in forces:
            return

        outputString = f""" <AmoebaVdwForce type="{scalars['vdwtype']}" radiusrule="{scalars['radiusrule']}" radiustype="{scalars['radiustype']}" radiussize="{scalars['radiussize']}" epsilonrule="{scalars['epsilonrule']}" vdw-13-scale="{scalars['vdw-13-scale']}" vdw-14-scale="{scalars['vdw-14-scale']}" vdw-15-scale="{scalars['vdw-15-scale']}" >"""
        tinkerXmlFile.write(f"{outputString}\n")
        if "vdw" in forces:
            for vdw in forces["vdw"]:
                sigma = float(vdw[1]) * 0.1
                epsilon = float(vdw[2]) * 4.184
                reduction = vdw[3] if len(vdw) > 3 else 1.0
                outputString = f"""  <Vdw class="{vdw[0]}" sigma="{sigma}" epsilon="{epsilon}" reduction="{reduction}"/>"""
                tinkerXmlFile.write(f"{outputString}\n")

        if "vdwpr" in forces:
            for pair in forces["vdwpr"]:
                sigma = float(pair[2]) * 0.1
                epsilon = float(pair[3]) * 4.184
                outputString = f"""  <Pair class1="{pair[0]}" class2="{pair[1]}" sigma="{sigma}" epsilon="{epsilon}"/>"""
                tinkerXmlFile.write(f"{outputString}\n")
        tinkerXmlFile.write(" </AmoebaVdwForce>\n")

    @staticmethod
    def _writeAmoebaMultipoleForce(tinkerXmlFile, forces, recognizedScalars):
        scalarList = {
            "mpole12Scale": recognizedScalars["mpole-12-scale"],
            "mpole13Scale": recognizedScalars["mpole-13-scale"],
            "mpole14Scale": recognizedScalars["mpole-14-scale"],
            "mpole15Scale": recognizedScalars["mpole-15-scale"],
            "polar12Scale": recognizedScalars["polar-12-scale"],
            "polar13Scale": recognizedScalars["polar-13-scale"],
            "polar14Scale": recognizedScalars["polar-14-scale"],
            "polar15Scale": recognizedScalars["polar-15-scale"],
            "polar14Intra": recognizedScalars["polar-14-intra"],
            "direct11Scale": recognizedScalars["direct-11-scale"],
            "direct12Scale": recognizedScalars["direct-12-scale"],
            "direct13Scale": recognizedScalars["direct-13-scale"],
            "direct14Scale": recognizedScalars["direct-14-scale"],
            "mutual11Scale": recognizedScalars["mutual-11-scale"],
            "mutual12Scale": recognizedScalars["mutual-12-scale"],
            "mutual13Scale": recognizedScalars["mutual-13-scale"],
            "mutual14Scale": recognizedScalars["mutual-14-scale"],
        }

        outputString = " <AmoebaMultipoleForce "
        for key in sorted(scalarList.keys()):
            outputString += f' {key}="{scalarList[key]}" '
        outputString += " > "
        tinkerXmlFile.write(f"{outputString}\n")

        multipoleArray = forces["multipole"]
        bohr = 0.52917720859
        dipoleConversion = 0.1 * bohr
        quadrupoleConversion = 0.01 * bohr * bohr / 3.0
        for multipoleInfo in multipoleArray:
            axisInfo = multipoleInfo[0]
            multipoles = multipoleInfo[1]
            outputString = f'  <Multipole type="{axisInfo[0]}" '
            axisInfoLen = len(axisInfo)

            if axisInfoLen > 1:
                outputString += f'kz="{axisInfo[1]}" '

            if axisInfoLen > 2:
                outputString += f'kx="{axisInfo[2]}" '

            if axisInfoLen > 3:
                outputString += f'ky="{axisInfo[3]}" '

            outputString += f'c0="{multipoles[0]}" d1="{str(dipoleConversion * float(multipoles[1]))}" d2="{str(dipoleConversion * float(multipoles[2]))}" d3="{str(dipoleConversion * float(multipoles[3]))}" q11="{str(quadrupoleConversion * float(multipoles[4]))}" q21="{str(quadrupoleConversion * float(multipoles[5]))}" q22="{str(quadrupoleConversion * float(multipoles[6]))}" q31="{str(quadrupoleConversion * float(multipoles[7]))}" q32="{str(quadrupoleConversion * float(multipoles[8]))}" q33="{str(quadrupoleConversion * float(multipoles[9]))}"  />'
            tinkerXmlFile.write(f"{outputString}\n")

        polarizeArray = forces["polarize"]
        polarityConversion = 0.001
        m = {}
        for polarize in polarizeArray:
            m[polarize[0]] = []
            outputString = f'  <Polarize type="{polarize[0]}" polarizability="{str(polarityConversion * float(polarize[1]))}" thole="{polarize[2]}" '
            for ii in range(3, len(polarize)):
                outputString += f'pgrp{ii - 2}="{polarize[ii]}" '
                m[polarize[0]].append(polarize[ii])

            outputString += "/>"
            tinkerXmlFile.write(f"{outputString}\n")
            print(m[polarize[0]])
        for t in sorted(m):
            for k in m[t]:
                if t not in m[k]:
                    print(t, k)

        tinkerXmlFile.write(" </AmoebaMultipoleForce>\n")

    @staticmethod
    def _writeAmoebaGeneralizedKirkwoodForce(
        gkXmlFile, forces, atomTypes, bioTypes, recognizedScalars
    ):
        solventDielectric = 78.3
        soluteDielectric = 1.0
        includeCavityTerm = 1
        probeRadius = 0.14
        surfaceAreaFactor = -6.0 * 3.1415926535 * 0.0216 * 1000.0 * 0.4184
        outputString = f' <AmoebaGeneralizedKirkwoodForce solventDielectric="{solventDielectric}" soluteDielectric="{soluteDielectric}" includeCavityTerm="{includeCavityTerm}" probeRadius="{probeRadius}" surfaceAreaFactor="{surfaceAreaFactor}">'
        gkXmlFile.write(f"{outputString}\n")

        for type in sorted(atomTypes):
            print(f"atom type={type}  {atomTypes[type]}")

        for type in sorted(bioTypes):
            print(f"bio type={type}  {bioTypes[type]}")

        multipoleArray = forces["multipole"]
        for multipoleInfo in multipoleArray:
            axisInfo = multipoleInfo[0]
            multipoles = multipoleInfo[1]
            type = int(axisInfo[0])
            shct = 0.8
            if type in atomTypes:
                element = atomTypes[type][1]
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
                        f"Warning no overlap scale factor for type={type} element={element}"
                    )
            else:
                print(f"Warning no overlap scale factor for type={type} ")

            outputString = f'  <GeneralizedKirkwood type="{axisInfo[0]}" charge="{multipoles[0]}" shct="{shct}"  /> '
            gkXmlFile.write(f"{outputString}\n")
        gkXmlFile.write(" </AmoebaGeneralizedKirkwoodForce>\n")

    @staticmethod
    def _writeAmoebaWcaDispersionForce(gkXmlFile, forces, scalars):
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

        outputString = f' <AmoebaWcaDispersionForce epso="{epso * 4.184}" epsh="{epsh * 4.184}" rmino="{rmino * 0.1}" rminh="{rminh * 0.1}" awater="{1000.0 * awater}" slevy="{slevy}"  dispoff="{0.1 * dispoff}" shctd="{shctd}" >'
        gkXmlFile.write(f"{outputString}\n")
        vdws = forces["vdw"]
        convert = 0.1
        if scalars["radiustype"] == "SIGMA":
            convert *= 1.122462048309372

        if scalars["radiussize"] == "DIAMETER":
            convert *= 0.5

        for vdw in vdws:
            sigma = float(vdw[1])
            sigma *= convert
            epsilon = float(vdw[2]) * 4.184
            outputString = f'  <WcaDispersion class="{vdw[0]}" radius="{sigma}" epsilon="{epsilon}" /> '
            gkXmlFile.write(f"{outputString}\n")
        gkXmlFile.write(" </AmoebaWcaDispersionForce>\n")

    @staticmethod
    def _writeAmoebaUreyBradleyForce(tinkerXmlFile, forces):
        if "ureybrad" not in forces:
            return

        cubic = 0.0
        quartic = 0.0

        outputString = (
            f' <AmoebaUreyBradleyForce cubic="{cubic}" quartic="{quartic}"  >'
        )
        tinkerXmlFile.write(f"{outputString}\n")
        ubs = forces["ureybrad"]
        for ub in ubs:
            k = float(ub[3]) * 4.184 * 100.0
            d = float(ub[4]) * 0.1
            outputString = f'  <UreyBradley class1="{ub[0]}" class2="{ub[1]}" class3="{ub[2]}" k="{k}" d="{d}" /> '
            tinkerXmlFile.write(f"{outputString}\n")
        tinkerXmlFile.write(" </AmoebaUreyBradleyForce>\n")

    @staticmethod
    def _writeForceFieldEnd(tinkerXmlFile, gkXmlFile):
        tinkerXmlFile.write("</ForceField>\n")
        gkXmlFile.write("</ForceField>\n")
        tinkerXmlFile.close()
        gkXmlFile.close()

    @classmethod
    def createXmlFiles(
        cls, filename, atomTypes, bioTypes, forces, residueDict, scalars
    ):
        """
        Create the XML file.

        Parameters
        ----------


        Returns
        -------
        tinkXmlFile, gkXmlFile : tuple(file, file)
            The XML files.
        """
        if "forcefield" not in scalars:
            scalars["forcefield"] = filename.rsplit(".", 1)[0]

        tinkerXmlFileName = scalars["forcefield"]
        tinkerXmlFileName += ".xml"
        tinkerXmlFile = open(tinkerXmlFileName, "w")
        print(f"Opened {tinkerXmlFileName}.")

        gkXmlFileName = scalars["forcefield"]
        gkXmlFileName += "_gk.xml"
        gkXmlFile = open(gkXmlFileName, "w")
        print(f"Opened {gkXmlFileName}.")

        # Write the header
        today = datetime.date.today().isoformat()
        sourceFile = os.path.basename(filename)
        header = f""" <Info>
        <Source>{sourceFile}</Source>
        <DateGenerated>{today}</DateGenerated>
        <Reference></Reference>
        </Info>
    """

        gkXmlFile.write("<ForceField>\n")
        gkXmlFile.write(header)
        tinkerXmlFile.write("<ForceField>\n")
        tinkerXmlFile.write(header)

        cls._writeAtomTypes(tinkerXmlFile, atomTypes)
        cls._writeResidues(tinkerXmlFile, residueDict)
        # if bioTypes:
        #    cls._writeEndCaps(tinkerXmlFile, bioTypes)
        # cls._writeIons(tinkerXmlFile, ions)

        radian = 57.2957795130
        cls._writeAmoebaBondForce(tinkerXmlFile, forces, scalars)
        cls._writeAmoebaAngleForce(tinkerXmlFile, forces, scalars)
        cls._writeAmoebaOutOfPlaneBendForce(tinkerXmlFile, forces, scalars, radian)
        cls._writeAmoebaTorsionForce(tinkerXmlFile, forces, scalars, radian)
        cls._writeAmoebaPiTorsionForce(tinkerXmlFile, forces)
        cls._writeAmoebaStretchTorsionForce(tinkerXmlFile, forces)
        cls._writeAmoebaAngleTorsionForce(tinkerXmlFile, forces)
        cls._writeAmoebaStretchBendForce(tinkerXmlFile, forces, radian)
        cls._writeAmoebaTorsionTorsionForce(tinkerXmlFile, forces)
        cls._writeAmoebaVdwForce(tinkerXmlFile, forces, scalars)
        cls._writeAmoebaMultipoleForce(tinkerXmlFile, forces, scalars)
        cls._writeAmoebaGeneralizedKirkwoodForce(
            gkXmlFile, forces, atomTypes, bioTypes, scalars
        )
        cls._writeAmoebaWcaDispersionForce(gkXmlFile, forces, scalars)
        cls._writeAmoebaUreyBradleyForce(tinkerXmlFile, forces)
        cls._writeForceFieldEnd(tinkerXmlFile, gkXmlFile)

        return tinkerXmlFile, gkXmlFile
