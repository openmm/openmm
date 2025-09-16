"""
tinkerfiles.py: A reader of Tinker files (.xyz, .prm, .key, .prm).

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

import math
import re
import shlex
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from openmm.app.internal.unitcell import computePeriodicBoxVectors
from openmm.unit import nanometers
from openmm.vec3 import Vec3
from openmm.app import element as elem

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

    def updateFromAtomType(self, atomTypeData: 'TinkerAtomType') -> None:
        """
        Update atom data from atom type information.

        Parameters
        ----------
        atomTypeData : TinkerAtomType
            TinkerAtomType object containing atom type data
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

    Attributes
    ----------
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
    """TinkerFiles parses Tinker files (.xyz, .prm, .key), constructs a Topology, and (optionally) an OpenMM System from it.
    This class only supports the AMOEBA force field.  It cannot create a System from Tinker files that use other force fields."""

    @staticmethod
    def _initialize_class() -> Tuple[Dict[str, Any], Dict[str, str], Dict[str, float], Dict[str, float]]:
        """
        Initialize class variables upon class creation.
        
        Returns
        -------
        Tuple[Dict[str, Any], Dict[str, str], Dict[str, float], Dict[str, float]]
            A tuple containing RECOGNIZED_FORCES, RECOGNIZED_SCALARS, WCA_PARAMS, and GK_PARAMS
        """

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
            bohr = 0.52917720859
            dipoleConversion = 0.1*bohr
            quadrupoleConversion = 0.01*bohr*bohr/3.0
            fields = allLines[lineIndex]
            type = int(fields[1])
            axis = [int(x) for x in fields[2:-1]]
            charge = float(fields[-1])
            lineIndex += 1
            fields = allLines[lineIndex]
            dipole = [float(x)*dipoleConversion for x in fields[:3]]
            lineIndex += 1
            fields = allLines[lineIndex]
            quadrupole = [fields[0]]
            lineIndex += 1
            fields = allLines[lineIndex]
            quadrupole.append(fields[0])
            quadrupole.append(fields[1])
            lineIndex += 1
            fields = allLines[lineIndex]
            quadrupole.append(fields[0])
            quadrupole.append(fields[1])
            quadrupole.append(fields[2])
            lineIndex += 1
            quadrupole = [float(x)*quadrupoleConversion for x in quadrupole]
            forces["multipole"].append([type, axis, charge, dipole, quadrupole])
            return lineIndex

        def addPolarize(
            lineIndex: int, allLines: List[List[str]], forces: Dict[str, Any]
        ) -> int:
            """
            Parse and store polarization data from the key file.

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
            if "polarize" not in forces:
                forces["polarize"] = []
            fields = allLines[lineIndex]
            type = int(fields[1])
            polarizability = float(fields[2])
            thole = float(fields[3])
            group = [int(x) for x in fields[4:]]
            lineIndex += 1
            forces["polarize"].append([type, polarizability, thole, group])
            return lineIndex

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
            "anglep": 1,
            "strbnd": 1,
            "ureybrad": 1,
            "opbend": 1,
            "torsion": 1,
            "pitors": 1,
            "strtors": 1,
            "angtors": 1,
            "vdw": 1,
            "vdwpair": 1,
            "vdwpr": 1,
            "multipole": addMultipole,
            "polarize": addPolarize,
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
            "pitorsunit": "1.0",
            "angtorunit": "1.0",
            "strtorunit": "1.0",
            "torsionunit": "1.0",
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

        WCA_PARAMS = {
            "epso": 0.1100 * 4.184,
            "epsh": 0.0135 * 4.184,
            "rmino": 1.7025 * 0.1,
            "rminh": 1.3275 * 0.1,
            "awater": 0.033428 * 1000.0,
            "slevy": 1.0,
            "dispoff": 0.26 * 0.1,
            "shctd": 0.81,
        }
        
        GK_PARAMS = {
            "solventDielectric": 78.3,
            "soluteDielectric": 1.0,
            "includeCavityTerm": 1,
            "probeRadius": 0.14,
            "surfaceAreaFactor": -6.0 * 3.1415926535 * 0.0216 * 1000.0 * 0.4184,
        }

        return RECOGNIZED_FORCES, RECOGNIZED_SCALARS, WCA_PARAMS, GK_PARAMS

    # Call the initialization method
    RECOGNIZED_FORCES, RECOGNIZED_SCALARS, WCA_PARAMS, GK_PARAMS = _initialize_class()

    def __init__(
        self,
        xyzFile: str,
        parameterFiles: Union[str, List[str]],
        periodicBoxVectors: Optional[Tuple[Vec3, Vec3, Vec3]] = None,
        unitCellDimensions: Optional[Vec3] = None,
    ):
        """
        Load a set of Tinker files, including one .xyz files and one or more .key or .prm files.

        Parameters
        ----------
        xyzFile : str
            Path to the Tinker .xyz file.
        parameterFiles : str or list[str]
            Paths to the Tinker .key and .prm files.
        periodicBoxVectors : Optional[Tuple[Vec3, Vec3, Vec3]]
            The periodic box vectors.
        unitCellDimensions : Optional[Vec3]
            The unit cell dimensions.
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

        # Load the xyz file
        self.atoms, self.boxVectors, self.positions = TinkerFiles._loadXyzFile(xyzFile)
        self.positions = self.positions * nanometers

        # Load the .key or .prm file(s)
        parameterFiles = [parameterFiles] if isinstance(parameterFiles, str) else parameterFiles
        for paramFile in parameterFiles:
            atomTypes, forces, scalars = TinkerFiles._loadKeyFile(paramFile)
            self._atomTypes.update(atomTypes)
            self._scalars.update(scalars)

            # Extend forces lists if force type already exists
            for key, value in forces.items():
                if key in self._forces:
                    self._forces[key].extend(value)
                else:
                    self._forces[key] = value

        # Update atoms with atom type information
        for atom in self.atoms:
            if atom.atomType not in self._atomTypes:
                raise ValueError(f"Atom type {atom.atomType} not found. Check that the parameter files are correct.")
            atom.updateFromAtomType(self._atomTypes[atom.atomType])

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

    def createSystem(
        self,
        nonbondedMethod=ff.NoCutoff,
        nonbondedCutoff=1.0*nanometers,
        vdwCutoff=None,
        removeCMMotion: bool = True,
        hydrogenMass=None,
        polarization: str = "mutual",
        mutualInducedTargetEpsilon: float = 0.00001,
        implicitSolvent: bool = False,
        ewaldErrorTolerance = 0.0005,
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
        vdwCutoff : distance=None
            An optional alternate cutoff to use for vdw interactions.  If this is omitted, nonbondedCutoff
            is used for vdw interactions as well as for other nonbonded interactions.
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
        ewaldErrorTolerance : float=0.0005
            The error tolerance to use if nonbondedMethod is Ewald, PME, or LJPME.

        Returns
        -------
        openmm.System
            The created OpenMM System.
        """
        from openmm.app.internal.amoebaforces import (
            AmoebaBondForceBuilder,
            AmoebaUreyBradleyForceBuilder,
            AmoebaAngleForceBuilder,
            AmoebaOutOfPlaneBendForceBuilder,
            AmoebaInPlaneAngleForceBuilder,
            AmoebaStretchBendForceBuilder,
            AmoebaTorsionForceBuilder,
            AmoebaPiTorsionForceBuilder,
            AmoebaStretchTorsionForceBuilder,
            AmoebaAngleTorsionForceBuilder,
            AmoebaTorsionTorsionForceBuilder,
            AmoebaWcaDispersionForceBuilder,
            AmoebaGeneralizedKirkwoodForceBuilder,
            AmoebaMultipoleForceBuilder,
            AmoebaVdwForceBuilder,
            MultipoleParams,
            PolarizationParams
        )
        import openmm as mm

        sys = mm.System()
        for atom in self.atoms:
            sys.addParticle(atom.mass)

        # Mapping from atom index to atom class used by various force builders
        atomClasses = {}
        for atom_idx in range(len(self.atoms)):
            atomClasses[atom_idx] = self.atoms[atom_idx].atomClass

        # List of bonds as (atom1_index, atom2_index) used by various force builders
        bonds = [(at1.index, at2.index) for at1, at2 in self.topology.bonds()]

        # Add AmoebaBondForce
        if "bond" in self._forces:
            bondParams = {(at1, at2): {"k": float(k), "r0": float(r0)} for at1, at2, k, r0 in self._forces["bond"]}
            bondForceBuilder = AmoebaBondForceBuilder(float(self._scalars["bond-cubic"])*10, float(self._scalars["bond-quartic"])*100)
            for (class1, class2), params in bondParams.items():
                bondForceBuilder.registerParams((class1, class2), params["r0"]*0.1, params["k"]*100.0*4.184)
            bondForce = bondForceBuilder.getForce(sys)
            bondForceBuilder.addBonds(bondForce, atomClasses, bonds)

        if "angle" in self._forces or "opbend" in self._forces or "anglep" in self._forces:
            # Find all unique angles
            bondedToAtom = [[] for _ in range(self.topology.getNumAtoms())]
            for (at1, at2) in bonds:
                bondedToAtom[at1].append(at2)
                bondedToAtom[at2].append(at1)

            uniqueAngles = set()
            for atom in range(len(bondedToAtom)):
                neighbors = bondedToAtom[atom]
                for i, n1 in enumerate(neighbors):
                    for n2 in neighbors[i+1:]:  
                        angle = (min(n1, n2), atom, max(n1, n2))
                        uniqueAngles.add(angle)
            angles = sorted(list(uniqueAngles))

            # Need to initialize here as used to create out-of-plane angles
            if "opbend" in self._forces:
                opbendParams = {(at1, at2, at3, at4): {"k": float(k)} for at1, at2, at3, at4, k in self._forces["opbend"]}
            else:
                opbendParams = {}

            # Classify angles into in-plane, out-of-plane and generic.
            # If middle atom has covalency of 3 and the types of the middle atom and the partner atom (atom bonded to
            # middle atom, but not in angle) match types1 and types2, then three out-of-plane bend angles are generated.
            # Three in-plane angle are also generated.
            # If the conditions are not satisfied, the angle is marked as 'generic' angle (not a in-plane angle).
            genericAngles = []
            inPlaneAngles = []
            outOfPlaneAngles = []
            skipAtoms = {}
            for angle in angles:
                middleAtom = angle[1]
                middleClass = atomClasses[middleAtom] 
                middleValence = len(self.atoms[middleAtom].bonds)
                if middleValence == 3 and middleAtom not in skipAtoms:
                    partners = []
                    for partner in self.atoms[middleAtom].bonds:
                        partnerClass = atomClasses[partner]
                        for opBendClass in opbendParams:
                            if middleClass == opBendClass[1] and partnerClass == opBendClass[0]:
                                partners.append(partner)
                                break
                    if len(partners) == 3:
                        outOfPlaneAngles.append([partners[0], middleAtom, partners[1], partners[2]])
                        outOfPlaneAngles.append([partners[2], middleAtom, partners[0], partners[1]])
                        outOfPlaneAngles.append([partners[1], middleAtom, partners[2], partners[0]])

                        skipAtoms[middleAtom] = set(partners[:3])

                        angleList = list(angle[:3])
                        for atomIndex in partners:
                            if atomIndex not in angleList:
                                angleList.append(atomIndex)

                        inPlaneAngles.append(angleList)
                    else:
                        angleList = list(angle[:3])
                        for atomIndex in partners:
                            if atomIndex not in angleList:
                                angleList.append(atomIndex)
                        genericAngles.append(angleList)
                elif middleValence == 3 and middleAtom in skipAtoms:
                    angleList = list(angle[:3])
                    for atomIndex in skipAtoms[middleAtom]:
                        if atomIndex not in angleList:
                            angleList.append(atomIndex)
                    inPlaneAngles.append(angleList)
                else:
                    genericAngles.append(list(angle))

            # Sanity checks for angle classification
            assert len(inPlaneAngles) + len(genericAngles)== len(angles), (
                f"Angle classification mismatch:\n"
                f"  in-plane: {len(inPlaneAngles)}\n"
                f"  generic: {len(genericAngles)}\n"
                f"  total classified: {len(inPlaneAngles) + len(genericAngles)}\n"
                f"  expected total: {len(angles)}"
            )
            assert len(inPlaneAngles) == len(outOfPlaneAngles), (
                f"Number of in-plane and out-of-plane angles do not match:\n"
                f"  in-plane: {len(inPlaneAngles)}\n"
                f"  out-of-plane: {len(outOfPlaneAngles)}"
            )

            # Add AmoebaOutOfPlaneBendForce
            if outOfPlaneAngles:
                outOfPlaneBendForceBuilder = AmoebaOutOfPlaneBendForceBuilder(float(self._scalars["opbend-cubic"]), float(self._scalars["opbend-quartic"]), float(self._scalars["opbend-pentic"]), float(self._scalars["opbend-sextic"]))
                for (class1, class2, class3, class4), params in opbendParams.items():
                    outOfPlaneBendForceBuilder.registerParams((class1, class2, class3, class4), (params["k"]*4.184*(math.pi/180.0)**2,))
                outOfPlaneBendForce = outOfPlaneBendForceBuilder.getForce(sys)
                outOfPlaneBendForceBuilder.addOutOfPlaneBends(outOfPlaneBendForce, atomClasses, outOfPlaneAngles)

            # Add AmoebaAngleForce
            idealAngles = {}
            angleParams = {(at1, at2, at3): {"k": float(k), "theta0": [float(theta) for theta in theta]} for at1, at2, at3, k, *theta in self._forces["angle"]}
            if genericAngles:
                angleForceBuilder = AmoebaAngleForceBuilder(self._scalars["angle-cubic"], self._scalars["angle-quartic"], self._scalars["angle-pentic"], self._scalars["angle-sextic"])
                processedParams = set()
                processedAngles = []
                for angle in genericAngles:
                    class1 = atomClasses[angle[0]]
                    class2 = atomClasses[angle[1]]
                    class3 = atomClasses[angle[2]]
                    params = angleParams.get((class1, class2, class3)) or angleParams.get((class3, class2, class1))
                    if params:
                        if len(params["theta0"]) > 1:
                            # Get k-index by counting number of non-angle hydrogens on the central atom
                            partners = [at for at in bondedToAtom[angle[1]] if at not in angle]
                            nHyd = sum(1 for i in partners if self.atoms[i].atomicNumber == 1)
                            if nHyd < len(params["theta0"]):
                                theta0 = params["theta0"][nHyd]
                            else:
                                raise ValueError(f"Angle parameters out of range for atom classes {class1}-{class2}-{class3}")
                        else:
                            theta0 = params["theta0"][0]
                        angleTuple = (angle[0], angle[1], angle[2])
                        idealAngles[angleTuple] = theta0 # Store ideal angle for stretch-bend
                        angleForceBuilder.registerParams(angleTuple, (theta0, params["k"]*4.184*(math.pi/180)**2))
                        processedAngles.append(angleTuple)
                if processedAngles:
                    angleForce = angleForceBuilder.getForce(sys)
                    angleForceBuilder.addAngles(angleForce, processedAngles)

            # Add AmoebaInPlaneAngleForce
            if "anglep" in self._forces:
                inPlaneAngleParams = {(at1, at2, at3): {"k": float(k), "theta0": [float(theta) for theta in theta]} for at1, at2, at3, k, *theta in self._forces["anglep"]}
            else:
                inPlaneAngleParams = {}
    
            inPlaneAngleParams.update(angleParams) # Poltype e.g. does not use the "anglep" keyword, but uses "angle" keyword for in-plane angles

            if inPlaneAngles:
                inPlaneAngleForceBuilder = AmoebaInPlaneAngleForceBuilder(self._scalars["angle-cubic"], self._scalars["angle-quartic"], self._scalars["angle-pentic"], self._scalars["angle-sextic"])
                processedParams = set()
                for angle in inPlaneAngles:
                    class1 = atomClasses[angle[0]]
                    class2 = atomClasses[angle[1]]
                    class3 = atomClasses[angle[2]]
                    params = inPlaneAngleParams.get((class1, class2, class3)) or inPlaneAngleParams.get((class3, class2, class1))
                    if params:
                        class4 = atomClasses[angle[3]]
                        paramKey = (class1, class2, class3, class4)
                        if len(params["theta0"]) > 1:
                            # Get k-index by counting number of non-angle hydrogens on the central atom
                            partners = [at for at in bondedToAtom[angle[1]] if at not in angle[:]]
                            nHyd = sum(1 for i in partners if self.atoms[i].atomicNumber == 1)
                            if nHyd < len(params["theta0"]):
                                theta0 = params["theta0"][nHyd]
                            else:
                                raise ValueError(f"Angle parameters out of range for atom classes {class1}-{class2}-{class3}")
                        else:
                            theta0 = params["theta0"][0]
                        idealAngles[tuple(angle)] = theta0 # Store ideal angle for stretch-bend
                        if paramKey not in processedParams:
                            inPlaneAngleForceBuilder.registerParams(paramKey, (theta0, params["k"]*4.184*(math.pi/180)**2))
                            processedParams.add(paramKey)
                if processedParams:
                    inPlaneAngleForce = inPlaneAngleForceBuilder.getForce(sys)
                    inPlaneAngleForceBuilder.addInPlaneAngles(inPlaneAngleForce, atomClasses, inPlaneAngles)

        # Add AmoebaStretchBendForce
        if "strbnd" in self._forces:
            stretchBendParams = {(at1, at2, at3): {"k1": float(k1), "k2": float(k2)} for at1, at2, at3, k1, k2 in self._forces["strbnd"]}
            stretchBendForceBuilder = AmoebaStretchBendForceBuilder()
            processedAngles = []
            for angle in genericAngles + inPlaneAngles:
                class1 = atomClasses[angle[0]]
                class2 = atomClasses[angle[1]]
                class3 = atomClasses[angle[2]]

                params = stretchBendParams.get((class1, class2, class3))
                if params is not None:
                    swap = False
                else:
                    params = stretchBendParams.get((class3, class2, class1))
                    swap = params is not None
           
                if params:
                    idealAngle = idealAngles.get(tuple(angle), None) 
                    if idealAngle is None:
                        continue

                    if not swap:
                        # 1-2-3
                        bondAB = bondParams.get((class1, class2)) or bondParams.get((class2, class1))
                        bondCB = bondParams.get((class2, class3)) or bondParams.get((class3, class2))
                        if bondAB and bondCB:
                            bondAB, bondCB = bondAB['r0'], bondCB['r0']
                            k1, k2 = params["k1"], params["k2"]
                    else:
                        # 3-2-1
                        bondAB = bondParams.get((class3, class2)) or bondParams.get((class2, class3))
                        bondCB = bondParams.get((class2, class1)) or bondParams.get((class1, class2))
                        if bondAB and bondCB:
                            bondAB, bondCB = bondCB['r0'], bondAB['r0']
                            k1, k2 = params["k2"], params["k1"]

                    if bondAB and bondCB:
                        # Because for the same class triplet, depending on the number of hydrogens on the central atom,
                        # the ideal angle can be different, we need to register parameters for each angle.
                        paramKey = (angle[0], angle[1], angle[2])
                        stretchBendForceBuilder.registerParams(paramKey, (bondAB*0.1, bondCB*0.1, idealAngle*math.pi/180, k1*41.84*math.pi/180, k2*41.84*math.pi/180))
                        processedAngles.append((angle[0], angle[1], angle[2]))
        
            if processedAngles:
                stretchBendForce = stretchBendForceBuilder.getForce(sys)
                stretchBendForceBuilder.addStretchBends(stretchBendForce, processedAngles)
     
        # Add AMOEBA Urey-Bradley force
        if "ureybrad" in self._forces:
            ureyBradleyParams = {(at1, at2, at3): {"k": float(k), "d": float(d)} for at1, at2, at3, k, d in self._forces["ureybrad"]}
            ureyBradleyForceBuilder = AmoebaUreyBradleyForceBuilder()

            for (class1, class2, class3), params in ureyBradleyParams.items():
                ureyBradleyForceBuilder.registerParams((class1, class2, class3), (params["k"]*4.184*100.0, params["d"]*0.1))

            processedAngles = []
            for at1, at2, at3 in angles:
                class1 = atomClasses[at1]
                class2 = atomClasses[at2]
                class3 = atomClasses[at3]
                params = ureyBradleyParams.get((class1, class2, class3)) or ureyBradleyParams.get((class3, class2, class1))
                if params:
                    processedAngles.append((at1, at2, at3))
            if processedAngles:
                ureyBradleyForce = ureyBradleyForceBuilder.getForce(sys)
                ureyBradleyForceBuilder.addUreyBradleys(ureyBradleyForce, atomClasses, processedAngles)

        # Find all unique proper torsions in the system
        uniquePropers = set()
        for angle in angles:
            for atom in self.atoms[angle[0]].bonds:
                if atom not in angle:
                    if atom < angle[2]:
                        uniquePropers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniquePropers.add((angle[2], angle[1], angle[0], atom))
            for atom in self.atoms[angle[2]].bonds:
                if atom not in angle:
                    if atom > angle[0]:
                        uniquePropers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniquePropers.add((atom, angle[2], angle[1], angle[0]))
        propers = sorted(list(uniquePropers))

        # Add AmoebaTorsionForce
        if "torsion" in self._forces:
            torsionScale = float(self._scalars["torsionunit"])*4.184
            torsionParams = {(at1, at2, at3, at4): {"t1": [float(t11)*torsionScale, float(t12)*math.pi/180.0, int(t13)],
                                                    "t2": [float(t21)*torsionScale, float(t22)*math.pi/180.0, int(t23)],
                                                    "t3": [float(t31)*torsionScale, float(t32)*math.pi/180.0, int(t33)]}
                            for at1, at2, at3, at4, t11, t12, t13, t21, t22, t23, t31, t32, t33 in self._forces["torsion"]}

            torsionForceBuilder = AmoebaTorsionForceBuilder()
            for (class1, class2, class3, class4), params in torsionParams.items():
                torsionForceBuilder.registerParams((class1, class2, class3, class4), (params["t1"], params["t2"], params["t3"]))

            processedTorsions = []
            for at1, at2, at3, at4 in propers:
                class1 = atomClasses[at1]
                class2 = atomClasses[at2]
                class3 = atomClasses[at3]
                class4 = atomClasses[at4]
                params = torsionParams.get((class1, class2, class3, class4)) or torsionParams.get((class4, class3, class2, class1))
                if params:
                    processedTorsions.append((at1, at2, at3, at4))

            if processedTorsions:
                torsionForce = torsionForceBuilder.getForce(sys)
                atomClassesDict = {}
                for atom_idx in range(len(self.atoms)):
                    atomClassesDict[atom_idx] = atomClasses[atom_idx]
                torsionForceBuilder.addTorsions(torsionForce, atomClasses, processedTorsions)
    
        # Add AmoebaPiTorsionForce
        if "pitors" in self._forces:
            piTorsionParams = {(at1, at2): {"k": float(k)} for at1, at2, k in self._forces["pitors"]}
            piTorsionForceBuilder = AmoebaPiTorsionForceBuilder()

            piTorsionScale = float(self._scalars["pitorsunit"])*4.184
            for (class1, class2), params in piTorsionParams.items():
                piTorsionForceBuilder.registerParams((class1, class2), (params["k"]*piTorsionScale,))

            processedPiTorsions = []
            for (at1, at2) in bonds:
                valence1 = len(self.atoms[at1].bonds)
                valence2 = len(self.atoms[at2].bonds)
                if valence1 == 3 and valence2 == 3:
                    class1 = atomClasses[at1]
                    class2 = atomClasses[at2]
                    params = piTorsionParams.get((class1, class2)) or piTorsionParams.get((class2, class1))
                    if params:
                        piTorsionAtom3 = at1
                        piTorsionAtom4 = at2
                        # piTorsionAtom1, piTorsionAtom2 are the atoms bonded to atom1, excluding atom2
                        # piTorsionAtom5, piTorsionAtom6 are the atoms bonded to atom2, excluding atom1
                        piTorsionAtom1, piTorsionAtom2 = [bond for bond in self.atoms[at1].bonds if bond != piTorsionAtom4]
                        piTorsionAtom5, piTorsionAtom6 = [bond for bond in self.atoms[at2].bonds if bond != piTorsionAtom3]
                        processedPiTorsions.append((piTorsionAtom1, piTorsionAtom2, piTorsionAtom3, piTorsionAtom4, piTorsionAtom5, piTorsionAtom6))
            
            if processedPiTorsions:
                piTorsionForce = piTorsionForceBuilder.getForce(sys)
                piTorsionForceBuilder.addPiTorsions(piTorsionForce, atomClassesDict, processedPiTorsions)

        # Add AmoebaStretchTorsionForce
        if "strtors" in self._forces:
            stretchTorsionScale = float(self._scalars["strtorunit"])*10*4.184
            stretchTorsionParams = {(p[0], p[1], p[2], p[3]): [float(x)*stretchTorsionScale for x in p[4:]] for p in self._forces["strtors"]}
            stretchTorsionForceBuilder = AmoebaStretchTorsionForceBuilder()
            for classes, params in stretchTorsionParams.items():
                stretchTorsionForceBuilder.registerParams(classes, params)
            processedStretchTorsions = []
            for at1, at2, at3, at4 in propers:
                class1 = atomClasses[at1]
                class2 = atomClasses[at2]
                class3 = atomClasses[at3]
                class4 = atomClasses[at4]
                params = stretchTorsionParams.get((class1, class2, class3, class4)) or stretchTorsionParams.get((class4, class3, class2, class1))
                if params:
                    processedStretchTorsions.append((at1, at2, at3, at4))
            if processedStretchTorsions:
                stretchTorsionForce = stretchTorsionForceBuilder.getForce(sys)
                stretchTorsionForceBuilder.addStretchTorsions(sys, stretchTorsionForce, atomClasses, processedStretchTorsions)

        # Add AmoebaAngleTorsionForce
        if "angtors" in self._forces:
            angleTorsionScale = float(self._scalars["angtorunit"])*4.184
            angleTorsionParams = {(p[0], p[1], p[2], p[3]): [float(x)*angleTorsionScale for x in p[4:]] for p in self._forces["angtors"]}
            angleTorsionForceBuilder = AmoebaAngleTorsionForceBuilder()
            for classes, params in angleTorsionParams.items():
                angleTorsionForceBuilder.registerParams(classes, params)
            processedAngleTorsions = []
            for at1, at2, at3, at4 in propers:
                class1 = atomClasses[at1]
                class2 = atomClasses[at2]
                class3 = atomClasses[at3]
                class4 = atomClasses[at4]
                params = angleTorsionParams.get((class1, class2, class3, class4)) or angleTorsionParams.get((class4, class3, class2, class1))
                if params:
                    processedAngleTorsions.append((at1, at2, at3, at4))
            if processedAngleTorsions:
                angleTorsionForce = angleTorsionForceBuilder.getForce(sys)
                angleTorsionForceBuilder.addAngleTorsions(sys, angleTorsionForce, atomClasses, processedAngleTorsions)

        # Add AmoebaTorsionTorsionForce
        if "tortors" in self._forces:
            torsionTorsionForce = AmoebaTorsionTorsionForceBuilder()
            for gridIndex, (tortorInfo, gridData) in enumerate(self._forces["tortors"]):
                nx = int(tortorInfo[5]) 
                ny = int(tortorInfo[6])  
                grid = np.array(gridData, dtype=np.float64).reshape((nx, ny, -1))
                grid[:, :, 2] *= 4.184  
                torsionTorsionForce.registerGridData(gridIndex, grid)
            force = torsionTorsionForce.getForce(sys)
            AmoebaTorsionTorsionForceBuilder.createTorsionTorsionInteractions(force, angles, self.atoms, self._forces["tortors"])

        # Add AmoebaVdwForce
        if "vdw" in self._forces:
            vdwForceBuilder = AmoebaVdwForceBuilder(self._scalars['vdwtype'], self._scalars['radiusrule'], self._scalars['radiustype'],
                                            self._scalars['radiussize'], self._scalars['epsilonrule'], float(self._scalars['vdw-13-scale']),
                                            float(self._scalars['vdw-14-scale']), float(self._scalars['vdw-15-scale']))
            for params in self._forces['vdw']:
                reduction = 0.0 if len(params) < 4 else float(params[3])
                vdwForceBuilder.registerClassParams(int(params[0]), float(params[1])*0.1, float(params[2])*4.184, reduction)
            pairs = self._forces.get('vdwpair', None)
            if pairs is None:
                pairs = self._forces.get('vdwpr', None)
            if pairs is not None:
                for params in pairs:
                    vdwForceBuilder.registerPairParams(int(params[0]), int(params[1]), float(params[2])*0.1, float(params[3])*4.184)
            vdwForce = vdwForceBuilder.getForce(sys, nonbondedMethod, nonbondedCutoff if vdwCutoff is None else vdwCutoff, True)
            atomClasses = [int(atom.atomClass) for atom in self.atoms]
            vdwForceBuilder.addParticles(vdwForce, atomClasses, list(self.topology.atoms()), bonds)

        # Add AmoebaMultipoleForce
        if "multipole" in self._forces:
            multipoleForceBuilder = AmoebaMultipoleForceBuilder()
            for atType, axis, charge, dipole, quadrupole in self._forces["multipole"]:
                kIndices = [atType] + axis
                q = [quadrupole[i] for i in (0, 1, 3, 1, 2, 4, 3, 4, 5)]
                multipoleForceBuilder.registerMultipoleParams(atType, MultipoleParams(kIndices, charge, dipole, q))
            for atType, polarizability, thole, group in self._forces["polarize"]:
                multipoleForceBuilder.registerPolarizationParams(atType, PolarizationParams(polarizability*0.001, thole, group))
            multipoleForce = multipoleForceBuilder.getForce(sys, nonbondedMethod, nonbondedCutoff, ewaldErrorTolerance, polarization, mutualInducedTargetEpsilon, 60)
            atomTypes = [int(atom.atomType) for atom in self.atoms]
            multipoleForceBuilder.addMultipoles(multipoleForce, atomTypes, list(self.topology.atoms()), bonds)

        # Add AmoebaGeneralizedKirkwoodForce
        if implicitSolvent and "multipole" in self._forces:
            gkForceBuilder = AmoebaGeneralizedKirkwoodForceBuilder(**self.GK_PARAMS)
            force = gkForceBuilder.getForce(sys, implicitSolvent=implicitSolvent)
            multipoleForce = [f for f in sys.getForces() if isinstance(f, mm.AmoebaMultipoleForce)][0]
            for atomIndex in range(0, multipoleForce.getNumMultipoles()):
                multipoleParameters = multipoleForce.getMultipoleParameters(atomIndex)
                gkForceBuilder.registerParams(multipoleParameters[0])
            gkForce = gkForceBuilder.getForce(sys, implicitSolvent=implicitSolvent)
            gkForceBuilder.addParticles(gkForce, list(self.topology.atoms()), bonds)

        # Add AmoebaWcaDispersionForce
        if implicitSolvent and "vdw" in self._forces:
            wcaDispersionForceBuilder = AmoebaWcaDispersionForceBuilder(**self.WCA_PARAMS)
            convert = 0.1
            if self._scalars["radiustype"] == "SIGMA":
                convert *= 1.122462048309372
            if self._scalars["radiussize"] == "DIAMETER":
                convert *= 0.5
            # Register VdW parameters for each atom class
            for vdw in self._forces["vdw"]:
                atomClass = int(vdw[0])
                sigma = float(vdw[1])*convert
                epsilon = float(vdw[2])*4.184
                wcaDispersionForceBuilder.registerParams(atomClass, sigma, epsilon)

            wcaDispersionForce = wcaDispersionForceBuilder.getForce(sys)
            atomClasses = [int(atom.atomClass) for atom in self.atoms]
            wcaDispersionForceBuilder.addParticles(wcaDispersionForce, atomClasses)

        # Set periodic boundary conditions
        boxVectors = self.topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
        elif nonbondedMethod not in [ff.NoCutoff, ff.CutoffNonPeriodic]:
            raise ValueError('Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')

        # Adjust masses of hydrogens
        if hydrogenMass is not None:
            for at1, at2 in self.topology.bonds():
                if at1.element == elem.hydrogen:
                    (at1, at2) = (at2, at1)
                if at2.element == elem.hydrogen and at1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(at2.index)
                    sys.setParticleMass(at2.index, hydrogenMass)
                    sys.setParticleMass(at1.index, sys.getParticleMass(at1.index)-transferMass)

        # Add CMMotionRemover
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        return sys


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
    def _processGenericMolecule(molecule: List[int], atoms: List[TinkerAtom]) -> List[Tuple[List[int], str]]:
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
        List[Tuple[List[int], str]]
            A list containing [atomIndices, residueName] for the molecule.
        """
        return [(molecule, "UNK")]

    @staticmethod
    def _processWaterMolecule(molecule: List[int], atoms: List[TinkerAtom]) -> Union[List[Tuple[List[int], str]], bool]:
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
        Union[List[Tuple[List[int], str]], bool]
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

        return [(molecule, "HOH")]

    @staticmethod
    def _processPeptideChain(molecule: List[int], atoms: List[TinkerAtom]) -> List[Tuple[List[int], str]]:
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
        List[Tuple[List[int], str]] 
            If peptide chain is found, returns list of (residueAtoms, residueLabel) tuples.
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
    ) -> List[Tuple[List[int], str]]:
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
        List[Tuple[List[int], str]]
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
    def _processIonMolecule(molecule: List[int], atoms: List[TinkerAtom]) -> Union[List[Tuple[List[int], str]], bool]:
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
        Union[List[Tuple[List[int], str]], bool]
            If an ion is identified, returns a list containing [atomIndices, ionResidueName].
            Otherwise, returns False.
        """
        if len(molecule) != 1:
            return False

        atom = atoms[molecule[0]]
        resName = TinkerFiles._getIonResidueName(atom.atomicNumber)
        return [(molecule, resName)]

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
    def _findBitorsions(atoms: List[TinkerAtom]) -> List[Tuple[int, int, int, int, int]]:
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
    def _addResidueToTopology(topology: top.Topology, resAtoms: List[int], resName: str, chain: top.Chain, atoms: List[TinkerAtom]) -> None:
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
    def _getPeptideResidueAtoms(atoms: List[TinkerAtom], atomIndex: Optional[int] = None) -> List[Tuple[List[int], str]]:
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
            if atomIndex is not None and atomIndex not in [ib, ic, id]:
                continue

            # Check for N/C-terminal cap groups (ACE, NME)
            capAtoms = TinkerFiles._checkCapGroup(ia, ib, ic, id, ie, atoms)
            if capAtoms is not None:
                allAtoms = sorted(list(capAtoms))
                seenAtoms.update(allAtoms)
                residues.append((allAtoms, "ACE" if atoms[ia].atomicNumber == 6 else "NME"))
                continue
      
            if any(atom in seenAtoms for atom in [ib, ic, id]):
                continue

            # Backbone pattern
            if TinkerFiles._checkBackbonePattern(ia, ib, ic, id, ie, atoms):
                residueLabel, sideChain = TinkerFiles._identifyAminoAcid(ic, ib, id, atoms)
                if residueLabel:
                    residueAtoms = TinkerFiles._collectResidueAtoms(ic, ib, id, sideChain, atoms)
                    seenAtoms.update(residueAtoms)
                    residues.append((residueAtoms, residueLabel))
                continue

            # Terminal pattern
            if TinkerFiles._checkTerminalPattern(ia, ib, ic, id, ie, atoms):
                residueLabel, sideChain = TinkerFiles._identifyAminoAcid(ic, ib, id, atoms)
                if residueLabel:
                    residueAtoms = TinkerFiles._collectResidueAtoms(ic, ib, id, sideChain, atoms)
                    seenAtoms.update(residueAtoms)
                    residues.append((residueAtoms, residueLabel))
                continue

        return sorted(residues, key=lambda x: x[0][0])

    @staticmethod
    def _checkBackbonePattern(ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]) -> bool:
        """
        Check if a 5-atom sequence matches a peptide backbone pattern.
        
        Parameters
        ----------
        ia, ib, ic, id, ie : int
            Indices of the 5 atoms in the sequence
        atoms : List[TinkerAtom]
            List of atoms in the system
            
        Returns
        -------
        bool
            True if the sequence matches a peptide backbone pattern, False otherwise
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
    def _checkTerminalPattern(ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]) -> bool:
        """
        Check if a 5-atom sequence matches a terminal residue pattern.
        
        Parameters
        ----------
        ia, ib, ic, id, ie : int
            Indices of the 5 atoms in the sequence
        atoms : List[TinkerAtom]
            List of atoms in the system
            
        Returns
        -------
        bool
            True if the sequence matches a terminal residue pattern, False otherwise
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
    def _checkCapGroup(ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]) -> Optional[List[int]]:
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
        ) or (
            len(atomA.bonds) == 3
            and len(atomB.bonds) == 4
            and len(atomC.bonds) == 3
            and len(atomD.bonds) == 3
            and len(atomE.bonds) == 4
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
            and len(atomD.bonds) == 3
            and len(atomE.bonds) == 4
            and atomA.atomicNumber == 7
            and atomB.atomicNumber == 6
            and atomC.atomicNumber == 6
            and atomD.atomicNumber == 7
            and atomE.atomicNumber == 6
        ) or (
            len(atomA.bonds) == 4
            and len(atomB.bonds) == 3
            and len(atomC.bonds) == 3
            and len(atomD.bonds) == 4
            and len(atomE.bonds) == 3
            and atomA.atomicNumber == 6
            and atomB.atomicNumber == 7
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
    def _identifyAminoAcid(ic: int, ib: int, id: int, atoms: List[TinkerAtom]) -> Tuple[str, Dict[str, int]]:
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
                    label = "CYS" # Deprotonated cysteine, CYS
                else:
                    for bond in atoms[sideChain["isg"]].bonds:
                        if atoms[bond].atomicNumber == 1:
                            label = "CYS"

                        elif atoms[bond].atomicNumber == 16:
                            sideChain["isg"] = 0
                            label = "CYS" # Disulfide bonded cysteine, CYX
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
                    label = "ASP" # Protonated aspartic acid, ASH
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
                    label = "GLU" # Protonated glutamic acid, GLH
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
                        label = "HIS" # Protonated histidine, HIE
                        return label, sideChain
                    if len(atoms[sideChain["ine"]].bonds) == 2:
                        label = "HIS" # Protonated histidine, HID
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
                    label = "LYS" # LYD
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
                        label = "TYR" # "TYD"
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
    def _collectResidueAtoms(ic: int, ib: int, id: int, sideChain: Dict[str, int], atoms: List[TinkerAtom]) -> List[int]:
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
    def _getNucleicAcidResidueAtoms(atoms: List[TinkerAtom]) -> List[Tuple[List[int], str]]:
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
    def _checkNucleotidePattern(ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]) -> bool:
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
    def _identifyNucleotide(ia: int, ib: int, ic: int, id: int, ie: int, atoms: List[TinkerAtom]) -> Tuple[str, Dict[str, int]]:
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
                label = "G" if not deoxy else "DG"
            elif nucleotideAtoms["in9"] != 0:
                label = "A" if not deoxy else "DA"

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
            label = "U"

            # Check if cytosine (has N-H)
            if (
                nucleotideAtoms["in5"] != 0
                and len(atoms[nucleotideAtoms["in5"]].bonds) == 2
            ):
                # Cytosine identified by N5 with 2 bonds
                label = "C" if not deoxy else "DC"

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
                            label = "T" if not deoxy else "DT"
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
    def _loadKeyFile(keyFile: str) -> Tuple[Dict[str, TinkerAtomType], Dict[str, Any], Dict[str, str]]:
        """
        Load a TINKER .key or .prm file.

        Parameters
        ----------
        keyFile : str
            The path to the .key or .prm file.

        Returns
        -------
        Tuple[Dict[str, TinkerAtomType], Dict[str, Any], Dict[str, str]]
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