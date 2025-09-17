"""
amoebaforces.py: AMOEBA force field classes.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2025 Stanford University and the Authors.
Authors: Joao Morado, Peter Eastman
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

__author__ = "Joao Morado"

import math
from collections import defaultdict, namedtuple
import openmm as mm
import openmm.app.forcefield as ff
from openmm.app import element
from openmm.unit import amu, radian

from typing import Optional, List, Dict, Tuple, Union, Any, Callable


class BaseAmoebaForceBuilder:
    """Base class for AMOEBA force builders"""

    def __init__(self) -> None:
        self.params = []

    @staticmethod
    def _getAtomicNumber(atom: Any) -> int:
        """
        Get the atomic number of an atom.

        Parameters
        ----------
        atom : Any
            Atom object that has either an element attribute or atomicNumber attribute.

        Returns
        -------
        int
            The atomic number of the atom.

        Raises
        ------
        ValueError
            If the atom does not have an associated element or atomic number.
        """
        if hasattr(atom, 'element'):
            return atom.element.atomic_number
        elif hasattr(atom, 'atomicNumber'):
            return atom.atomicNumber
        else:
            raise ValueError(f"Atom {atom} does not have an associated element or atomic number.")

    def _findExistingForce(self, sys: mm.System, forceType: type, energyFunction: Optional[str] = None, name: Optional[str] = None) -> Optional[Any]:
        """
        Find existing force in system by type and optionally by energy function or name.

        Parameters
        ----------
        sys : mm.System
            The OpenMM System to search in.
        forceType : type
            The type of force to search for.
        energyFunction : Optional[str], default=None
            The energy function to match (if applicable).
        name : Optional[str], default=None
            The name to match (if applicable).

        Returns
        -------
        Optional[Any]
            The matching force object if found, None otherwise.
        """
        existing = [f for f in sys.getForces() if isinstance(f, forceType)]
        if energyFunction:
            existing = [f for f in existing if hasattr(f, 'getEnergyFunction') and f.getEnergyFunction() == energyFunction]
        if name:
            existing = [f for f in existing if hasattr(f, 'getName') and f.getName() == name]
        return existing[0] if existing else None

    def _createOrGetForce(self, sys: mm.System, forceType: type, creatorFunc: Callable[[], Any], energyFunction: Optional[str] = None, name: Optional[str] = None) -> Any:
        """
        Create new force or get existing one.

        Parameters
        ----------
        sys : mm.System
            The OpenMM System to add the force to.
        forceType : type
            The type of force to create or get.
        creatorFunc : Callable[[], Any]
            Function that creates a new force instance.
        energyFunction : Optional[str], default=None
            The energy function to match for existing forces.
        name : Optional[str], default=None
            The name to match for existing forces.

        Returns
        -------
        Any
            The force object (either existing or newly created).
        """
        existing = self._findExistingForce(sys, forceType, energyFunction, name)
        if existing:
            return existing

        force = creatorFunc()
        sys.addForce(force)
        return force

    def _matchParams(self, atomTypes: Tuple[Any, ...], paramTypes: Tuple[Any, ...], reverseMatch: bool = True) -> bool:
        """
        Match atom types with parameter types, optionally allowing reverse matching.

        Parameters
        ----------
        atomTypes : Tuple[Any, ...]
            Tuple of atom types to match.
        paramTypes : Tuple[Any, ...]
            Tuple of parameter types to match against.
        reverseMatch : bool, default=True
            Whether to allow reverse matching.

        Returns
        -------
        bool
            True if the atom types match the parameter types.
        """
        if tuple(atomTypes) == tuple(paramTypes):
            return True
        if reverseMatch and tuple(atomTypes) == tuple(paramTypes[::-1]):
            return True
        return False

    def _findMatchingParams(self, param_list: List[Tuple[Tuple[Any, ...], Any]], atomTypes: Tuple[Any, ...], reverseMatch: bool = True) -> Optional[Any]:
        """
        Find matching parameters from a parameter list.

        Parameters
        ----------
        param_list : List[Tuple[Tuple[Any, ...], Any]]
            List of (param_type, params) tuples to search through.
        atomTypes : Tuple[Any, ...]
            Tuple of atom types to match.
        reverseMatch : bool, default=True
            Whether to allow reverse matching.

        Returns
        -------
        Optional[Any]
            The matching parameters if found, None otherwise.
        """
        for param_type, params in param_list:
            if self._matchParams(atomTypes, param_type, reverseMatch):
                return params
        return None


class AmoebaBondForceBuilder(BaseAmoebaForceBuilder):
    """Builder for Bond force for AMOEBA force field"""

    def __init__(self, cubic, quartic):
        super().__init__()
        self.cubic = cubic
        self.quartic = quartic
        self.bondParams = []

    def registerParams(self, bondType, r0, k0):
        """Register bond parameters"""
        self.bondParams.append((bondType, (r0, k0)))

    def getForce(self, sys):
        energy = f"k*(d^2 + {self.cubic}*d^3 + {self.quartic}*d^4); d=r-r0"

        def createForce():
            force = mm.CustomBondForce(energy)
            force.addPerBondParameter("r0")
            force.addPerBondParameter("k")
            force.setName("AmoebaBond")
            return force

        return self._createOrGetForce(sys, mm.CustomBondForce, createForce, energyFunction=energy)

    def addBonds(self, force, atomClasses, bonds):
        """Add bonds to the force"""
        for atom1, atom2 in bonds:
            atomTypes = (atomClasses[atom1], atomClasses[atom2])
            params = self._findMatchingParams(self.bondParams, atomTypes)
            if params is not None and params[1] != 0:
                force.addBond(atom1, atom2, params)

class AmoebaAngleForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for AngleForce force for AMOEBA force field.

    Attributes
    ----------
    cubic : float
        Cubic term coefficient.
    quartic : float
        Quartic term coefficient.
    pentic : float
        Pentic (fifth-order) term coefficient.
    sextic : float
        Sextic (sixth-order) term coefficient.
    angleParams : List[Tuple[Tuple[Any, ...], Any]]
        List of angle parameters as tuples of (angleType, params).
    """

    def __init__(self, cubic: float, quartic: float, pentic: float, sextic: float) -> None:
        """
        Initialize the AmoebaAngleForceBuilder.

        Parameters
        ----------
        cubic : float
            Cubic term coefficient.
        quartic : float
            Quartic term coefficient.
        pentic : float
            Pentic (fifth-order) term coefficient.
        sextic : float
            Sextic (sixth-order) term coefficient.
        """
        super().__init__()
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic
        self.angleParams = []

    def registerParams(self, angleType: Tuple[Any, ...], params: Tuple[float, ...]) -> None:
        """
        Register angle parameters.

        Parameters
        ----------
        angleType : Tuple[Any, ...]
            A tuple of atom classes defining the angle type.
        params : Tuple[float, ...]
            The angle parameters (equilibrium angle, force constant, etc.).
        """
        self.angleParams.append((angleType, params))

    def getForce(self, sys: mm.System) -> mm.CustomAngleForce:
        """
        Get or create the AmoebaAngleForce in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.

        Returns
        -------
        mm.CustomAngleForce
            The AmoebaAngleForce instance.
        """
        energy = "k*(d^2 + %s*d^3 + %s*d^4 + %s*d^5 + %s*d^6); d=%.15g*theta-theta0" % (
            self.cubic,
            self.quartic,
            self.pentic,
            self.sextic,
            180 / math.pi,
        )

        def createForce():
            force = mm.CustomAngleForce(energy)
            force.addPerAngleParameter("theta0")
            force.addPerAngleParameter("k")
            force.setName("AmoebaAngle")
            return force

        return self._createOrGetForce(sys, mm.CustomAngleForce, createForce, energyFunction=energy)

    def addAngles(self, force: mm.CustomAngleForce, angles: List[Tuple[int, int, int]]) -> None:
        """
        Add angles to the force.

        Parameters
        ----------
        force : mm.CustomAngleForce
            The AmoebaAngleForce instance to add angles to.
        angles : List[Tuple[int, int, int]]
            List of angle indices as tuples of (atom1, atom2, atom3).
        """
        for atom1, atom2, atom3 in angles:
            for angleType, params in self.angleParams:
                if params[1] != 0:
                    if self._matchParams((atom1, atom2, atom3), angleType):
                        force.addAngle(atom1, atom2, atom3, params)
                        break


class AmoebaInPlaneAngleForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for InPlaneAngleForce force for AMOEBA force field.


    Attributes
    ----------
    cubic : float
        Cubic term coefficient.
    quartic : float
        Quartic term coefficient.
    pentic : float
        Pentic (fifth-order) term coefficient.
    sextic : float
        Sextic (sixth-order) term coefficient.
    inPlaneAngleParams : List[Tuple[Tuple[Any, ...], Any]]
        List of in-plane angle parameters.
    """

    def __init__(self, cubic: float, quartic: float, pentic: float, sextic: float) -> None:
        """
        Initialize the AmoebaInPlaneAngleForceBuilder.

        Parameters
        ----------
        cubic : float
            Cubic term coefficient.
        quartic : float
            Quartic term coefficient.
        pentic : float
            Pentic (fifth-order) term coefficient.
        sextic : float
            Sextic (sixth-order) term coefficient.
        """
        super().__init__()
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic
        self.inPlaneAngleParams = []

    def registerParams(self, inPlaneAngleType: Tuple[Any, ...], params: Tuple[float, ...]) -> None:
        """
        Register in-plane angle parameters.

        Parameters
        ----------
        inPlaneAngleType : Tuple[Any, ...]
            A tuple of atom classes defining the in-plane angle type.
        params : Tuple[float, ...]
            The in-plane angle parameters (equilibrium angle, force constant, etc.).
        """
        self.inPlaneAngleParams.append((inPlaneAngleType, params))

    def getForce(self, sys: mm.System) -> mm.CustomCompoundBondForce:
        """
        Get or create the AmoebaInPlaneAngleForce in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.

        Returns
        -------
        mm.CustomCompoundBondForce
            The AmoebaInPlaneAngleForce instance.
        """
        energy = """k*(d^2 + %s*d^3 + %s*d^4 + %s*d^5 + %s*d^6); d=theta-theta0;
            theta = %.15g*pointangle(x1, y1, z1, projx, projy, projz, x3, y3, z3);
            projx = x2-nx*dot; projy = y2-ny*dot; projz = z2-nz*dot;
            dot = nx*(x2-x3) + ny*(y2-y3) + nz*(z2-z3);
            nx = px/norm; ny = py/norm; nz = pz/norm;
            norm = sqrt(px*px + py*py + pz*pz);
            px = (d1y*d2z-d1z*d2y); py = (d1z*d2x-d1x*d2z); pz = (d1x*d2y-d1y*d2x);
            d1x = x1-x4; d1y = y1-y4; d1z = z1-z4;
            d2x = x3-x4; d2y = y3-y4; d2z = z3-z4""" % (
            self.cubic,
            self.quartic,
            self.pentic,
            self.sextic,
            180.0 / math.pi,
        )

        def createForce():
            force = mm.CustomCompoundBondForce(4, energy)
            force.addPerBondParameter("theta0")
            force.addPerBondParameter("k")
            force.setName("AmoebaInPlaneAngle")
            return force

        return self._createOrGetForce(sys, mm.CustomCompoundBondForce, createForce, energyFunction=energy)
    
    def addInPlaneAngles(self, force: mm.CustomCompoundBondForce, atomClasses: List[Any], inPlaneAngles: List[Tuple[int, int, int, int]]) -> None:
        """
        Add in-plane angles to the force.

        Parameters
        ----------
        force : mm.CustomCompoundBondForce
            The AmoebaInPlaneAngleForce instance to add angles to.
        atomClasses : List[Any]
            List of atom classes indexed by atom index.
        inPlaneAngles : List[Tuple[int, int, int, int]]
            List of in-plane angle indices as tuples of (atom1, atom2, atom3, atom4).
        """
        for atom1, atom2, atom3, atom4 in inPlaneAngles:
            for inPlaneAngleType, params in self.inPlaneAngleParams:
                if params[1] != 0:
                    angleClasses = (atomClasses[atom1], atomClasses[atom2], atomClasses[atom3])
                    if self._matchParams(angleClasses, inPlaneAngleType[:3], reverseMatch=True):
                        force.addBond([atom1, atom2, atom3, atom4], params)
                        break


class AmoebaOutOfPlaneBendForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for OutOfPlaneBendForce force for AMOEBA force field.

    Attributes
    ----------
    cubic : float
        Cubic term coefficient.
    quartic : float
        Quartic term coefficient.
    pentic : float
        Pentic (fifth-order) term coefficient.
    sextic : float
        Sextic (sixth-order) term coefficient.
    outOfPlaneBendParams : List[Tuple[Tuple[Any, ...], Any]]
        List of out-of-plane bend parameters.
    """

    def __init__(self, cubic: float, quartic: float, pentic: float, sextic: float) -> None:
        """
        Initialize the AmoebaOutOfPlaneBendForceBuilder.

        Parameters
        ----------
        cubic : float
            Cubic term coefficient.
        quartic : float
            Quartic term coefficient.
        pentic : float
            Pentic (fifth-order) term coefficient.
        sextic : float
            Sextic (sixth-order) term coefficient.
        """
        super().__init__()
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic
        self.outOfPlaneBendParams = []

    def _matchParams(self, atomClasses: Tuple[Any, ...], paramClasses: Tuple[Any, ...]) -> bool:
        """
        Match atom classes with parameter classes for out-of-plane bends.

        Parameters
        ----------
        atomClasses : Tuple[Any, ...]
            Tuple of atom classes (at1, at2, at3, at4).
        paramClasses : Tuple[Any, ...]
            Tuple of parameter classes (p1, p2, p3, p4).

        Returns
        -------
        bool
            True if the atom classes match the parameter classes.
        """
        at1, at2, at3, at4 = atomClasses
        p1, p2, p3, p4 = paramClasses
        return  at2 == p2 and at4 == p1 and {at1, at3} == {at1 if p3 == "0" else p3, at3 if p4 == "0" else p4}
    
    def registerParams(self, outOfPlaneBendType: Tuple[Any, ...], params: Tuple[float, ...]) -> None:
        """
        Register out-of-plane bend parameters.

        Parameters
        ----------
        outOfPlaneBendType : Tuple[Any, ...]
            A tuple of atom classes defining the out-of-plane bend type.
        params : Tuple[float, ...]
            The out-of-plane bend parameters (force constant, etc.).
        """
        self.outOfPlaneBendParams.append((outOfPlaneBendType, params))

    def getForce(self, sys: mm.System) -> mm.CustomCompoundBondForce:
        """
        Get or create the AmoebaOutOfPlaneBendForce in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.

        Returns
        -------
        mm.CustomCompoundBondForce
            The AmoebaOutOfPlaneBendForce instance.
        """
        energy = """k*(theta^2 + %s*theta^3 + %s*theta^4 + %s*theta^5 + %s*theta^6);
                    theta = %.15g*pointangle(x2, y2, z2, x4, y4, z4, projx, projy, projz);
                    projx = x2-nx*dot; projy = y2-ny*dot; projz = z2-nz*dot;
                    dot = nx*(x2-x3) + ny*(y2-y3) + nz*(z2-z3);
                    nx = px/norm; ny = py/norm; nz = pz/norm;
                    norm = sqrt(px*px + py*py + pz*pz);
                    px = (d1y*d2z-d1z*d2y); py = (d1z*d2x-d1x*d2z); pz = (d1x*d2y-d1y*d2x);
                    d1x = x1-x4; d1y = y1-y4; d1z = z1-z4;
                    d2x = x3-x4; d2y = y3-y4; d2z = z3-z4""" % (
            self.cubic,
            self.quartic,
            self.pentic,
            self.sextic,
            180.0 / math.pi,
        )

        def createForce():
            force = mm.CustomCompoundBondForce(4, energy)
            force.addPerBondParameter("k")
            force.setName("AmoebaOutOfPlaneBend")
            return force

        return self._createOrGetForce(sys, mm.CustomCompoundBondForce, createForce, energyFunction=energy)

    def addOutOfPlaneBends(self, force: mm.CustomCompoundBondForce, atomClasses: List[Any], outOfPlaneBends: List[Tuple[int, int, int, int]]) -> None:
        """
        Add out-of-plane bends to the force.

        Parameters
        ----------
        force : mm.CustomCompoundBondForce
            The AmoebaOutOfPlaneBendForce instance to add bends to.
        atomClasses : List[Any]
            List of atom classes indexed by atom index.
        outOfPlaneBends : List[Tuple[int, int, int, int]]
            List of out-of-plane bend indices as tuples of (atom1, atom2, atom3, atom4).
        """
        for atom1, atom2, atom3, atom4 in outOfPlaneBends:
            for outOfPlaneBendType, params in self.outOfPlaneBendParams:
                if params[0] != 0:
                    angleClasses = (atomClasses[atom1], atomClasses[atom2], atomClasses[atom3], atomClasses[atom4])
                    if self._matchParams(angleClasses, outOfPlaneBendType):
                        force.addBond([atom1, atom2, atom3, atom4], params)
                        break


class AmoebaStretchBendForceBuilder(BaseAmoebaForceBuilder):
    """Builder for StretchBendForce force for AMOEBA force field"""

    def __init__(self):
        super().__init__()
        self.stretchBendParams = []

    def registerParams(self, stretchBendType, params):
        """Register stretch-bend parameters"""
        self.stretchBendParams.append((stretchBendType, params))

    def getForce(self, sys):
        energy = (
            "(k1*(distance(p1,p2)-r12) + k2*(distance(p2,p3)-r23))*(%.15g*(angle(p1,p2,p3)-theta0))"
            % (180 / math.pi)
        )

        def createForce():
            force = mm.CustomCompoundBondForce(3, energy)
            force.addPerBondParameter("r12")
            force.addPerBondParameter("r23")
            force.addPerBondParameter("theta0")
            force.addPerBondParameter("k1")
            force.addPerBondParameter("k2")
            force.setName("AmoebaStretchBend")
            return force

        return self._createOrGetForce(sys, mm.CustomCompoundBondForce, createForce, energyFunction=energy)

    def addStretchBends(self, force, angles):
        """Add stretch-bend terms to the force"""
        for atom1, atom2, atom3 in angles:
            for stretchBendType, params in self.stretchBendParams:
                angle = (atom1, atom2, atom3)
                if self._matchParams(angle, stretchBendType, reverseMatch=True):
                    force.addBond(angle, params)
                    break

class AmoebaStretchTorsionForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for StretchTorsionForce force for AMOEBA force field.

    Attributes
    ----------
    stretchTorsionParams : List[Tuple[Tuple[Any, ...], Any]]
        List of stretch-torsion parameters.
    """

    def __init__(self) -> None:
        """Initialize the AmoebaStretchTorsionForceBuilder."""
        super().__init__()
        self.stretchTorsionParams = []

    def registerParams(self, stretchTorsionType: Tuple[Any, ...], params: Tuple[float, ...]) -> None:
        """
        Register stretch-torsion parameters.

        Parameters
        ----------
        stretchTorsionType : Tuple[Any, ...]
            A tuple of atom classes defining the stretch-torsion type.
        params : Tuple[float, ...]
            The stretch-torsion parameters.
        """
        self.stretchTorsionParams.append((stretchTorsionType, params))

    def getForce(self, sys: mm.System) -> mm.CustomCompoundBondForce:
        """
        Get or create the AmoebaStretchTorsionForce in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.

        Returns
        -------
        mm.CustomCompoundBondForce
            The AmoebaStretchTorsionForce instance.
        """
        energy = """v11*(distance(p1,p2)-length1)*phi1 +
                    v12*(distance(p1,p2)-length1)*phi2 +
                    v13*(distance(p1,p2)-length1)*phi3 +
                    v21*(distance(p2,p3)-length2)*phi1 +
                    v22*(distance(p2,p3)-length2)*phi2 +
                    v23*(distance(p2,p3)-length2)*phi3 +
                    v31*(distance(p3,p4)-length3)*phi1 +
                    v32*(distance(p3,p4)-length3)*phi2 +
                    v33*(distance(p3,p4)-length3)*phi3;
                    phi1=1+cos(phi+phase1); phi2=1+cos(2*phi+phase2); phi3=1+cos(3*phi+phase3);
                    phi=dihedral(p1,p2,p3,p4)"""

        def createForce():
            force = mm.CustomCompoundBondForce(4, energy)
            for param in ('v11', 'v12', 'v13', 'v21', 'v22', 'v23', 'v31', 'v32', 'v33'):
                force.addPerBondParameter(param)
            for i in range(3):
                force.addPerBondParameter(f'length{i+1}')
            for i in range(3):
                force.addPerBondParameter(f'phase{i+1}')
            force.setName('AmoebaStretchTorsion')
            return force

        return self._createOrGetForce(sys, mm.CustomCompoundBondForce, createForce, energyFunction=energy)

    def addStretchTorsions(self, sys: mm.System, force: mm.CustomCompoundBondForce, atomClasses: List[Any], torsions: List[Tuple[int, int, int, int]]) -> None:
        """
        Add stretch-torsion terms to the force.

        Parameters
        ----------
        sys : mm.System
            The OpenMM System containing bond and torsion forces.
        force : mm.CustomCompoundBondForce
            The AmoebaStretchTorsionForce instance to add terms to.
        atomClasses : List[Any]
            List of atom classes indexed by atom index.
        torsions : List[Tuple[int, int, int, int]]
            List of torsion indices as tuples of (atom1, atom2, atom3, atom4).
        """

        # Record parameters for bonds and torsions so we can look them up quickly.

        bondForce = [f for f in sys.getForces() if type(f) == mm.CustomBondForce and f.getName() == 'AmoebaBond'][0]
        torsionForce = [f for f in sys.getForces() if type(f) == mm.PeriodicTorsionForce][0]
        bondLength = {}
        torsionPhase = defaultdict(lambda: [0.0, math.pi, 0.0])
        for i in range(bondForce.getNumBonds()):
            p1, p2, params = bondForce.getBondParameters(i)
            bondLength[(p1, p2)] = params[0]
            bondLength[(p2, p1)] = params[0]
        for i in range(torsionForce.getNumTorsions()):
            p1, p2, p3, p4, periodicity, phase, k = torsionForce.getTorsionParameters(i)
            if periodicity < 4:
                phase = phase.value_in_unit(radian)
                torsionPhase[(p1, p2, p3, p4)][periodicity-1] = phase
                torsionPhase[(p4, p3, p2, p1)][periodicity-1] = phase

        # Add stretch-torsions.

        for torsion in torsions:
            atom1, atom2, atom3, atom4 = torsion
            atomTypes = (atomClasses[atom1], atomClasses[atom2], atomClasses[atom3], atomClasses[atom4])
            for torsionType, params in self.stretchTorsionParams:
                if self._matchParams(atomTypes, torsionType):
                    if atomTypes == tuple(reversed(torsionType)):
                        atom1, atom2, atom3, atom4 = atom4, atom3, atom2, atom1
                    params = list(params)
                    params.append(bondLength[(atom1, atom2)])
                    params.append(bondLength[(atom2, atom3)])
                    params.append(bondLength[(atom3, atom4)])
                    params += torsionPhase[torsion]
                    force.addBond((atom1, atom2, atom3, atom4), params)
                    break


class AmoebaAngleTorsionForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for AngleTorsionForce force for AMOEBA force field.

    Attributes
    ----------
    angleTorsionParams : List[Tuple[Tuple[Any, ...], Any]]
        List of angle-torsion parameters.
    """

    def __init__(self) -> None:
        """Initialize the AmoebaAngleTorsionForceBuilder."""
        super().__init__()
        self.angleTorsionParams = []

    def registerParams(self, angleTorsionType: Tuple[Any, ...], params: Tuple[float, ...]) -> None:
        """
        Register angle-torsion parameters.

        Parameters
        ----------
        angleTorsionType : Tuple[Any, ...]
            A tuple of atom classes defining the angle-torsion type.
        params : Tuple[float, ...]
            The angle-torsion parameters.
        """
        self.angleTorsionParams.append((angleTorsionType, params))

    def getForce(self, sys: mm.System) -> mm.CustomCompoundBondForce:
        """
        Get or create the AmoebaAngleTorsionForce in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.

        Returns
        -------
        mm.CustomCompoundBondForce
            The AmoebaAngleTorsionForce instance.
        """
        energy = """v11*(angle(p1,p2,p3)-angle1)*phi1 +
                    v12*(angle(p1,p2,p3)-angle1)*phi2 +
                    v13*(angle(p1,p2,p3)-angle1)*phi3 +
                    v21*(angle(p2,p3,p4)-angle2)*phi1 +
                    v22*(angle(p2,p3,p4)-angle2)*phi2 +
                    v23*(angle(p2,p3,p4)-angle2)*phi3;
                    phi1=1+cos(phi+phase1); phi2=1+cos(2*phi+phase2); phi3=1+cos(3*phi+phase3);
                    phi=dihedral(p1,p2,p3,p4)"""

        def createForce():
            force = mm.CustomCompoundBondForce(4, energy)
            for param in ('v11', 'v12', 'v13', 'v21', 'v22', 'v23'):
                force.addPerBondParameter(param)
            for i in range(2):
                force.addPerBondParameter(f'angle{i+1}')
            for i in range(3):
                force.addPerBondParameter(f'phase{i+1}')
            force.setName('AmoebaAngleTorsion')
            return force

        return self._createOrGetForce(sys, mm.CustomCompoundBondForce, createForce, energyFunction=energy)

    def addAngleTorsions(self, sys: mm.System, force: mm.CustomCompoundBondForce, atomClasses: List[Any], torsions: List[Tuple[int, int, int, int]]) -> None:
        """
        Add angle-torsion terms to the force.

        Parameters
        ----------
        sys : mm.System
            The OpenMM System containing angle and torsion forces.
        force : mm.CustomCompoundBondForce
            The AmoebaAngleTorsionForce instance to add terms to.
        atomClasses : List[Any]
            List of atom classes indexed by atom index.
        torsions : List[Tuple[int, int, int, int]]
            List of torsion indices as tuples of (atom1, atom2, atom3, atom4).
        """

        # Record parameters for angles and torsions so we can look them up quickly.

        angleForce = [f for f in sys.getForces() if type(f) == mm.CustomAngleForce and f.getName() == 'AmoebaAngle'][0]
        inPlaneAngleForce = [f for f in sys.getForces() if type(f) == mm.CustomCompoundBondForce and f.getName() == 'AmoebaInPlaneAngle'][0]
        torsionForce = [f for f in sys.getForces() if type(f) == mm.PeriodicTorsionForce][0]
        equilAngle = {}
        torsionPhase = defaultdict(lambda: [0.0, math.pi, 0.0])
        angleScale = math.pi/180
        for i in range(angleForce.getNumAngles()):
            p1, p2, p3, params = angleForce.getAngleParameters(i)
            equilAngle[(p1, p2, p3)] = params[0]*angleScale
            equilAngle[(p3, p2, p1)] = params[0]*angleScale
        for i in range(inPlaneAngleForce.getNumBonds()):
            particles, params = inPlaneAngleForce.getBondParameters(i)
            equilAngle[tuple(particles[:3])] = params[0]*angleScale
            equilAngle[tuple(reversed(particles[:3]))] = params[0]*angleScale
        for i in range(torsionForce.getNumTorsions()):
            p1, p2, p3, p4, periodicity, phase, k = torsionForce.getTorsionParameters(i)
            if periodicity < 4:
                phase = phase.value_in_unit(radian)
                torsionPhase[(p1, p2, p3, p4)][periodicity-1] = phase
                torsionPhase[(p4, p3, p2, p1)][periodicity-1] = phase

        # Add angle-torsions.

        for torsion in torsions:
            atom1, atom2, atom3, atom4 = torsion
            atomTypes = (atomClasses[atom1], atomClasses[atom2], atomClasses[atom3], atomClasses[atom4])
            for torsionType, params in self.angleTorsionParams:
                if self._matchParams(atomTypes, torsionType):
                    if atomTypes == tuple(reversed(torsionType)):
                        atom1, atom2, atom3, atom4 = atom4, atom3, atom2, atom1
                    params = list(params)
                    params.append(equilAngle[(atom1, atom2, atom3)])
                    params.append(equilAngle[(atom2, atom3, atom4)])
                    params += torsionPhase[torsion]
                    force.addBond((atom1, atom2, atom3, atom4), params)
                    break


class AmoebaTorsionForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for PeriodicTorsionForce force for AMOEBA force field.

    Attributes
    ----------
    torsionParams : List[Tuple[Tuple[Any, ...], Any]]
        List of torsion parameters.
    """

    def __init__(self) -> None:
        """Initialize the AmoebaTorsionForceBuilder."""
        super().__init__()
        self.torsionParams = []

    def registerParams(self, torsionType: Tuple[Any, ...], params: Tuple[Tuple[float, float], ...]) -> None:
        """
        Register torsion parameters.

        Parameters
        ----------
        torsionType : Tuple[Any, ...]
            A tuple of atom classes defining the torsion type.
        params : Tuple[Tuple[float, float], ...]
            The torsion parameters as tuples of (force_constant, phase_angle) for each periodicity.
        """
        self.torsionParams.append((torsionType, params))

    def getForce(self, sys: mm.System) -> mm.PeriodicTorsionForce:
        """
        Get or create the PeriodicTorsionForce in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.

        Returns
        -------
        mm.PeriodicTorsionForce
            The PeriodicTorsionForce instance.
        """
        return self._createOrGetForce(sys, mm.PeriodicTorsionForce, mm.PeriodicTorsionForce)

    def addTorsions(self, force: mm.PeriodicTorsionForce, atomClasses: List[Any], torsions: List[Tuple[int, int, int, int]]) -> None:
        """
        Add torsions to the force.

        Parameters
        ----------
        force : mm.PeriodicTorsionForce
            The PeriodicTorsionForce instance to add torsions to.
        atomClasses : List[Any]
            List of atom classes indexed by atom index.
        torsions : List[Tuple[int, int, int, int]]
            List of torsion indices as tuples of (atom1, atom2, atom3, atom4).
        """
        for atom1, atom2, atom3, atom4 in torsions:
            for torsionType, params in self.torsionParams:
                atomTypes = (atomClasses[atom1], atomClasses[atom2],
                             atomClasses[atom3], atomClasses[atom4])
                if self._matchParams(atomTypes, torsionType):
                    t1, t2, t3 = params
                    if t1[0] != 0:
                        force.addTorsion(atom1, atom2, atom3, atom4, 1, t1[1], t1[0])
                    if t2[0] != 0:
                        force.addTorsion(atom1, atom2, atom3, atom4, 2, t2[1], t2[0])
                    if t3[0] != 0:
                        force.addTorsion(atom1, atom2, atom3, atom4, 3, t3[1], t3[0])
                    break


class AmoebaPiTorsionForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for PiTorsionForce force for AMOEBA force field.

    Attributes
    ----------
    piTorsionParams : List[Tuple[Tuple[Any, ...], Any]]
        List of pi-torsion parameters.
    """

    def __init__(self) -> None:
        """Initialize the AmoebaPiTorsionForceBuilder."""
        super().__init__()
        self.piTorsionParams = []

    def registerParams(self, piTorsionType: Tuple[Any, ...], params: Tuple[float, ...]) -> None:
        """
        Register pi-torsion parameters.

        Parameters
        ----------
        piTorsionType : Tuple[Any, ...]
            A tuple of atom classes defining the pi-torsion type.
        params : Tuple[float, ...]
            The pi-torsion parameters (force constant, etc.).
        """
        self.piTorsionParams.append((piTorsionType, params))

    def getForce(self, sys: mm.System) -> mm.CustomCompoundBondForce:
        """
        Get or create the AmoebaPiTorsionForce in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.

        Returns
        -------
        mm.CustomCompoundBondForce
            The AmoebaPiTorsionForce instance.
        """
        energy = """2*k*sin(phi)^2;
                    phi = pointdihedral(x3+c1x, y3+c1y, z3+c1z, x3, y3, z3, x4, y4, z4, x4+c2x, y4+c2y, z4+c2z);
                    c1x = (d14y*d24z-d14z*d24y); c1y = (d14z*d24x-d14x*d24z); c1z = (d14x*d24y-d14y*d24x);
                    c2x = (d53y*d63z-d53z*d63y); c2y = (d53z*d63x-d53x*d63z); c2z = (d53x*d63y-d53y*d63x);
                    d14x = x1-x4; d14y = y1-y4; d14z = z1-z4;
                    d24x = x2-x4; d24y = y2-y4; d24z = z2-z4;
                    d53x = x5-x3; d53y = y5-y3; d53z = z5-z3;
                    d63x = x6-x3; d63y = y6-y3; d63z = z6-z3"""

        def createForce():
            force = mm.CustomCompoundBondForce(6, energy)
            force.addPerBondParameter("k")
            force.setName("AmoebaPiTorsion")
            return force

        return self._createOrGetForce(sys, mm.CustomCompoundBondForce, createForce, energyFunction=energy)

    def addPiTorsions(self, force: mm.CustomCompoundBondForce, atomClasses: List[Any], piTorsions: List[Tuple[int, int, int, int, int, int]]) -> None:
        """
        Add pi-torsions to the force.

        Parameters
        ----------
        force : mm.CustomCompoundBondForce
            The AmoebaPiTorsionForce instance to add pi-torsions to.
        atomClasses : List[Any]
            List of atom classes indexed by atom index.
        piTorsions : List[Tuple[int, int, int, int, int, int]]
            List of pi-torsion indices as tuples of six atoms.
        """
        for atom1, atom2, atom3, atom4, atom5, atom6 in piTorsions:
            for piTorsionType, params in self.piTorsionParams:
                if params[0] != 0:
                    types = (atomClasses[atom3], atomClasses[atom4])
                    if self._matchParams(types, piTorsionType, reverseMatch=True):
                        force.addBond([atom1, atom2, atom3, atom4, atom5, atom6], params)
                        break


class AmoebaUreyBradleyForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for UreyBradleyForce force for AMOEBA force field.

    Attributes
    ----------
    ureyBradleyParams : List[Tuple[Tuple[Any, ...], Any]]
        List of Urey-Bradley parameters.
    """

    def __init__(self) -> None:
        """Initialize the AmoebaUreyBradleyForceBuilder."""
        super().__init__()
        self.ureyBradleyParams = []

    def registerParams(self, ureyBradleyType: Tuple[Any, ...], params: Tuple[float, float]) -> None:
        """
        Register Urey-Bradley parameters.

        Parameters
        ----------
        ureyBradleyType : Tuple[Any, ...]
            A tuple of atom classes defining the Urey-Bradley type.
        params : Tuple[float, float]
            The Urey-Bradley parameters (force constant, equilibrium distance).
        """
        self.ureyBradleyParams.append((ureyBradleyType, params))

    def getForce(self, sys: mm.System) -> mm.HarmonicBondForce:
        """
        Get or create the HarmonicBondForce for Urey-Bradley terms in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.

        Returns
        -------
        mm.HarmonicBondForce
            The HarmonicBondForce instance for Urey-Bradley terms.
        """
        return self._createOrGetForce(sys, mm.HarmonicBondForce, mm.HarmonicBondForce)

    def addUreyBradleys(self, force: mm.HarmonicBondForce, atomClasses: List[Any], angles: List[Tuple[int, int, int]], anglesConstraints: Optional[List[bool]] = None, flexibleConstraints: bool = False) -> None:
        """
        Add Urey-Bradley terms to the force.

        Parameters
        ----------
        force : mm.HarmonicBondForce
            The HarmonicBondForce instance to add Urey-Bradley terms to.
        atomClasses : List[Any]
            List of atom classes indexed by atom index.
        angles : List[Tuple[int, int, int]]
            List of angle indices as tuples of (atom1, atom2, atom3).
        anglesConstraints : Optional[List[bool]], default=None
            List of flags indicating if a given angle is constrained.
        flexibleConstraints : bool, default=False
            If True, constrained angles will still be added to the system.
        """
        for i, (atom1, atom2, atom3) in enumerate(angles):
            isConstrained = False if anglesConstraints is None else anglesConstraints[i]
            if not isConstrained or flexibleConstraints:
                ubType = (atomClasses[atom1], atomClasses[atom2], atomClasses[atom3])
                params = self._findMatchingParams(self.ureyBradleyParams, ubType, reverseMatch=True)
                if params is not None:
                    k, d = params
                    force.addBond(atom1, atom3, d, 2 * k)


class AmoebaTorsionTorsionForceBuilder(BaseAmoebaForceBuilder):
    """Builder for TorsionTorsion force for AMOEBA force field"""

    def __init__(self):
        super().__init__()
        self.torsionTorsionParams = []
        self.gridData = {}

    def registerParams(self, torsionTorsionType, params):
        """Register torsion-torsion parameters"""
        self.torsionTorsionParams.append((torsionTorsionType, params))

    def registerGridData(self, gridIndex, grid):
        """Register grid data for torsion-torsion interactions"""
        self.gridData[gridIndex] = grid

    def getForce(self, sys):
        def createForce():
            force = mm.AmoebaTorsionTorsionForce()
            return force

        force = self._createOrGetForce(sys, mm.AmoebaTorsionTorsionForce, createForce)

        for gridIndex, grid in self.gridData.items():
            force.setTorsionTorsionGrid(gridIndex, grid)

        return force

    @staticmethod
    def setTorsionTorsionGrid(force, gridIndex, grid):
        """Set the torsion-torsion grid for the force."""
        force.setTorsionTorsionGrid(gridIndex, grid)

    @staticmethod
    def addTorsionTorsion(force, atom1, atom2, atom3, atom4, atom5, chiralIndex, gridIndex):
        """Add a torsion-torsion interaction to the force."""
        force.addTorsionTorsion(atom1, atom2, atom3, atom4, atom5, chiralIndex, gridIndex)

    @staticmethod
    def addTorsionTorsionInteractions(force, data, types1, types2, types3, types4, types5, gridIndex, sys):
        """Add torsion-torsion interactions based on TINKER subroutine bitors()"""
        for angle in data.angles:
            # search for bitorsions; based on TINKER subroutine bitors()
            ib = angle[0]
            ic = angle[1]
            id = angle[2]

            for bondIndex in data.atomBonds[ib]:
                bondedAtom1 = data.bonds[bondIndex].atom1
                bondedAtom2 = data.bonds[bondIndex].atom2
                if (bondedAtom1 != ib):
                    ia = bondedAtom1
                else:
                    ia = bondedAtom2

                if (ia != ic and ia != id):
                    for bondIndex2 in data.atomBonds[id]:
                        bondedAtom1 = data.bonds[bondIndex2].atom1
                        bondedAtom2 = data.bonds[bondIndex2].atom2
                        if (bondedAtom1 != id):
                            ie = bondedAtom1
                        else:
                            ie = bondedAtom2

                        if (ie != ic and ie != ib and ie != ia):
                            # found candidate set of atoms
                            # check if types match in order or reverse order
                            type1 = data.atomType[data.atoms[ia]]
                            type2 = data.atomType[data.atoms[ib]]
                            type3 = data.atomType[data.atoms[ic]]
                            type4 = data.atomType[data.atoms[id]]
                            type5 = data.atomType[data.atoms[ie]]

                            for i in range(len(types1)):
                                types1_i = types1[i]
                                types2_i = types2[i]
                                types3_i = types3[i]
                                types4_i = types4[i]
                                types5_i = types5[i]

                                # match in order
                                if (type1 in types1_i and type2 in types2_i and type3 in types3_i and
                                    type4 in types4_i and type5 in types5_i):
                                    chiralAtomIndex = AmoebaTorsionTorsionForceBuilder._getChiralAtomIndex(data, ib, ic, id)
                                    force.addTorsionTorsion(ia, ib, ic, id, ie, chiralAtomIndex, gridIndex[i])

                                # match in reverse order
                                elif (type5 in types1_i and type4 in types2_i and type3 in types3_i and
                                      type2 in types4_i and type1 in types5_i):
                                    chiralAtomIndex = AmoebaTorsionTorsionForceBuilder._getChiralAtomIndex(data, ib, ic, id)
                                    force.addTorsionTorsion(ie, id, ic, ib, ia, chiralAtomIndex, gridIndex[i])


    @staticmethod
    def createTorsionTorsionInteractions(force, angles, atoms, tortorParams):
        """
        Create torsion-torsion interactions based on the angles and tortor parameters.

        Parameters
        ----------
        force : mm.AmoebaTorsionTorsionForce
            The OpenMM AmoebaTorsionTorsionForce object
        angles : list
            List of angles in the system
        atoms : list
            List of TinkerAtom objects
        tortorParams : list
            List of torsion-torsion parameter sets
        """
        for angle in angles:
            # angle = (atom1, atom2, atom3) where atom2 is the central atom
            ib = angle[0]  # first atom of angle
            ic = angle[1]  # central atom of angle
            id = angle[2]  # last atom of angle

            # Find atoms bonded to ib (excluding ic and id)
            for ia in atoms[ib].bonds:
                if ia != ic and ia != id:
                    # Find atoms bonded to id (excluding ic and ib and ia)
                    for ie in atoms[id].bonds:
                        if ie != ic and ie != ib and ie != ia:
                            # We have a potential torsion-torsion pattern: ia-ib-ic-id-ie
                            # Check if atom types match any tortor parameters
                            for gridIndex, (tortorInfo, gridData) in enumerate(tortorParams):
                                class1_param = tortorInfo[0]
                                class2_param = tortorInfo[1]
                                class3_param = tortorInfo[2]
                                class4_param = tortorInfo[3]
                                class5_param = tortorInfo[4]

                                type1 = atoms[ia].atomClass
                                type2 = atoms[ib].atomClass
                                type3 = atoms[ic].atomClass
                                type4 = atoms[id].atomClass
                                type5 = atoms[ie].atomClass

                                # Check forward direction
                                if (type1 == class1_param and type2 == class2_param and
                                    type3 == class3_param and type4 == class4_param and type5 == class5_param):
                                    # Find chiral atom index
                                    chiralIndex = AmoebaTorsionTorsionForceBuilder._getChiralAtomIndex(atoms, ib, ic, id)
                                    AmoebaTorsionTorsionForceBuilder.addTorsionTorsion(force, ia, ib, ic, id, ie, chiralIndex, gridIndex)

                                # Check reverse direction
                                elif (type5 == class1_param and type4 == class2_param and
                                      type3 == class3_param and type2 == class4_param and type1 == class5_param):
                                    # Find chiral atom index
                                    chiralIndex = AmoebaTorsionTorsionForceBuilder._getChiralAtomIndex(atoms, ib, ic, id)
                                    AmoebaTorsionTorsionForceBuilder.addTorsionTorsion(force, ie, id, ic, ib, ia, chiralIndex, gridIndex)

    @staticmethod
    def _getChiralAtomIndex(data, atomB, atomC, atomD):
        """
        Get the chiral atom index based on the TINKER algorithm.

        Parameters
        ----------
        data : ForceFieldData or list
            ForceFieldData object (for ForceField usage) or list of TinkerAtom objects (for TinkerFiles usage)
        atomB : int
            Index of atom B
        atomC : int
            Index of atom C (central atom)
        atomD : int
            Index of atom D

        Returns
        -------
        int
            Index of the chiral atom, or -1 if not found
        """
        chiralAtomIndex = -1

        # Check if we're dealing with ForceField data (has atomBonds attribute) or TinkerAtom list
        if hasattr(data, 'atomBonds'):
            # ForceField usage - data is ForceFieldData object
            atoms = data.atoms
            atomBonds = data.atomBonds
            bonds = data.bonds
            atomType = data.atomType

            # Get atoms bonded to atomC
            bondedAtoms = []
            for bondIndex in atomBonds[atomC]:
                bond = bonds[bondIndex]
                if bond.atom1 == atomC:
                    bondedAtoms.append(bond.atom2)
                else:
                    bondedAtoms.append(bond.atom1)

            # If atomC has four bonds, find the two bonds that do not include atomB and atomD
            if len(bondedAtoms) == 4:
                atomE = -1
                atomF = -1

                for bondedAtom in bondedAtoms:
                    if bondedAtom != atomB and bondedAtom != atomD:
                        if atomE == -1:
                            atomE = bondedAtom
                        else:
                            atomF = bondedAtom

                # Raise error if atoms E or F not found
                if atomE == -1 or atomF == -1:
                    raise ValueError(f"getChiralAtomIndex: error getting bonded partners of atomC={atomC}")

                # Check for different type between atoms E & F
                typeE = atomType[atoms[atomE]]
                typeF = atomType[atoms[atomF]]

                if typeE and typeF:
                    # Compare atom types as strings/identifiers
                    if str(typeE) > str(typeF):
                        chiralAtomIndex = atomE
                    elif str(typeF) > str(typeE):
                        chiralAtomIndex = atomF

                # If types are same, check masses
                if chiralAtomIndex == -1:
                    massE = atoms[atomE].element.mass.value_in_unit(amu) if atoms[atomE].element else 0.0
                    massF = atoms[atomF].element.mass.value_in_unit(amu) if atoms[atomF].element else 0.0
                    if massE > massF:
                        chiralAtomIndex = atomE
                    elif massF > massE:
                        chiralAtomIndex = atomF
        else:
            # TinkerFiles usage - data is list of TinkerAtom objects
            atoms = data

            # If atomC has four bonds, find the two bonds that do not include atomB and atomD
            # Set chiralAtomIndex to one of these, if they are not the same atom (type/mass)
            atomC_bonds = atoms[atomC].bonds
            if len(atomC_bonds) == 4:
                atomE = -1
                atomF = -1

                for bondedAtom in atomC_bonds:
                    if bondedAtom != atomB and bondedAtom != atomD:
                        if atomE == -1:
                            atomE = bondedAtom
                        else:
                            atomF = bondedAtom

                # Raise error if atoms E or F not found
                if atomE == -1 or atomF == -1:
                    raise ValueError(f"getChiralAtomIndex: error getting bonded partners of atomC={atomC}")

                # Check for different type/mass between atoms E & F
                typeE = atoms[atomE].atomClass
                typeF = atoms[atomF].atomClass

                if typeE and typeF:
                    if int(typeE) > int(typeF):
                        chiralAtomIndex = atomE
                    elif int(typeF) > int(typeE):
                        chiralAtomIndex = atomF

                # If types are same, check masses
                if chiralAtomIndex == -1:
                    massE = atoms[atomE].mass if atoms[atomE].mass else 0.0
                    massF = atoms[atomF].mass if atoms[atomF].mass else 0.0
                    if massE > massF:
                        chiralAtomIndex = atomE
                    elif massF > massE:
                        chiralAtomIndex = atomF

        return chiralAtomIndex


class AmoebaWcaDispersionForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for WCA Dispersion force for AMOEBA force field.

    Attributes
    ----------
    epso : float
        Oxygen epsilon parameter.
    epsh : float
        Hydrogen epsilon parameter.
    rmino : float
        Oxygen minimum radius parameter.
    rminh : float
        Hydrogen minimum radius parameter.
    awater : float
        Water A parameter.
    slevy : float
        Slevy parameter.
    dispoff : float
        Dispersion offset parameter.
    shctd : float
        SHCTD parameter.
    classParams : Dict[Any, Tuple[float, float]]
        Dictionary mapping atom classes to (radius, epsilon) parameters.
    """

    def __init__(self, epso: float, epsh: float, rmino: float, rminh: float, awater: float, slevy: float, dispoff: float, shctd: float) -> None:
        """
        Initialize the AmoebaWcaDispersionForceBuilder.

        Parameters
        ----------
        epso : float
            Oxygen epsilon parameter.
        epsh : float
            Hydrogen epsilon parameter.
        rmino : float
            Oxygen minimum radius parameter.
        rminh : float
            Hydrogen minimum radius parameter.
        awater : float
            Water A parameter.
        slevy : float
            Slevy parameter.
        dispoff : float
            Dispersion offset parameter.
        shctd : float
            SHCTD parameter.
        """
        super().__init__()
        self.epso = epso
        self.epsh = epsh
        self.rmino = rmino
        self.rminh = rminh
        self.awater = awater
        self.slevy = slevy
        self.dispoff = dispoff
        self.shctd = shctd
        self.classParams = {}

    def registerParams(self, atomClass: Any, radius: float, epsilon: float) -> None:
        """
        Register atom parameters for WCA dispersion.

        Parameters
        ----------
        atomClass : Any
            The atom class identifier.
        radius : float
            The van der Waals radius.
        epsilon : float
            The well depth parameter.
        """
        self.classParams[atomClass] = (radius, epsilon)

    def getForce(self, sys: mm.System) -> mm.AmoebaWcaDispersionForce:
        """
        Get or create the AmoebaWcaDispersionForce in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.

        Returns
        -------
        mm.AmoebaWcaDispersionForce
            The AmoebaWcaDispersionForce instance.
        """
        def createForce():
            force = mm.AmoebaWcaDispersionForce()
            force.setEpso(self.epso)    
            force.setEpsh(self.epsh)     
            force.setRmino(self.rmino)     
            force.setRminh(self.rminh)   
            force.setDispoff(self.dispoff)
            force.setSlevy(self.slevy)
            force.setAwater(self.awater) 
            force.setShctd(self.shctd)
            return force

        return self._createOrGetForce(sys, mm.AmoebaWcaDispersionForce, createForce)

    def addParticles(self, force: mm.AmoebaWcaDispersionForce, atomClasses: List[Any]) -> None:
        """
        Add particles to the WCA dispersion force.

        Parameters
        ----------
        force : mm.AmoebaWcaDispersionForce
            The AmoebaWcaDispersionForce instance to add particles to.
        atomClasses : List[Any]
            List of atom classes indexed by atom index.
        """
        for atomCls in atomClasses:
            values = self.classParams[atomCls]
            force.addParticle(*values)


class AmoebaGeneralizedKirkwoodForceBuilder(BaseAmoebaForceBuilder):
    """
    Builder for Generalized Kirkwood force for AMOEBA force field.

    Attributes
    ----------
    solventDielectric : float
        The dielectric constant of the solvent.
    soluteDielectric : float
        The dielectric constant of the solute.
    includeCavityTerm : bool
        Whether to include the cavity term in the calculation.
    probeRadius : float
        The probe radius for surface area calculations.
    surfaceAreaFactor : float
        The surface area factor for cavity term calculations.
    atomParams : List[float]
        List of atomic charges.
    """

    def __init__(self, solventDielectric: float, soluteDielectric: float, includeCavityTerm: bool, probeRadius: float, surfaceAreaFactor: float) -> None:
        """
        Initialize the AmoebaGeneralizedKirkwoodForceBuilder.

        Parameters
        ----------
        solventDielectric : float
            The dielectric constant of the solvent.
        soluteDielectric : float
            The dielectric constant of the solute.
        includeCavityTerm : bool
            Whether to include the cavity term in the calculation.
        probeRadius : float
            The probe radius for surface area calculations.
        surfaceAreaFactor : float
            The surface area factor for cavity term calculations.
        """
        super().__init__()
        self.solventDielectric = solventDielectric
        self.soluteDielectric = soluteDielectric
        self.includeCavityTerm = includeCavityTerm
        self.probeRadius = probeRadius
        self.surfaceAreaFactor = surfaceAreaFactor
        self.atomParams = []

    def registerParams(self, charge: float) -> None:
        """
        Register atom parameters for GK force.

        Parameters
        ----------
        charge : float
            The atomic charge.
        """
        self.atomParams.append(charge)

    def getForce(self, sys: mm.System, implicitSolvent: bool = False) -> Optional[mm.AmoebaGeneralizedKirkwoodForce]:
        """
        Get or create the AmoebaGeneralizedKirkwoodForce in the system.

        Parameters
        ----------
        sys : mm.System
            The system to get the force from or add it to.
        implicitSolvent : bool, default=False
            Whether implicit solvent is being used.

        Returns
        -------
        Optional[mm.AmoebaGeneralizedKirkwoodForce]
            The AmoebaGeneralizedKirkwoodForce instance, or None if not using implicit solvent.
        """
        if not implicitSolvent:
            return None

        def createForce():
            force = mm.AmoebaGeneralizedKirkwoodForce()
            force.setSolventDielectric(self.solventDielectric)
            force.setSoluteDielectric(self.soluteDielectric)
            force.setIncludeCavityTerm(self.includeCavityTerm)
            force.setProbeRadius(self.probeRadius)
            force.setSurfaceAreaFactor(self.surfaceAreaFactor)
            return force
        return self._createOrGetForce(sys, mm.AmoebaGeneralizedKirkwoodForce, createForce)

    def addParticles(self, force: mm.AmoebaGeneralizedKirkwoodForce, atoms: List[Any], bonds: List[Tuple[int, int]], radiusType: str = "Bondi") -> None:
        """
        Add particles to the Generalized Kirkwood force.

        Parameters
        ----------
        force : mm.AmoebaGeneralizedKirkwoodForce
            The AmoebaGeneralizedKirkwoodForce instance to add particles to.
        atoms : List[Any]
            List of atom objects.
        bonds : List[Tuple[int, int]]
            List of bond indices as tuples of (atom1, atom2).
        radiusType : str, default="Bondi"
            The type of atomic radius to use ("Bondi" or "amoeba").
        """
        # Create bonded particle sets
        bondedParticleSets = [set() for _ in range(len(atoms))]
        for atom1, atom2 in bonds:
            bondedParticleSets[atom1].add(atom2)
            bondedParticleSets[atom2].add(atom1)

        for i, atom in enumerate(atoms):
            bondedAtomicNumbers = []
            for bondedIndex in bondedParticleSets[i]:
                bondedAtom = atoms[bondedIndex]
                bondedAtomicNumbers.append(self._getAtomicNumber(bondedAtom))
            atomicNumber = self._getAtomicNumber(atom)

            if radiusType.lower() == "bondi":
                radius = self.getBondiRadius(atomicNumber)
            elif radiusType.lower() == "amoeba":
                radius = self.getAtomicRadius(atomicNumber, bondedAtomicNumbers)
            else:
                raise ValueError(f"Unknown radius type for GK: {radiusType}")

            charge = self.atomParams[i]
            shct = 0.69 # self.getOverlapScaleFactor(atomicNumber)
            force.addParticle(charge, radius, shct)

    @staticmethod
    def getOverlapScaleFactor(atomicNumber: int) -> float:
        """
        Get the overlap scale factor based on atomic number for Generalized Kirkwood.

        Parameters
        ----------
        atomicNumber : int
            The atomic number.

        Returns
        -------
        float
            The overlap scale factor.

        Raises
        ------
        ValueError
            If no overlap scale factor is available for the given atomic number.
        """
        shctMap = {
            1: 0.85,   # H
            6: 0.72,   # C
            7: 0.79,   # N
            8: 0.85,   # O
            9: 0.88,   # F
            15: 0.86,  # P
            16: 0.96,  # S
            26: 0.88,  # Fe
        }
        if atomicNumber not in shctMap:
            raise ValueError(f"No overlap scale factor available for atomic number {atomicNumber}")
        else:
            return shctMap[atomicNumber]

    @staticmethod
    def getAtomicRadius(atomicNumber: int, bondedAtomicNumbers: Optional[List[int]] = None) -> float:
        """
        Get the atomic radius based on atomic number and bonded atoms for Generalized Kirkwood.

        Parameters
        ----------
        atomicNumber : int
            The atomic number.
        bondedAtomicNumbers : Optional[List[int]], default=None
            List of atomic numbers of bonded atoms.

        Returns
        -------
        float
            The atomic radius.

        Raises
        ------
        ValueError
            If no radius is available for the given atomic number.
        """
        if bondedAtomicNumbers is None:
            bondedAtomicNumbers = []

        if atomicNumber == 1:  # H
            radius = 0.132
            if 7 in bondedAtomicNumbers:  # bonded to N
                radius = 0.11
            elif 8 in bondedAtomicNumbers:  # bonded to O
                radius = 0.105
        elif atomicNumber == 3:  # Li
            radius = 0.15
        elif atomicNumber == 6:  # C
            radius = 0.20
            if len(bondedAtomicNumbers) == 3:
                radius = 0.205
            elif len(bondedAtomicNumbers) == 4:
                for bondedAtomicNumber in bondedAtomicNumbers:
                    if bondedAtomicNumber in [7, 8]:  # N or O
                        radius = 0.175
                        break
        elif atomicNumber == 7:  # N
            radius = 0.16
        elif atomicNumber == 8:  # O
            radius = 0.155
            if len(bondedAtomicNumbers) == 2:
                radius = 0.145
        elif atomicNumber == 9:  # F
            radius = 0.154
        elif atomicNumber == 10:  # Ne
            radius = 0.146
        elif atomicNumber == 11:  # Na
            radius = 0.209
        elif atomicNumber == 12:  # Mg
            radius = 0.179
        elif atomicNumber == 14:  # Si
            radius = 0.189
        elif atomicNumber == 15:  # P
            radius = 0.196
        elif atomicNumber == 16:  # S
            radius = 0.186
        elif atomicNumber == 17:  # Cl
            radius = 0.182
        elif atomicNumber == 18:  # Ar
            radius = 0.179
        elif atomicNumber == 19:  # K
            radius = 0.223
        elif atomicNumber == 20:  # Ca
            radius = 0.191
        elif atomicNumber == 35:  # Br
            radius = 2.00
        elif atomicNumber == 36:  # Kr
            radius = 0.190
        elif atomicNumber == 37:  # Rb
            radius = 0.226
        elif atomicNumber == 53:  # I
            radius = 0.237
        elif atomicNumber == 54:  # Xe
            radius = 0.207
        elif atomicNumber == 55:  # Cs
            radius = 0.263
        elif atomicNumber == 56:  # Ba
            radius = 0.230
        else:
            raise ValueError(f"No GK radius available for atomic number {atomicNumber}")

        return radius

    @staticmethod
    def getBondiRadius(atomicNumber: int, scaleFactor: float = 1.03) -> float:
        """
        Get Bondi radius based on atomic number with optional scaling for Generalized Kirkwood.

        Parameters
        ----------
        atomicNumber : int
            The atomic number.
        scaleFactor : float, default=1.03
            Scaling factor to apply to the Bondi radius.

        Returns
        -------
        float
            The scaled Bondi radius.

        Raises
        ------
        ValueError
            If no Bondi radius is available for the given atomic number.
        """
        bondiMap = {
            0: 0.00,
            1: 0.12,    # H
            2: 0.14,    # He
            5: 0.18,    # B
            6: 0.170,   # C
            7: 0.155,   # N
            8: 0.152,   # O
            9: 0.147,   # F
            10: 0.154,  # Ne
            14: 0.210,  # Si
            15: 0.180,  # P
            16: 0.180,  # S
            17: 0.175,  # Cl
            18: 0.188,  # Ar
            34: 0.190,  # Se
            35: 0.185,  # Br
            36: 0.202,  # Kr
            53: 0.198,  # I
            54: 0.216,  # Xe
        }

        if atomicNumber not in bondiMap:
            raise ValueError(f"No Bondi radius available for atomic number {atomicNumber}")

        return bondiMap[atomicNumber] * scaleFactor


class MultipoleParams(namedtuple('MultipoleParams', ['kIndices', 'charge', 'dipole', 'quadrupole', 'axisType'])):
    def __new__(cls, kIndices, charge, dipole, quadrupole):
        # Tinker encodes the axis type based on the number of indices, and whether they are positive or negative.
        # Work out the axis type and correct atom type indices.

        while len(kIndices) < 4:
            kIndices.append(0)
        kz, kx, ky = kIndices[1:4]
        axisType = mm.AmoebaMultipoleForce.ZThenX
        if kz == 0:
            axisType = mm.AmoebaMultipoleForce.NoAxisType
        if kz != 0 and kx == 0:
            axisType = mm.AmoebaMultipoleForce.ZOnly
        if kz < 0 or kx < 0:
            axisType = mm.AmoebaMultipoleForce.Bisector
        if kx < 0 and ky < 0:
            axisType = mm.AmoebaMultipoleForce.ZBisect
        if kz < 0 and kx < 0 and ky  < 0:
            axisType = mm.AmoebaMultipoleForce.ThreeFold
        kIndices[1] = abs(kz)
        kIndices[2] = abs(kx)
        kIndices[3] = abs(ky)

        return tuple.__new__(cls, (kIndices, charge, dipole, quadrupole, axisType))

PolarizationParams = namedtuple('PolarizationParams', ['polarizability', 'thole', 'groupAtomTypes'])


class AmoebaMultipoleForceBuilder(BaseAmoebaForceBuilder):
    """
    Multipole force for AMOEBA force field.

    Attributes
    ----------
    multipoleParams : Dict[Any, List[MultipoleParams]]
        Dictionary mapping atom types to lists of multipole parameters.
    polarizationParams : Dict[Any, PolarizationParams]
        Dictionary mapping atom types to polarization parameters.
    """

    def __init__(self) -> None:
        """Initialize the AmoebaMultipoleForceBuilder."""
        super().__init__()
        self.multipoleParams = defaultdict(list)
        self.polarizationParams = defaultdict(dict)

    def registerMultipoleParams(self, type: Any, params: MultipoleParams) -> None:
        """
        Register multipole parameters for an atom type.

        Parameters
        ----------
        type : Any
            The atom type identifier.
        params : MultipoleParams
            The multipole parameters.
        """
        self.multipoleParams[type].append(params)

    def registerPolarizationParams(self, type: Any, params: PolarizationParams) -> None:
        """
        Register polarization parameters for an atom type.

        Parameters
        ----------
        type : Any
            The atom type identifier.
        params : PolarizationParams
            The polarization parameters.
        """
        self.polarizationParams[type] = params

    def getForce(self, sys: mm.System, nonbondedMethod: Any, nonbondedCutoff: float, ewaldErrorTolerance: float, polarization: str, mutualInducedTargetEpsilon: float, mutualInducedMaxIterations: int) -> mm.AmoebaMultipoleForce:
        """
        Get the AmoebaMultipoleForce. If there is not already one present in the System, create a new one and add it.

        Parameters
        ----------
        sys : mm.System
            The OpenMM System to add the force to.
        nonbondedMethod : Any
            The nonbonded method to use.
        nonbondedCutoff : float
            The cutoff distance for nonbonded interactions.
        ewaldErrorTolerance : float
            The error tolerance for Ewald summation.
        polarization : str
            The polarization type ('direct', 'extrapolated', or 'mutual').
        mutualInducedTargetEpsilon : float
            Target epsilon for mutual induced dipole iterations.
        mutualInducedMaxIterations : int
            Maximum iterations for mutual induced dipole calculations.

        Returns
        -------
        mm.AmoebaMultipoleForce
            The AmoebaMultipoleForce instance.
        """
        def createForce():
            force = mm.AmoebaMultipoleForce()
            methodMap = {ff.NoCutoff: mm.AmoebaMultipoleForce.NoCutoff,
                         ff.PME: mm.AmoebaMultipoleForce.PME}
            if nonbondedMethod not in methodMap:
                raise ValueError('Invalid nonbonded method for AmoebaMultipoleForce')
            force.setNonbondedMethod(methodMap[nonbondedMethod])
            force.setCutoffDistance(nonbondedCutoff)
            force.setEwaldErrorTolerance(ewaldErrorTolerance)
            if polarization.lower() == 'direct':
                force.setPolarizationType(mm.AmoebaMultipoleForce.Direct)
            elif polarization.lower() == 'extrapolated':
                force.setPolarizationType(mm.AmoebaMultipoleForce.Extrapolated)
            elif polarization.lower() == 'mutual':
                force.setPolarizationType(mm.AmoebaMultipoleForce.Mutual)
            else:
                raise ValueError('Invalid polarization type for AmoebaMultipoleForce: '+polarization)
            force.setMutualInducedTargetEpsilon(mutualInducedTargetEpsilon)
            force.setMutualInducedMaxIterations(mutualInducedMaxIterations)
            return force

        return self._createOrGetForce(sys, mm.AmoebaMultipoleForce, createForce)

    def addMultipoles(self, force: mm.AmoebaMultipoleForce, atomTypes: List[Any], atoms: List[Any], bonds: List[Tuple[int, int]]) -> None:
        """
        Add multipoles to the AmoebaMultipoleForce.

        Parameters
        ----------
        force : mm.AmoebaMultipoleForce
            The AmoebaMultipoleForce instance to add multipoles to.
        atomTypes : List[Any]
            List of atom types indexed by atom index.
        atoms : List[Any]
            List of atom objects.
        bonds : List[Tuple[int, int]]
            List of bond indices as tuples of (atom1, atom2).
        """
        self.buildBondedParticleSets(len(atoms), bonds)
        for atomIndex, t in enumerate(atomTypes):
            if t in self.multipoleParams:
                multipoleList = self.multipoleParams[t]
                hit = False
                savedMultipoleParams = None
                for multipoleParams in multipoleList:
                    if hit:
                        break
                    kz, kx, ky = multipoleParams.kIndices[1:4]
                    bondedAtomIndices = self.bonded12ParticleSets[atomIndex]
                    zaxis = -1
                    xaxis = -1
                    yaxis = -1
                    for bondedAtomZIndex in bondedAtomIndices:
                       if hit:
                           break
                       bondedAtomZType = atomTypes[bondedAtomZIndex]
                       if bondedAtomZType == kz:
                          for bondedAtomXIndex in bondedAtomIndices:
                              if bondedAtomXIndex == bondedAtomZIndex or hit:
                                  continue
                              bondedAtomXType = atomTypes[bondedAtomXIndex]
                              if bondedAtomXType == kx:
                                  if ky == 0:
                                      zaxis = bondedAtomZIndex
                                      xaxis = bondedAtomXIndex
                                      if bondedAtomXType == bondedAtomZType and xaxis < zaxis:
                                          swapI = zaxis
                                          zaxis = xaxis
                                          xaxis = swapI
                                      else:
                                          for bondedAtomXIndex2 in bondedAtomIndices:
                                              bondedAtomX1Type = atomTypes[bondedAtomXIndex2]
                                              if bondedAtomX1Type == kx and bondedAtomXIndex2 != bondedAtomZIndex and bondedAtomXIndex2 < xaxis:
                                                  xaxis = bondedAtomXIndex2
                                      savedMultipoleParams = multipoleParams
                                      hit = True
                                  else:
                                      for bondedAtomYIndex in bondedAtomIndices:
                                          if bondedAtomYIndex == bondedAtomZIndex or bondedAtomYIndex == bondedAtomXIndex or hit:
                                              continue
                                          bondedAtomYType = atomTypes[bondedAtomYIndex]
                                          if bondedAtomYType == ky:
                                              zaxis = bondedAtomZIndex
                                              xaxis = bondedAtomXIndex
                                              yaxis = bondedAtomYIndex
                                              savedMultipoleParams = multipoleParams
                                              hit = True

                # assign multipole parameters via 1-2 and 1-3 connected atoms

                for multipoleParams in multipoleList:
                    if hit:
                        break
                    kz, kx, ky = multipoleParams.kIndices[1:4]

                    # assign multipole parameters
                    #    (1) get bonded partners
                    #    (2) match parameter types

                    bondedAtom12Indices = self.bonded12ParticleSets[atomIndex]
                    bondedAtom13Indices = self.bonded13ParticleSets[atomIndex]
                    zaxis = -1
                    xaxis = -1
                    yaxis = -1
                    for bondedAtomZIndex in bondedAtom12Indices:
                       if hit:
                           break
                       bondedAtomZType = atomTypes[bondedAtomZIndex]
                       if bondedAtomZType == kz:
                          for bondedAtomXIndex in bondedAtom13Indices:
                              if bondedAtomXIndex == bondedAtomZIndex or hit:
                                  continue
                              bondedAtomXType = atomTypes[bondedAtomXIndex]
                              if bondedAtomXType == kx and bondedAtomZIndex in self.bonded12ParticleSets[bondedAtomXIndex]:
                                  if ky == 0:
                                      zaxis = bondedAtomZIndex
                                      xaxis = bondedAtomXIndex

                                      # select xaxis w/ smallest index

                                      for bondedAtomXIndex2 in bondedAtom13Indices:
                                          bondedAtomX1Type = atomTypes[bondedAtomXIndex2]
                                          if bondedAtomX1Type == kx and bondedAtomXIndex2 != bondedAtomZIndex and bondedAtomZIndex in self.bonded12ParticleSets[bondedAtomXIndex2] and bondedAtomXIndex2 < xaxis:
                                              xaxis = bondedAtomXIndex2

                                      savedMultipoleParams = multipoleParams
                                      hit = True
                                  else:
                                      for bondedAtomYIndex in bondedAtom13Indices:
                                          if bondedAtomYIndex == bondedAtomZIndex or bondedAtomYIndex == bondedAtomXIndex or hit:
                                              continue
                                          bondedAtomYType = atomTypes[bondedAtomYIndex]
                                          if bondedAtomYType == ky and bondedAtomZIndex in self.bonded12ParticleSets[bondedAtomYIndex]:
                                              zaxis = bondedAtomZIndex
                                              xaxis = bondedAtomXIndex
                                              yaxis = bondedAtomYIndex
                                              savedMultipoleParams = multipoleParams
                                              hit = True

                # assign multipole parameters via only a z-defining atom

                for multipoleParams in multipoleList:
                    if hit:
                        break
                    kz, kx = multipoleParams.kIndices[1:3]
                    zaxis = -1
                    xaxis = -1
                    yaxis = -1
                    for bondedAtomZIndex in bondedAtom12Indices:
                        if hit:
                            break
                        bondedAtomZType = atomTypes[bondedAtomZIndex]
                        if kx == 0 and kz == bondedAtomZType:
                            zaxis = bondedAtomZIndex
                            savedMultipoleParams = multipoleParams
                            hit = True

                # assign multipole parameters via no connected atoms

                for multipoleParams in multipoleList:
                    if hit:
                        break
                    kz = multipoleParams.kIndices[1]
                    zaxis = -1
                    xaxis = -1
                    yaxis = -1
                    if kz == 0:
                        savedMultipoleParams = multipoleParams
                        hit = True                # add particle if there was a hit

                if hit:
                    polarizationParams = self.polarizationParams[t]
                    pdamp = 0 if polarizationParams.thole == 0 else pow(polarizationParams.polarizability, 1.0/6.0)
                    newIndex = force.addMultipole(savedMultipoleParams.charge, savedMultipoleParams.dipole, savedMultipoleParams.quadrupole, savedMultipoleParams.axisType,
                                                                 zaxis, xaxis, yaxis, polarizationParams.thole, pdamp, polarizationParams.polarizability)
                    if atomIndex == newIndex:
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent12, tuple(self.bonded12ParticleSets[atomIndex]))
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent13, tuple(self.bonded13ParticleSets[atomIndex]))
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent14, tuple(self.bonded14ParticleSets[atomIndex]))
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent15, tuple(self.bonded15ParticleSets[atomIndex]))
                    else:
                        atom = atoms[atomIndex]
                        raise ValueError("Atom %d (%s of %s %d) is out of sync!." %(atomIndex, atom.name, atom.residue.name, atom.residue.index))
                else:
                    atom = atoms[atomIndex]
                    raise ValueError('Atom %d (%s %s %d) was not assigned multipole parameters' % (atomIndex, atom.name, atom.residue.name, atom.residue.index))
            else:
                atom = atoms[atomIndex]
                raise ValueError('No multipole type for atom %d (%s %s %d)' % (atomIndex, atom.name, atom.residue.name, atom.residue.index))

        # Set up the polarization groups.

        self.setPolarGroups(force, atomTypes)

    def buildBondedParticleSets(self, numAtoms, bonds):
        """Identify sets of particles that are separated by various numbers of bonds."""

        # 1-2

        bonded12ParticleSets = [set() for _ in range(numAtoms)]
        for atom1, atom2 in bonds:
            bonded12ParticleSets[atom1].add(atom2)
            bonded12ParticleSets[atom2].add(atom1)

        # 1-3

        bonded13ParticleSets = []
        for i in range(numAtoms):
            bonded13Set = set()
            bonded12ParticleSet = bonded12ParticleSets[i]
            for j in bonded12ParticleSet:
                bonded13Set = bonded13Set.union(bonded12ParticleSets[j])

            # remove 1-2 and self from set

            bonded13Set = bonded13Set - bonded12ParticleSet
            selfSet = set()
            selfSet.add(i)
            bonded13Set = bonded13Set - selfSet
            bonded13Set = set(sorted(bonded13Set))
            bonded13ParticleSets.append(bonded13Set)

        # 1-4

        bonded14ParticleSets = []
        for i in range(numAtoms):
            bonded14Set = set()
            bonded13ParticleSet = bonded13ParticleSets[i]
            for j in bonded13ParticleSet:
                bonded14Set = bonded14Set.union(bonded12ParticleSets[j])

            # remove 1-3, 1-2 and self from set

            bonded14Set = bonded14Set - bonded12ParticleSets[i]
            bonded14Set = bonded14Set - bonded13ParticleSet
            selfSet = set()
            selfSet.add(i)
            bonded14Set = bonded14Set - selfSet
            bonded14Set = set(sorted(bonded14Set))
            bonded14ParticleSets.append(bonded14Set)

        # 1-5

        bonded15ParticleSets = []
        for i in range(numAtoms):
            bonded15Set = set()
            bonded14ParticleSet = bonded14ParticleSets[i]
            for j in bonded14ParticleSet:
                bonded15Set = bonded15Set.union(bonded12ParticleSets[j])

            # remove 1-4, 1-3, 1-2 and self from set

            bonded15Set = bonded15Set - bonded12ParticleSets[i]
            bonded15Set = bonded15Set - bonded13ParticleSets[i]
            bonded15Set = bonded15Set - bonded14ParticleSet
            selfSet = set()
            selfSet.add(i)
            bonded15Set = bonded15Set - selfSet
            bonded15Set = set(sorted(bonded15Set))
            bonded15ParticleSets.append(bonded15Set)

        self.bonded12ParticleSets = bonded12ParticleSets
        self.bonded13ParticleSets = bonded13ParticleSets
        self.bonded14ParticleSets = bonded14ParticleSets
        self.bonded15ParticleSets = bonded15ParticleSets

    def setPolarGroups(self, force, atomTypes):
        """Assign covalent maps to atoms based on polarization groups."""
        groupSetsForAtom = [[] for _ in range(force.getNumMultipoles())]
        polarizationGroupForAtom = [set() for _ in range(force.getNumMultipoles())]
        for atomIndex in range(len(atomTypes)):

            # assign multipole parameters via only 1-2 connected atoms

            groupTypes = self.polarizationParams[atomTypes[atomIndex]].groupAtomTypes
            bondedAtomIndices = self.bonded12ParticleSets[atomIndex]
            polarizationGroupForAtom[atomIndex].add(atomIndex)
            for bondedAtomIndex in bondedAtomIndices:
                bondedAtomType = atomTypes[bondedAtomIndex]
                if bondedAtomType in groupTypes:
                    polarizationGroupForAtom[atomIndex].add(bondedAtomIndex)
                    polarizationGroupForAtom[bondedAtomIndex].add(atomIndex)

        # pgrp11

        for atomIndex in range(len(atomTypes)):
            if len(groupSetsForAtom[atomIndex]) > 0:
                continue
            group = set()
            visited = set()
            notVisited = set()
            for pgrpAtomIndex in polarizationGroupForAtom[atomIndex]:
                group.add(pgrpAtomIndex)
                notVisited.add(pgrpAtomIndex)
            visited.add(atomIndex)
            while len(notVisited) > 0:
                nextAtom = notVisited.pop()
                if nextAtom not in visited:
                   visited.add(nextAtom)
                   for ii in polarizationGroupForAtom[nextAtom]:
                       group.add(ii)
                       if ii not in visited:
                           notVisited.add(ii)
            pGroup = group
            for pgrpAtomIndex in group:
                groupSetsForAtom[pgrpAtomIndex].append(pGroup)

        for atomIndex in range(len(atomTypes)):
            groupSetsForAtom[atomIndex][0] = sorted(groupSetsForAtom[atomIndex][0])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent11, groupSetsForAtom[atomIndex][0])

        # pgrp12

        for atomIndex in range(len(atomTypes)):
            if len(groupSetsForAtom[atomIndex]) > 1:
                continue
            pgrp11 = set(groupSetsForAtom[atomIndex][0])
            pgrp12 = set()
            for pgrpAtomIndex in pgrp11:
                for bonded12 in self.bonded12ParticleSets[pgrpAtomIndex]:
                    pgrp12 = pgrp12.union(groupSetsForAtom[bonded12][0])
            pgrp12 = pgrp12 - pgrp11
            for pgrpAtomIndex in pgrp11:
                groupSetsForAtom[pgrpAtomIndex].append(pgrp12)

        for atomIndex in range(len(atomTypes)):
            groupSetsForAtom[atomIndex][1] = sorted(groupSetsForAtom[atomIndex][1])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent12, groupSetsForAtom[atomIndex][1])

        # pgrp13

        for atomIndex in range(len(atomTypes)):
            if len(groupSetsForAtom[atomIndex]) > 2:
                continue
            pgrp11 = set(groupSetsForAtom[atomIndex][0])
            pgrp12 = set(groupSetsForAtom[atomIndex][1])
            pgrp13 = set()
            for pgrpAtomIndex in pgrp12:
                for bonded12 in self.bonded12ParticleSets[pgrpAtomIndex]:
                    pgrp13 = pgrp13.union(groupSetsForAtom[bonded12][0])
            pgrp13 = pgrp13 - pgrp12
            pgrp13 = pgrp13 - set(pgrp11)
            for pgrpAtomIndex in pgrp11:
                groupSetsForAtom[pgrpAtomIndex].append(pgrp13)

        for atomIndex in range(len(atomTypes)):
            groupSetsForAtom[atomIndex][2] = sorted(groupSetsForAtom[atomIndex][2])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent13, groupSetsForAtom[atomIndex][2])

        # pgrp14

        for atomIndex in range(len(atomTypes)):
            if len(groupSetsForAtom[atomIndex]) > 3:
                continue
            pgrp11 = set(groupSetsForAtom[atomIndex][0])
            pgrp12 = set(groupSetsForAtom[atomIndex][1])
            pgrp13 = set(groupSetsForAtom[atomIndex][2])
            pgrp14 = set()
            for pgrpAtomIndex in pgrp13:
                for bonded12 in self.bonded12ParticleSets[pgrpAtomIndex]:
                    pgrp14 = pgrp14.union(groupSetsForAtom[bonded12][0])
            pgrp14 = pgrp14 - pgrp13
            pgrp14 = pgrp14 - pgrp12
            pgrp14 = pgrp14 - set(pgrp11)
            for pgrpAtomIndex in pgrp11:
                groupSetsForAtom[pgrpAtomIndex].append(pgrp14)

        for atomIndex in range(len(atomTypes)):
            groupSetsForAtom[atomIndex][3] = sorted(groupSetsForAtom[atomIndex][3])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent14, groupSetsForAtom[atomIndex][3])


class AmoebaVdwForceBuilder(BaseAmoebaForceBuilder):
    """
    AmoebaVdwForce for AMOEBA force field.

    Attributes
    ----------
    type : str
        The VdW potential type.
    radiusrule : str
        The sigma combining rule.
    radiustype : str
        The radius type.
    radiussize : str
        The radius size definition.
    epsilonrule : str
        The epsilon combining rule.
    vdw13Scale : float
        The 1-3 scaling factor.
    classParams : Dict[Any, Tuple[float, float, float]]
        Dictionary mapping atom classes to (sigma, epsilon, reduction) parameters.
    pairParams : List[Tuple[Any, Any, float, float]]
        List of pair-specific parameters.
    """

    def __init__(self, type: str, radiusrule: str, radiustype: str, radiussize: str, epsilonrule: str, vdw13Scale: float, vdw14Scale: float, vdw15Scale: float) -> None:
        """
        Initialize the AmoebaVdwForceBuilder.

        Parameters
        ----------
        type : str
            The VdW potential type.
        radiusrule : str
            The sigma combining rule.
        radiustype : str
            The radius type.
        radiussize : str
            The radius size definition.
        epsilonrule : str
            The epsilon combining rule.
        vdw13Scale : float
            The 1-3 scaling factor.
        vdw14Scale : float
            The 1-4 scaling factor.
        vdw15Scale : float
            The 1-5 scaling factor.

        Raises
        ------
        ValueError
            If unsupported scaling factors are provided.
        """
        super().__init__()
        if vdw13Scale != 0.0 and vdw13Scale != 1.0:
            raise ValueError('AmoebaVdwForce: the only supported values for vdw-13-scale are 0 or 1')
        if vdw14Scale != 1.0:
            raise ValueError('AmoebaVdwForce: the only supported value for vdw-14-scale is 1')
        if vdw15Scale != 1.0:
            raise ValueError('AmoebaVdwForce: the only supported value for vdw-15-scale is 1')
        self.type = type.upper()
        self.radiusrule = radiusrule.upper()
        self.radiustype = radiustype
        self.radiussize = radiussize
        self.epsilonrule = epsilonrule.upper()
        self.vdw13Scale = vdw13Scale
        self.classParams = {}
        self.pairParams = []

    def registerClassParams(self, atomClass: Any, sigma: float, epsilon: float, reduction: float) -> None:
        """
        Register class-based VdW parameters.

        Parameters
        ----------
        atomClass : Any
            The atom class identifier.
        sigma : float
            The sigma parameter.
        epsilon : float
            The epsilon parameter.
        reduction : float
            The reduction factor.
        """
        self.classParams[atomClass] = (sigma, epsilon, reduction)

    def registerPairParams(self, class1: Any, class2: Any, sigma: float, epsilon: float) -> None:
        """
        Register pair-specific VdW parameters.

        Parameters
        ----------
        class1 : Any
            The first atom class identifier.
        class2 : Any
            The second atom class identifier.
        sigma : float
            The sigma parameter for this pair.
        epsilon : float
            The epsilon parameter for this pair.
        """
        self.pairParams.append((class1, class2, sigma, epsilon))

    def getForce(self, sys: mm.System, nonbondedMethod: Any, nonbondedCutoff: float, useDispersionCorrection: bool) -> mm.AmoebaVdwForce:
        """
        Get the AmoebaVdwForce. If there is not already one present in the System, create a new one and add it.

        Parameters
        ----------
        sys : mm.System
            The OpenMM System to add the force to.
        nonbondedMethod : Any
            The nonbonded method to use.
        nonbondedCutoff : float
            The cutoff distance for nonbonded interactions.
        useDispersionCorrection : bool
            Whether to use dispersion correction.

        Returns
        -------
        mm.AmoebaVdwForce
            The AmoebaVdwForce instance.
        """
        def createForce():
            force = mm.AmoebaVdwForce()
            force.setCutoff(nonbondedCutoff)
            force.setUseDispersionCorrection(useDispersionCorrection)
            if (nonbondedMethod == ff.PME):
                force.setNonbondedMethod(mm.AmoebaVdwForce.CutoffPeriodic)
            potentialMap = {'BUFFERED-14-7':0, 'LENNARD-JONES':1}
            allowedSigma = ['ARITHMETIC', 'GEOMETRIC', 'CUBIC-MEAN']
            allowedEpsilon = ['ARITHMETIC', 'GEOMETRIC', 'HARMONIC', 'W-H', 'HHG']
            if self.type not in potentialMap:
                raise ValueError("AmoebaVdwForce: potential type %s not recognized; valid values are %s." % (self.type, ', '.join(potentialMap.keys())))
            if self.radiusrule not in allowedSigma:
                raise ValueError("AmoebaVdwForce: sigma combining rule %s not recognized; valid values are %s." % (self.radiusrule, ', '.join(allowedSigma)))
            if self.epsilonrule not in allowedEpsilon:
                raise ValueError("AmoebaVdwForce: epsilon combining rule %s not recognized; valid values are %s." % (self.epsilonrule, ', '.join(allowedEpsilon)))
            force.setPotentialFunction(potentialMap[self.type])
            force.setSigmaCombiningRule(self.radiusrule)
            force.setEpsilonCombiningRule(self.epsilonrule)
            return force

        return self._createOrGetForce(sys, mm.AmoebaVdwForce, createForce)

    def addParticles(self, force: mm.AmoebaVdwForce, atomClasses: List[Any], atoms: List[Any], bonds: List[Tuple[int, int]]) -> None:
        """
        Add particles to the AmoebaVdwForce.

        Parameters
        ----------
        force : mm.AmoebaVdwForce
            The AmoebaVdwForce instance to add particles to.
        atomClasses : List[Any]
            List of atom classes indexed by atom index.
        atoms : List[Any]
            List of atom objects.
        bonds : List[Tuple[int, int]]
            List of bond indices as tuples of (atom1, atom2).
        """

        # Define types

        sigmaScale = 1
        if self.radiustype == 'SIGMA':
            sigmaScale = 1.122462048309372
        if self.radiussize == 'DIAMETER':
            sigmaScale = 0.5
        classToTypeMap = {}
        for className in self.classParams:
            sigma, epsilon, _ = self.classParams[className]
            classToTypeMap[className] = force.addParticleType(sigma*sigmaScale, epsilon)

        # Record what other atoms each atom is bonded to.

        bondedParticleSets = [set() for _ in range(len(atoms))]
        for atom1, atom2 in bonds:
            bondedParticleSets[atom1].add(atom2)
            bondedParticleSets[atom2].add(atom1)

        # add particles to force

        for i, atom in enumerate(atoms):
            className = atomClasses[i]
            if className not in self.classParams:
                raise ValueError(f"AmoebaVdwForce: No VdW parameters found for atom class {className}. "
                               f"Available classes: {list(self.classParams.keys())}")
            _, _, reduction = self.classParams[className]
            # ivIndex = index of bonded partner for hydrogens; otherwise ivIndex = particle index
            ivIndex = i
            if atom.element == element.hydrogen and len(bondedParticleSets[i]) == 1:
                ivIndex = list(bondedParticleSets[i])[0]
            if className not in classToTypeMap:
                raise ValueError(f"AmoebaVdwForce: No type mapping found for atom class {className}")
            force.addParticle(ivIndex, classToTypeMap[className], reduction)

        # Add pairs

        for c1, c2, sigma, epsilon in self.pairParams:
            force.addTypePair(classToTypeMap[c1], classToTypeMap[c2], sigma, epsilon)

        # set combining rules

        # set particle exclusions: self, 1-2 and 1-3 bonds
        # (1) collect in bondedParticleSets[i], 1-2 indices for all bonded partners of particle i
        # (2) add 1-2,1-3 and self to exclusion set

        for i in range(len(atoms)):
            # 1-2 partners

            exclusionSet = bondedParticleSets[i].copy()

            # 1-3 partners

            if self.vdw13Scale == 0.0:
                for bondedParticle in bondedParticleSets[i]:
                    exclusionSet = exclusionSet.union(bondedParticleSets[bondedParticle])

            # self

            exclusionSet.add(i)
            force.setParticleExclusions(i, tuple(exclusionSet))
