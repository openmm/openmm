"""
amoebaforces.py: AMOEBA force field classes.

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

import math
from collections import defaultdict, namedtuple
import openmm as mm
import openmm.app.forcefield as ff


class AmoebaBondForce:
    """Bond force for AMOEBA force field"""

    def __init__(self, cubic, quartic):
        self.cubic = cubic
        self.quartic = quartic

    def getForce(self, sys):
        energy = f"k*(d^2 + {self.cubic}*d^3 + {self.quartic}*d^4); d=r-r0"
        existing = [
            f
            for f in sys.getForces()
            if isinstance(f, mm.CustomBondForce) and f.getEnergyFunction() == energy
        ]
        if len(existing) == 0:
            force = mm.CustomBondForce(energy)
            force.addPerBondParameter("r0")
            force.addPerBondParameter("k")
            force.setName("AmoebaBond")
        else:
            force = existing[0]

        return force

    @staticmethod
    def addBond(force, atom1, atom2, length, k):
        force.addBond(atom1, atom2, [length, k])


class AmoebaTorsionForce:
    """PeriodicTorsionForce force for AMOEBA force field"""

    def __init__(self, torsionUnit):
        self.torsionUnit = torsionUnit

    def getForce(self, sys):
        existing = [
            f for f in sys.getForces() if isinstance(f, mm.PeriodicTorsionForce)
        ]
        if len(existing) == 0:
            force = mm.PeriodicTorsionForce()
        else:
            force = existing[0]
        return force

    @staticmethod
    def addTorsion(force, atom1, atom2, atom3, atom4, t1, t2, t3):
        if t1[0] != 0:
            force.addTorsion(atom1, atom2, atom3, atom4, 1, t1[1], t1[0])
        if t2[0] != 0:
            force.addTorsion(atom1, atom2, atom3, atom4, 2, t2[1], t2[0])
        if t3[0] != 0:
            force.addTorsion(atom1, atom2, atom3, atom4, 3, t3[1], t3[0])


class AmoebaPiTorsionForce(object):
    """PiTorsionForce force for AMOEBA force field"""

    def __init__(self, piTorsionUnit):
        self.piTorsionUnit = piTorsionUnit

    def getForce(self, sys):
        energy = """2*k*sin(phi)^2;
                    phi = pointdihedral(x3+c1x, y3+c1y, z3+c1z, x3, y3, z3, x4, y4, z4, x4+c2x, y4+c2y, z4+c2z);
                    c1x = (d14y*d24z-d14z*d24y); c1y = (d14z*d24x-d14x*d24z); c1z = (d14x*d24y-d14y*d24x);
                    c2x = (d53y*d63z-d53z*d63y); c2y = (d53z*d63x-d53x*d63z); c2z = (d53x*d63y-d53y*d63x);
                    d14x = x1-x4; d14y = y1-y4; d14z = z1-z4;
                    d24x = x2-x4; d24y = y2-y4; d24z = z2-z4;
                    d53x = x5-x3; d53y = y5-y3; d53z = z5-z3;
                    d63x = x6-x3; d63y = y6-y3; d63z = z6-z3"""
        existing = [
            f
            for f in sys.getForces()
            if isinstance(f, mm.CustomCompoundBondForce)
            and f.getEnergyFunction() == energy
        ]

        if len(existing) == 0:
            force = mm.CustomCompoundBondForce(6, energy)
            force.addPerBondParameter("k")
            force.setName("AmoebaPiTorsion")
        else:
            force = existing[0]

        return force

    @staticmethod
    def addPiTorsion(force, atom1, atom2, atom3, atom4, atom5, atom6, k):
        force.addBond([atom1, atom2, atom3, atom4, atom5, atom6], [k])


class AmoebaUreyBradleyForce:
    """UreyBradleyForce force for AMOEBA force field"""

    def __init__(self):
        pass

    def getForce(self, sys):
        existing = [f for f in sys.getForces() if isinstance(f, mm.HarmonicBondForce)]
        if len(existing) == 0:
            force = mm.HarmonicBondForce()
        else:
            force = existing[0]
        return force

    def addUreyBradley(self, force, atom1, atom3, k, d):
        force.addBond(atom1, atom3, d, 2 * k)


class AmoebaAngleForce:
    """AngleForce force for AMOEBA force field"""

    def __init__(self, cubic, quartic, pentic, sextic):
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic

    def getForce(self, sys):
        energy = "k*(d^2 + %s*d^3 + %s*d^4 + %s*d^5 + %s*d^6); d=%.15g*theta-theta0" % (
            self.cubic,
            self.quartic,
            self.pentic,
            self.sextic,
            180 / math.pi,
        )
        existing = [
            f
            for f in sys.getForces()
            if isinstance(f, mm.CustomAngleForce) and f.getEnergyFunction() == energy
        ]

        if len(existing) == 0:
            force = mm.CustomAngleForce(energy)
            force.addPerAngleParameter("theta0")
            force.addPerAngleParameter("k")
            force.setName("AmoebaAngle")
        else:
            force = existing[0]

        return force

    @staticmethod
    def addAngle(force, atom1, atom2, atom3, theta0, k):
        force.addAngle(atom1, atom2, atom3, [theta0, k])


class AmoebaInPlaneAngleForce:
    """InPlaneAngleForce force for AMOEBA force field"""

    def __init__(self, cubic, quartic, pentic, sextic):
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic

    def getForce(self, sys):
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
            180 / math.pi,
        )
        existing = [
            f
            for f in sys.getForces()
            if isinstance(f, mm.CustomCompoundBondForce)
            and f.getEnergyFunction() == energy
        ]
        if len(existing) == 0:
            force = mm.CustomCompoundBondForce(4, energy)
            force.addPerBondParameter("theta0")
            force.addPerBondParameter("k")
            force.setName("AmoebaInPlaneAngle")
        else:
            force = existing[0]

        return force

    @staticmethod
    def addInPlaneAngle(force, atom1, atom2, atom3, theta0, k):
        force.addBond(atom1, atom2, atom3, [theta0, k])


class AmoebaOutOfPlaneBendForce:
    """OutOfPlaneBendForce force for AMOEBA force field"""

    def __init__(self, cubic, quartic, pentic, sextic):
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic

    def getForce(self, sys):
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
            180 / math.pi,
        )
        existing = [
            f
            for f in sys.getForces()
            if isinstance(f, mm.CustomCompoundBondForce)
            and f.getEnergyFunction() == energy
        ]
        if len(existing) == 0:
            force = mm.CustomCompoundBondForce(4, energy)
            force.addPerBondParameter("k")
            force.setName("AmoebaOutOfPlaneBend")
        else:
            force = existing[0]

        return force

    @staticmethod
    def addOutOfPlaneBend(force, atom1, atom2, atom3, atom4, k):
        force.addBond([atom1, atom2, atom3, atom4], [k])


class AmoebaStretchBendForce:
    """StretchBendForce force for AMOEBA force field"""

    def __init__(self):
        pass

    def getForce(self, sys):
        energy = (
            "(k1*(distance(p1,p2)-r12) + k2*(distance(p2,p3)-r23))*(%.15g*(angle(p1,p2,p3)-theta0))"
            % (180 / math.pi)
        )
        existing = [
            f
            for f in sys.getForces()
            if isinstance(f, mm.CustomCompoundBondForce)
            and f.getEnergyFunction() == energy
        ]
        if len(existing) == 0:
            force = mm.CustomCompoundBondForce(3, energy)
            force.addPerBondParameter("r12")
            force.addPerBondParameter("r23")
            force.addPerBondParameter("theta0")
            force.addPerBondParameter("k1")
            force.addPerBondParameter("k2")
            force.setName("AmoebaStretchBend")
            sys.addForce(force)
        else:
            force = existing[0]

        return force

    @staticmethod
    def addStretchBend(force, atom1, atom2, atom3, r12, r23, theta0, k1, k2):
        force.addBond((atom1, atom2, atom3), (r12, r23, theta0, k1, k2))


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


class AmoebaMultipoleForceBuilder(object):
    """Multipole force for AMOEBA force field"""

    def __init__(self):
        self.multipoleParams = defaultdict(list)
        self.polarizationParams = defaultdict(dict)

    def registerMultipoleParams(self, type, params):
        self.multipoleParams[type].append(params)

    def registerPolarizationParams(self, type, params):
        self.polarizationParams[type] = params

    def getForce(self, sys, nonbondedMethod, nonbondedCutoff, ewaldErrorTolerance, polarization, mutualInducedTargetEpsilon, mutualInducedMaxIterations):
        """Get the AmoebaMultipoleForce.  If there is not already one present in the System, create a new one and add it."""
        existing = [f for f in sys.getForces() if isinstance(f, mm.AmoebaMultipoleForce)]
        if len(existing) == 0:
            force = mm.AmoebaMultipoleForce()
            sys.addForce(force)
        else:
            force = existing[0]
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

    def addMultipoles(self, force, atomTypes, atoms, bonds):
        """Add multipoles to the AmoebaMultipoleForce."""
        self.buildBondedParticleSets(len(atoms), bonds)
        for atomIndex, t in enumerate(atomTypes):
            if t in self.multipoleParams:
                multipoleList = self.multipoleParams[t]
                hit = False
                savedMultipoleParams = None

                # assign multipole parameters via only 1-2 connected atoms

                for multipoleParams in multipoleList:
                    if hit:
                        break
                    kz, kx, ky = multipoleParams.kIndices[1:4]

                    # assign multipole parameters
                    #    (1) get bonded partners
                    #    (2) match parameter types

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
                        hit = True

                # add particle if there was a hit

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
                    raise ValueError("Atom %d (%s of %s %d) was not assigned." %(atomIndex, atom.name, atom.residue.name, atom.residue.index))
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
