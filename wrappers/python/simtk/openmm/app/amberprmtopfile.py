"""
armberprmtopfile.py: Used for loading AMBER prmtop files.

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

from math import sqrt
from simtk.openmm.app import Topology
from simtk.openmm.app import PDBFile
from simtk.openmm.app.internal import amber_file_parser
from simtk.openmm.app.internal.singleton import Singleton
from . import forcefield as ff
from . import element as elem
import simtk.unit as u
import simtk.openmm as mm
from simtk.openmm.app.internal.unitcell import computePeriodicBoxVectors

# Enumerated values for implicit solvent model

class HCT(Singleton):
    def __repr__(self):
        return 'HCT'
HCT = HCT()

class OBC1(Singleton):
    def __repr__(self):
        return 'OBC1'
OBC1 = OBC1()

class OBC2(Singleton):
    def __repr__(self):
        return 'OBC2'
OBC2 = OBC2()

class GBn(Singleton):
    def __repr__(self):
        return 'GBn'
GBn = GBn()

class GBn2(Singleton):
    def __repr__(self):
        return 'GBn2'
GBn2 = GBn2()

def _strip_optunit(thing, unit):
    """
    Strips optional units, converting to specified unit type. If no unit
    present, it just returns the number
    """
    if u.is_quantity(thing):
        return thing.value_in_unit(unit)
    return thing

class AmberPrmtopFile(object):
    """AmberPrmtopFile parses an AMBER prmtop file and constructs a Topology and (optionally) an OpenMM System from it."""

    def __init__(self, file):
        """Load a prmtop file."""
        ## The Topology read from the prmtop file
        self.topology = top = Topology()
        self.elements = []

        # Load the prmtop file

        prmtop = amber_file_parser.PrmtopLoader(file)
        self._prmtop = prmtop

        # Add atoms to the topology

        PDBFile._loadNameReplacementTables()
        lastResidue = None
        c = top.addChain()
        for index in range(prmtop.getNumAtoms()):
            resNumber = prmtop.getResidueNumber(index)
            if resNumber != lastResidue:
                lastResidue = resNumber
                resName = prmtop.getResidueLabel(iAtom=index).strip()
                if resName in PDBFile._residueNameReplacements:
                    resName = PDBFile._residueNameReplacements[resName]
                r = top.addResidue(resName, c)
                if resName in PDBFile._atomNameReplacements:
                    atomReplacements = PDBFile._atomNameReplacements[resName]
                else:
                    atomReplacements = {}
            atomName = prmtop.getAtomName(index).strip()
            if atomName in atomReplacements:
                atomName = atomReplacements[atomName]

            # Get the element from the prmtop file if available
            if prmtop.has_atomic_number:
                try:
                    element = elem.Element.getByAtomicNumber(int(prmtop._raw_data['ATOMIC_NUMBER'][index]))
                except KeyError:
                    element = None
            else:
                # Try to guess the element from the atom name.

                upper = atomName.upper()
                if upper.startswith('CL'):
                    element = elem.chlorine
                elif upper.startswith('NA'):
                    element = elem.sodium
                elif upper.startswith('MG'):
                    element = elem.magnesium
                elif upper.startswith('ZN'):
                    element = elem.zinc
                else:
                    try:
                        element = elem.get_by_symbol(atomName[0])
                    except KeyError:
                        element = None

            top.addAtom(atomName, element, r)
            self.elements.append(element)

        # Add bonds to the topology

        atoms = list(top.atoms())
        for bond in prmtop.getBondsWithH():
            top.addBond(atoms[bond[0]], atoms[bond[1]])
        for bond in prmtop.getBondsNoH():
            top.addBond(atoms[bond[0]], atoms[bond[1]])

        # Set the periodic box size.

        if prmtop.getIfBox():
            box = prmtop.getBoxBetaAndDimensions()
            top.setPeriodicBoxVectors(computePeriodicBoxVectors(*(box[1:4] + box[0:1]*3)))

    def createSystem(self, nonbondedMethod=ff.NoCutoff, nonbondedCutoff=1.0*u.nanometer,
                     constraints=None, rigidWater=True, implicitSolvent=None,
                     implicitSolventSaltConc=0.0*(u.moles/u.liter),
                     implicitSolventKappa=None, temperature=298.15*u.kelvin,
                     soluteDielectric=1.0, solventDielectric=78.5,
                     removeCMMotion=True, hydrogenMass=None, ewaldErrorTolerance=0.0005,
                     switchDistance=0.0*u.nanometer, gbsaModel='ACE'):
        """Construct an OpenMM System representing the topology described by this
        prmtop file.

        Parameters
        ----------
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, or LJPME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        implicitSolvent : object=None
            If not None, the implicit solvent model to use.  Allowed values are
            HCT, OBC1, OBC2, GBn, or GBn2.
        implicitSolventSaltConc : float=0.0*unit.moles/unit.liter
            The salt concentration for GB calculations (modelled as a debye
            screening parameter). It is converted to the debye length (kappa)
            using the provided temperature and solventDielectric
        temperature : float=300*kelvin
            Temperature of the system. Only used to compute the Debye length
            from implicitSolventSoltConc
        implicitSolventKappa : float units of 1/length
            If this value is set, implicitSolventSaltConc will be ignored.
        soluteDielectric : float=1.0
            The solute dielectric constant to use in the implicit solvent model.
        solventDielectric : float=78.5
            The solvent dielectric constant to use in the implicit solvent
            model.
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep their
            total mass the same.
        ewaldErrorTolerance : float=0.0005
            The error tolerance to use if nonbondedMethod is Ewald, PME, or LJPME.
        switchDistance : float=0*nanometers
            The distance at which the potential energy switching function is
            turned on for Lennard-Jones interactions. If the switchDistance is 0
            or evaluates to boolean False, no switching function will be used.
            Values greater than nonbondedCutoff or less than 0 raise ValueError
        gbsaModel : str='ACE'
            The SA model used to model the nonpolar solvation component of GB
            implicit solvent models. If GB is active, this must be 'ACE' or None
            (the latter indicates no SA model will be used). Other values will
            result in a ValueError

        Returns
        -------
        System
            the newly created System
        """
        if self._prmtop.chamber:
            raise ValueError("CHAMBER-style topology file detected. CHAMBER "
                             "topologies are not supported -- use the native "
                             "CHARMM files directly.")
        methodMap = {ff.NoCutoff:'NoCutoff',
                     ff.CutoffNonPeriodic:'CutoffNonPeriodic',
                     ff.CutoffPeriodic:'CutoffPeriodic',
                     ff.Ewald:'Ewald',
                     ff.PME:'PME',
                     ff.LJPME:'LJPME'}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal value for nonbonded method')
        if not self._prmtop.getIfBox() and nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME, ff.LJPME):
            raise ValueError('Illegal nonbonded method for a non-periodic system')
        constraintMap = {None:None,
                         ff.HBonds:'h-bonds',
                         ff.AllBonds:'all-bonds',
                         ff.HAngles:'h-angles'}
        if constraints is None:
            constraintString = None
        elif constraints in constraintMap:
            constraintString = constraintMap[constraints]
        else:
            raise ValueError('Illegal value for constraints')
        if implicitSolvent is None:
            implicitString = None
        elif implicitSolvent is HCT:
            implicitString = 'HCT'
        elif implicitSolvent is OBC1:
            implicitString = 'OBC1'
        elif implicitSolvent is OBC2:
            implicitString = 'OBC2'
        elif implicitSolvent is GBn:
            implicitString = 'GBn'
        elif implicitSolvent is GBn2:
            implicitString = 'GBn2'
        else:
            raise ValueError('Illegal value for implicit solvent model')
        # If implicitSolventKappa is None, compute it from the salt concentration
        if implicitSolvent is not None and implicitSolventKappa is None:
            if u.is_quantity(implicitSolventSaltConc):
                implicitSolventSaltConc = implicitSolventSaltConc.value_in_unit(u.moles/u.liter)
            if u.is_quantity(temperature):
                temperature = temperature.value_in_unit(u.kelvin)
            # The constant is 1 / sqrt( epsilon_0 * kB / (2 * NA * q^2 * 1000) )
            # where NA is avogadro's number, epsilon_0 is the permittivity of
            # free space, q is the elementary charge (this number matches
            # Amber's kappa conversion factor)
            implicitSolventKappa = 50.33355 * sqrt(implicitSolventSaltConc / solventDielectric / temperature)
            # Multiply by 0.73 to account for ion exclusions, and multiply by 10
            # to convert to 1/nm from 1/angstroms
            implicitSolventKappa *= 7.3
        elif implicitSolvent is None:
            implicitSolventKappa = 0.0

        sys = amber_file_parser.readAmberSystem(self.topology, prmtop_loader=self._prmtop, shake=constraintString,
                        nonbondedCutoff=nonbondedCutoff, nonbondedMethod=methodMap[nonbondedMethod],
                        flexibleConstraints=False, gbmodel=implicitString, soluteDielectric=soluteDielectric,
                        solventDielectric=solventDielectric, implicitSolventKappa=implicitSolventKappa,
                        rigidWater=rigidWater, elements=self.elements, gbsaModel=gbsaModel)

        if hydrogenMass is not None:
            for atom1, atom2 in self.topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)
        for force in sys.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setEwaldErrorTolerance(ewaldErrorTolerance)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        if switchDistance and nonbondedMethod is not ff.NoCutoff:
            # make sure it's legal
            if (_strip_optunit(switchDistance, u.nanometer) >=
                    _strip_optunit(nonbondedCutoff, u.nanometer)):
                raise ValueError('switchDistance is too large compared '
                                 'to the cutoff!')
            if _strip_optunit(switchDistance, u.nanometer) < 0:
                # Detects negatives for both Quantity and float
                raise ValueError('switchDistance must be non-negative!')
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(switchDistance)

        return sys
