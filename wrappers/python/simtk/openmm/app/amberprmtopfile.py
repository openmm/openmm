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
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm.app import Topology
from simtk.openmm.app import PDBFile
from simtk.openmm.app.internal import amber_file_parser
import forcefield as ff
import element as elem
import simtk.unit as unit
import simtk.openmm as mm

# Enumerated values for implicit solvent model

HCT = object()
OBC1 = object()
OBC2 = object()
GBn = object()

class AmberPrmtopFile(object):
    """AmberPrmtopFile parses an AMBER prmtop file and constructs a Topology and (optionally) an OpenMM System from it."""
    
    def __init__(self, file):
        """Load a prmtop file."""
        top = Topology()
        ## The Topology read from the prmtop file
        self.topology = top
        
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

            # Try to guess the element.
            
            upper = atomName.upper()
            if upper.startswith('CL'):
                element = elem.chlorine
            elif upper.startswith('NA'):
                element = elem.sodium
            elif upper.startswith('MG'):
                element = elem.magnesium
            else:
                try:
                    element = elem.get_by_symbol(atomName[0])
                except KeyError:
                    element = None
            top.addAtom(atomName, element, r)
        
        # Add bonds to the topology
        
        atoms = list(top.atoms())
        for bond in prmtop.getBondsWithH():
            top.addBond(atoms[bond[0]], atoms[bond[1]])
        for bond in prmtop.getBondsNoH():
            top.addBond(atoms[bond[0]], atoms[bond[1]])
        
        # Set the periodic box size.
        
        if prmtop.getIfBox():
            top.setUnitCellDimensions(tuple(x.value_in_unit(unit.nanometer) for x in prmtop.getBoxBetaAndDimensions()[1:4])*unit.nanometer)

    def createSystem(self, nonbondedMethod=ff.NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=True, implicitSolvent=None, soluteDielectric=1.0, solventDielectric=78.5, removeCMMotion=True,
                     ewaldErrorTolerance=0.0005):
        """Construct an OpenMM System representing the topology described by this prmtop file.
        
        Parameters:
         - nonbondedMethod (object=NoCutoff) The method to use for nonbonded interactions.  Allowed values are
           NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
         - nonbondedCutoff (distance=1*nanometer) The cutoff distance to use for nonbonded interactions
         - constraints (object=None) Specifies which bonds angles should be implemented with constraints.
           Allowed values are None, HBonds, AllBonds, or HAngles.
         - rigidWater (boolean=True) If true, water molecules will be fully rigid regardless of the value passed for the constraints argument
         - implicitSolvent (object=None) If not None, the implicit solvent model to use.  Allowed values are HCT, OBC1, OBC2, or GBn.
         - soluteDielectric (float=1.0) The solute dielectric constant to use in the implicit solvent model.
         - solventDielectric (float=78.5) The solvent dielectric constant to use in the implicit solvent model.
         - removeCMMotion (boolean=True) If true, a CMMotionRemover will be added to the System
         - ewaldErrorTolerance (float=0.0005) The error tolerance to use if nonbondedMethod is Ewald or PME.
        Returns: the newly created System
        """
        methodMap = {ff.NoCutoff:'NoCutoff',
                     ff.CutoffNonPeriodic:'CutoffNonPeriodic',
                     ff.CutoffPeriodic:'CutoffPeriodic',
                     ff.Ewald:'Ewald',
                     ff.PME:'PME'}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal value for nonbonded method')
        if not self._prmtop.getIfBox() and nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME):
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
        elif implicitSolvent == HCT:
            implicitString = 'HCT'
        elif implicitSolvent == OBC1:
            implicitString = 'OBC1'
        elif implicitSolvent == OBC2:
            implicitString = 'OBC2'
        elif implicitSolvent == GBn:
            implicitString = 'GBn'
        else:
            raise ValueError('Illegal value for implicit solvent model')
        sys = amber_file_parser.readAmberSystem(prmtop_loader=self._prmtop, shake=constraintString, nonbondedCutoff=nonbondedCutoff,
                                                 nonbondedMethod=methodMap[nonbondedMethod], flexibleConstraints=False, gbmodel=implicitString,
                                                 soluteDielectric=soluteDielectric, solventDielectric=solventDielectric, rigidWater=rigidWater)
        for force in sys.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setEwaldErrorTolerance(ewaldErrorTolerance)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())
        return sys