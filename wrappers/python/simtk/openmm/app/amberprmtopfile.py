"""
armberprmtopfile.py: Used for loading AMBER prmtop files.
"""
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm.app import Topology
from simtk.openmm.app.internal import amber_file_parser
import forcefield as ff
import element as elem
import simtk.unit as unit

# Enumerated values for implicit solvent model

OBC = object()

class AmberPrmtopFile(object):
    """AmberPrmtopFile parses an AMBER prmtop file and constructs a Topology and (optionally) an OpenMM System from it."""
    
    def __init__(self, file):
        """Load a prmtop file."""
        top = Topology()
        self.topology = top
        
        # Load the prmtop file
        
        prmtop = amber_file_parser.PrmtopLoader(file)
        self.prmtop = prmtop

        # Add atoms to the topology

        lastResidue = None
        c = top.addChain()
        for index in range(prmtop.getNumAtoms()):
            resNumber = prmtop.getResidueNumber(index)
            if resNumber != lastResidue:
                lastResidue = resNumber
                resName = prmtop.getResidueLabel(iAtom=index)
                r = top.addResidue(resName, c)
            atomName = prmtop.getAtomName(index)

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
                    pass
            top.addAtom(atomName, element, r)
        
        # Add bonds to the topology
        
        for bond in prmtop.getBondsWithH():
            top.addBond(bond[0], bond[1])
        for bond in prmtop.getBondsNoH():
            top.addBond(bond[0], bond[1])

    def createSystem(self, nonbondedMethod=ff.NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=True, implicitSolvent=None):
        """Construct an OpenMM System representing the topology described by this prmtop file.
        
        Parameters:
         - nonbondedMethod (object=NoCutoff) The method to use for nonbonded interactions.  Allowed values are
           NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
         - nonbondedCutoff (distance=1*nanometer) The cutoff distance to use for nonbonded interactions
         - constraints (object=None) Specifies which bonds angles should be implemented with constraints.
           Allowed values are None, HBonds, AllBonds, or HAngles.
         - rigidWater (boolean=True) If true, water molecules will be fully rigid regardless of the value passed for the constraints argument
         - implicitSolvent (object=None) If not None, the implicit solvent model to use
        Returns: the newly created System
        """
        methodMap = {ff.NoCutoff:'NoCutoff',
                     ff.CutoffNonPeriodic:'CutoffNonPeriodic',
                     ff.CutoffPeriodic:'CutoffPeriodic',
                     ff.Ewald:'Ewald',
                     ff.PME:'PME'}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal value for nonbonded method')
        if not self.prmtop.getIfBox() and nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME):
            raise ValueError('Illegal nonbonded method for a non-periodic system')
        if nonbondedMethod == ff.NoCutoff:
            nonbondedCutoff = None
        constraintMap = {None:None,
                         ff.HBonds:'h-bonds',
                         ff.AllBonds:'all-bonds',
                         ff.HAngles:'h-angles'}
        if constraints == ff.NoConstraints:
            constraintString = None
        elif constraints in constraintMap:
            constraintString = constraintMap[constraints]
        else:
            raise ValueError('Illegal value for constraints')
        if implicitSolvent is None:
            implicitString = None
        elif implicitSolvent == OBC:
            implicitString = 'OBC'
        else:
            raise ValueError('Illegal value for implicit solvent model')
        return amber_file_parser.readAmberSystem(prmtop_loader=self.prmtop, shake=constraintString, nonbondedCutoff=nonbondedCutoff,
                                                 nonbondedMethod=methodMap[nonbondedMethod], flexibleConstraints=False, gbmodel=implicitString)