"""
gromacstopfile.py: Used for loading Gromacs top files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2013 Stanford University and the Authors.
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
import forcefield as ff
import element as elem
import amberprmtopfile as prmtop
import simtk.unit as unit
import simtk.openmm as mm
import math
import os

HBonds = ff.HBonds
AllBonds = ff.AllBonds
HAngles = ff.HAngles

OBC2 = prmtop.OBC2

class GromacsTopFile(object):
    """GromacsTopFile parses a Gromacs top file and constructs a Topology and (optionally) an OpenMM System from it."""
    
    class _MoleculeType(object):
        """Inner class to store information about a molecule type."""
        def __init__(self):
            self.atoms = []
            self.bonds = []
            self.angles = []
            self.dihedrals = []
            self.exclusions = []
            self.pairs = []
            self.cmaps = []
    
    def _processFile(self, file):
        append = ''
        for line in open(file):
            if line.strip().endswith('\\'):
                append = '%s %s' % (append, line[:line.rfind('\\')])
            else:
                self._processLine(append+' '+line, file)
                append = ''
    
    def _processLine(self, line, file):
        """Process one line from a file."""
        if ';' in line:
            line = line[:line.index(';')]
        stripped = line.strip()
        ignore = not all(self._ifStack)
        if stripped.startswith('*') or len(stripped) == 0:
            # A comment or empty line.
            return
        
        elif stripped.startswith('[') and not ignore:
            # The start of a category.
            if not stripped.endswith(']'):
                raise ValueError('Illegal line in .top file: '+line)
            self._currentCategory = stripped[1:-1].strip()
                
        elif stripped.startswith('#'):
            # A preprocessor command.
            fields = stripped.split()
            command = fields[0]
            if command == '#include' and not ignore:
                # Locate the file to include
                name = stripped[len(command):].strip(' \t"<>')
                searchDirs = self._includeDirs+(os.path.dirname(file),)
                for dir in searchDirs:
                    file = os.path.join(dir, name)
                    if os.path.isfile(file):
                        # We found the file, so process it.
                        self._processFile(file)
                        break
                else:
                    raise ValueError('Could not locate #include file: '+name)
            elif command == '#define' and not ignore:
                # Add a value to our list of defines.
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                name = fields[1]
                valueStart = stripped.find(name, len(command))+len(name)
                value = line[valueStart:].strip()
                self._defines[name] = value
            elif command == '#ifdef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                name = fields[1]
                self._ifStack.append(name in self._defines)
            elif command == '#ifndef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                name = fields[1]
                self._ifStack.append(name not in self._defines)
            elif command == '#endif':
                # Pop an entry off the if stack.
                if len(self._ifStack) == 0:
                    raise ValueError('Unexpected line in .top file: '+line)
                del(self._ifStack[-1])
                
        elif not ignore:
            # A line of data for the current category
            if self._currentCategory is None:
                raise ValueError('Unexpected line in .top file: '+line)
            if self._currentCategory == 'defaults':
                self._processDefaults(line)
            elif self._currentCategory == 'moleculetype':
                self._processMoleculeType(line)
            elif self._currentCategory == 'molecules':
                self._processMolecule(line)
            elif self._currentCategory == 'atoms':
                self._processAtom(line)
            elif self._currentCategory == 'bonds':
                self._processBond(line)
            elif self._currentCategory == 'angles':
                self._processAngle(line)
            elif self._currentCategory == 'dihedrals':
                self._processDihedral(line)
            elif self._currentCategory == 'exclusions':
                self._processExclusion(line)
            elif self._currentCategory == 'pairs':
                self._processPair(line)
            elif self._currentCategory == 'cmap':
                self._processCmap(line)
            elif self._currentCategory == 'atomtypes':
                self._processAtomType(line)
            elif self._currentCategory == 'bondtypes':
                self._processBondType(line)
            elif self._currentCategory == 'angletypes':
                self._processAngleType(line)
            elif self._currentCategory == 'dihedraltypes':
                self._processDihedralType(line)
            elif self._currentCategory == 'implicit_genborn_params':
                self._processImplicitType(line)
            elif self._currentCategory == 'pairtypes':
                self._processPairType(line)
            elif self._currentCategory == 'cmaptypes':
                self._processCmapType(line)
    
    def _processDefaults(self, line):
        """Process the [ defaults ] line."""
        fields = line.split()
        if len(fields) < 4:
            raise ValueError('Too few fields in [ defaults ] line: '+line);
        if fields[0] != '1':
            raise ValueError('Unsupported nonbonded type: '+fields[0])
        if fields[1] != '2':
            raise ValueError('Unsupported combination rule: '+fields[1])
        if fields[2].lower() == 'no':
            raise ValueError('gen_pairs=no is not supported')
        self._defaults = fields
    
    def _processMoleculeType(self, line):
        """Process a line in the [ moleculetypes ] category."""
        fields = line.split()
        if len(fields) < 1:
            raise ValueError('Too few fields in [ moleculetypes ] line: '+line);
        type = GromacsTopFile._MoleculeType()
        self._moleculeTypes[fields[0]] = type
        self._currentMoleculeType = type
    
    def _processMolecule(self, line):
        """Process a line in the [ molecules ] category."""
        fields = line.split()
        if len(fields) < 2:
            raise ValueError('Too few fields in [ molecules ] line: '+line);
        self._molecules.append((fields[0], int(fields[1])))
    
    def _processAtom(self, line):
        """Process a line in the [ atoms ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ atoms ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ atoms ] line: '+line);
        self._currentMoleculeType.atoms.append(fields)
    
    def _processBond(self, line):
        """Process a line in the [ bonds ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ bonds ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 3:
            raise ValueError('Too few fields in [ bonds ] line: '+line);
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ bonds ] line: '+line);
        self._currentMoleculeType.bonds.append(fields)
    
    def _processAngle(self, line):
        """Process a line in the [ angles ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ angles ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 4:
            raise ValueError('Too few fields in [ angles ] line: '+line);
        if fields[3] not in ('1', '5'):
            raise ValueError('Unsupported function type in [ angles ] line: '+line);
        self._currentMoleculeType.angles.append(fields)
    
    def _processDihedral(self, line):
        """Process a line in the [ dihedrals ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ dihedrals ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ dihedrals ] line: '+line);
        if fields[4] not in ('1', '2', '3', '4', '9'):
            raise ValueError('Unsupported function type in [ dihedrals ] line: '+line);
        self._currentMoleculeType.dihedrals.append(fields)
    
    def _processExclusion(self, line):
        """Process a line in the [ exclusions ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ exclusions ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 2:
            raise ValueError('Too few fields in [ exclusions ] line: '+line);
        self._currentMoleculeType.exclusions.append(fields)
    
    def _processPair(self, line):
        """Process a line in the [ pairs ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ pairs ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 3:
            raise ValueError('Too few fields in [ pairs ] line: '+line);
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ pairs ] line: '+line);
        self._currentMoleculeType.pairs.append(fields)
    
    def _processCmap(self, line):
        """Process a line in the [ cmaps ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ cmap ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 6:
            raise ValueError('Too few fields in [ pairs ] line: '+line);
        self._currentMoleculeType.cmaps.append(fields)
    
    def _processAtomType(self, line):
        """Process a line in the [ atomtypes ] category."""
        fields = line.split()
        if len(fields) < 7:
            raise ValueError('Too few fields in [ atomtypes ] line: '+line);
        self._atomTypes[fields[0]] = fields
    
    def _processBondType(self, line):
        """Process a line in the [ bondtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ bondtypes ] line: '+line);
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ bondtypes ] line: '+line);
        self._bondTypes[tuple(fields[:2])] = fields
    
    def _processAngleType(self, line):
        """Process a line in the [ angletypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError('Too few fields in [ angletypes ] line: '+line);
        if fields[3] not in ('1', '5'):
            raise ValueError('Unsupported function type in [ angletypes ] line: '+line);
        self._angleTypes[tuple(fields[:3])] = fields
    
    def _processDihedralType(self, line):
        """Process a line in the [ dihedraltypes ] category."""
        fields = line.split()
        if len(fields) < 7:
            raise ValueError('Too few fields in [ dihedraltypes ] line: '+line);
        if fields[4] not in ('1', '2', '3', '4', '9'):
            raise ValueError('Unsupported function type in [ dihedraltypes ] line: '+line);
        key = tuple(fields[:5])
        if fields[4] == '9' and key in self._dihedralTypes:
            # There are multiple dihedrals defined for these atom types.
            self._dihedralTypes[key].append(fields)
        else:
            self._dihedralTypes[key] = [fields]
    
    def _processImplicitType(self, line):
        """Process a line in the [ implicit_genborn_params ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError('Too few fields in [ implicit_genborn_params ] line: '+line);
        self._implicitTypes[fields[0]] = fields
    
    def _processPairType(self, line):
        """Process a line in the [ pairtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ pairtypes] line: '+line);
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ pairtypes ] line: '+line);
        self._pairTypes[tuple(fields[:2])] = fields
    
    def _processCmapType(self, line):
        """Process a line in the [ cmaptypes ] category."""
        fields = line.split()
        if len(fields) < 8 or len(fields) < 8+int(fields[6])*int(fields[7]):
            raise ValueError('Too few fields in [ cmaptypes ] line: '+line);
        if fields[5] != '1':
            raise ValueError('Unsupported function type in [ cmaptypes ] line: '+line);
        self._cmapTypes[tuple(fields[:5])] = fields
    
    def __init__(self, file, unitCellDimensions=None, includeDir='/usr/local/gromacs/share/gromacs/top', defines={}):
        """Load a top file.
        
        Parameters:
         - file (string) the name of the file to load
         - unitCellDimensions (Vec3=None) the dimensions of the crystallographic unit cell
         - includeDir (string=/usr/local/gromacs/share/gromacs/top) a directory in which to look for other files
           included from the top file
         - defines (map={}) preprocessor definitions that should be predefined when parsing the file
         """
        self._includeDirs = (os.path.dirname(file), includeDir)
        self._defines = defines
        
        # Parse the file.
        
        self._currentCategory = None
        self._ifStack = []
        self._moleculeTypes = {}
        self._molecules = []
        self._currentMoleculeType = None
        self._atomTypes = {}
        self._bondTypes= {}
        self._angleTypes = {}
        self._dihedralTypes = {}
        self._implicitTypes = {}
        self._pairTypes = {}
        self._cmapTypes = {}
        self._processFile(file)
        
        # Create the Topology from it.
        
        top = Topology()
        ## The Topology read from the prmtop file
        self.topology = top
        top.setUnitCellDimensions(unitCellDimensions)
        PDBFile._loadNameReplacementTables()
        for moleculeName, moleculeCount in self._molecules:
            if moleculeName not in self._moleculeTypes:
                raise ValueError("Unknown molecule type: "+moleculeName)
            moleculeType = self._moleculeTypes[moleculeName]
            
            # Create the specified number of molecules of this type.
            
            for i in range(moleculeCount):
                atoms = []
                lastResidue = None
                c = top.addChain()
                for index, fields in enumerate(moleculeType.atoms):
                    resNumber = fields[2]
                    if resNumber != lastResidue:
                        lastResidue = resNumber
                        resName = fields[3]
                        if resName in PDBFile._residueNameReplacements:
                            resName = PDBFile._residueNameReplacements[resName]
                        r = top.addResidue(resName, c)
                        if resName in PDBFile._atomNameReplacements:
                            atomReplacements = PDBFile._atomNameReplacements[resName]
                        else:
                            atomReplacements = {}
                    atomName = fields[4]
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
                    atoms.append(top.addAtom(atomName, element, r))
                
                # Add bonds to the topology
                
                for fields in moleculeType.bonds:
                    top.addBond(atoms[int(fields[0])-1], atoms[int(fields[1])-1])

    def createSystem(self, nonbondedMethod=ff.NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=True, implicitSolvent=None, soluteDielectric=1.0, solventDielectric=78.5, ewaldErrorTolerance=0.0005, removeCMMotion=True):
        """Construct an OpenMM System representing the topology described by this prmtop file.
        
        Parameters:
         - nonbondedMethod (object=NoCutoff) The method to use for nonbonded interactions.  Allowed values are
           NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
         - nonbondedCutoff (distance=1*nanometer) The cutoff distance to use for nonbonded interactions
         - constraints (object=None) Specifies which bonds and angles should be implemented with constraints.
           Allowed values are None, HBonds, AllBonds, or HAngles.
         - rigidWater (boolean=True) If true, water molecules will be fully rigid regardless of the value passed for the constraints argument
         - implicitSolvent (object=None) If not None, the implicit solvent model to use.  The only allowed value is OBC2.
         - soluteDielectric (float=1.0) The solute dielectric constant to use in the implicit solvent model.
         - solventDielectric (float=78.5) The solvent dielectric constant to use in the implicit solvent model.
         - ewaldErrorTolerance (float=0.0005) The error tolerance to use if nonbondedMethod is Ewald or PME.
         - removeCMMotion (boolean=True) If true, a CMMotionRemover will be added to the System
        Returns: the newly created System
        """
        # Create the System.
        
        sys = mm.System()
        boxSize = self.topology.getUnitCellDimensions()
        if boxSize is not None:
            sys.setDefaultPeriodicBoxVectors((boxSize[0], 0, 0), (0, boxSize[1], 0), (0, 0, boxSize[2]))
        elif nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME):
            raise ValueError('Illegal nonbonded method for a non-periodic system')
        nb = mm.NonbondedForce()
        sys.addForce(nb)
        if implicitSolvent is OBC2:
            gb = mm.GBSAOBCForce()
            gb.setSoluteDielectric(soluteDielectric)
            gb.setSolventDielectric(solventDielectric)
            sys.addForce(gb)
            nb.setReactionFieldDielectric(1.0)
        elif implicitSolvent is not None:
            raise ValueError('Illegal value for implicitSolvent')
        bonds = None
        angles = None
        periodic = None
        rb = None
        harmonicTorsion = None
        cmap = None
        mapIndices = {}
        bondIndices = []
        topologyAtoms = list(self.topology.atoms())
        exceptions = []
        fudgeQQ = float(self._defaults[4])
        
        # Loop over molecules and create the specified number of each type.
        
        for moleculeName, moleculeCount in self._molecules:
            moleculeType = self._moleculeTypes[moleculeName]
            for i in range(moleculeCount):
                
                # Record the types of all atoms.
                
                baseAtomIndex = sys.getNumParticles()
                atomTypes = [atom[1] for atom in moleculeType.atoms]
                try:
                    [self._atomTypes[t][1] for t in atomTypes]
                except KeyError as e:
                    raise ValueError('Unknown atom type: '+e.message)
                
                # Add atoms.
                
                for fields in moleculeType.atoms:
                    if len(fields) >= 8:
                        mass = float(fields[7])
                    else:
                        mass = float(self._atomTypes[fields[1]][2])
                    sys.addParticle(mass)
                
                # Add bonds.
                
                atomBonds = [{} for x in range(len(moleculeType.atoms))]
                for fields in moleculeType.bonds:
                    atoms = [int(x)-1 for x in fields[:2]]
                    types = tuple(atomTypes[i] for i in atoms)
                    if len(fields) >= 5:
                        params = fields[3:5]
                    elif types in self._bondTypes:
                        params = self._bondTypes[types][3:5]
                    elif types[::-1] in self._bondTypes:
                        params = self._bondTypes[types[::-1]][3:5]
                    else:
                        raise ValueError('No parameters specified for bond: '+fields[0]+', '+fields[1])
                    # Decide whether to use a constraint or a bond.
                    useConstraint = False
                    if rigidWater and topologyAtoms[baseAtomIndex+atoms[0]].residue.name == 'HOH':
                        useConstraint = True
                    if constraints in (AllBonds, HAngles):
                        useConstraint = True
                    elif constraints is HBonds:
                        elements = [topologyAtoms[baseAtomIndex+i].element for i in atoms]
                        if elem.hydrogen in elements:
                            useConstraint = True
                    # Add the bond or constraint.
                    length = float(params[0])
                    if useConstraint:
                        sys.addConstraint(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], length)
                    else:
                        if bonds is None:
                            bonds = mm.HarmonicBondForce()
                            sys.addForce(bonds)
                        bonds.addBond(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], length, float(params[1]))
                    # Record information that will be needed for constraining angles.
                    atomBonds[atoms[0]][atoms[1]] = length
                    atomBonds[atoms[1]][atoms[0]] = length
                
                # Add angles.
                
                degToRad = math.pi/180
                for fields in moleculeType.angles:
                    atoms = [int(x)-1 for x in fields[:3]]
                    types = tuple(atomTypes[i] for i in atoms)
                    if len(fields) >= 6:
                        params = fields[4:]
                    elif types in self._angleTypes:
                        params = self._angleTypes[types][4:]
                    elif types[::-1] in self._angleTypes:
                        params = self._angleTypes[types[::-1]][4:]
                    else:
                        raise ValueError('No parameters specified for angle: '+fields[0]+', '+fields[1]+', '+fields[2])
                    # Decide whether to use a constraint or a bond.
                    useConstraint = False
                    if rigidWater and topologyAtoms[baseAtomIndex+atoms[0]].residue.name == 'HOH':
                        useConstraint = True
                    if constraints is HAngles:
                        elements = [topologyAtoms[baseAtomIndex+i].element for i in atoms]
                        if elements[0] == elem.hydrogen and elements[2] == elem.hydrogen:
                            useConstraint = True
                        elif elements[1] == elem.oxygen and (elements[0] == elem.hydrogen or elements[2] == elem.hydrogen):
                            useConstraint = True
                    # Add the bond or constraint.
                    theta = float(params[0])*degToRad
                    if useConstraint:
                        # Compute the distance between atoms and add a constraint
                        if atoms[0] in atomBonds[atoms[1]] and atoms[2] in atomBonds[atoms[1]]:
                            l1 = atomBonds[atoms[1]][atoms[0]]
                            l2 = atomBonds[atoms[1]][atoms[2]]
                            length = math.sqrt(l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta))
                            sys.addConstraint(baseAtomIndex+atoms[0], baseAtomIndex+atoms[2], length)
                    else:
                        if angles is None:
                            angles = mm.HarmonicAngleForce()
                            sys.addForce(angles)
                        angles.addAngle(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], theta, float(params[1]))
                        if fields[3] == '5':
                            # This is a Urey-Bradley term, so add the bond.
                            if bonds is None:
                                bonds = mm.HarmonicBondForce()
                                sys.addForce(bonds)
                            k = float(params[3])
                            if k != 0:
                                bonds.addBond(baseAtomIndex+atoms[0], baseAtomIndex+atoms[2], float(params[2]), k)

                # Add torsions.
                
                for fields in moleculeType.dihedrals:
                    atoms = [int(x)-1 for x in fields[:4]]
                    types = tuple(atomTypes[i] for i in atoms)
                    dihedralType = fields[4]
                    reversedTypes = types[::-1]+(dihedralType,)
                    types = types+(dihedralType,)
                    if (dihedralType in ('1', '2', '4', '9') and len(fields) > 7) or (dihedralType == '3' and len(fields) > 10):
                        paramsList = [fields]
                    else:
                        # Look for a matching dihedral type.
                        paramsList = None
                        for key in self._dihedralTypes:
                            if all(a == b or a == 'X' for a, b in zip(key, types)) or all(a == b or a == 'X' for a, b in zip(key, reversedTypes)):
                                paramsList = self._dihedralTypes[key]
                                if 'X' not in key:
                                    break
                        if paramsList is None:
                            raise ValueError('No parameters specified for dihedral: '+fields[0]+', '+fields[1]+', '+fields[2]+', '+fields[3])
                    for params in paramsList:
                        if dihedralType in ('1', '4', '9'):
                            # Periodic torsion
                            k = float(params[6])
                            if k != 0:
                                if periodic is None:
                                    periodic = mm.PeriodicTorsionForce()
                                    sys.addForce(periodic)
                                periodic.addTorsion(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], int(params[7]), float(params[5])*degToRad, k)
                        elif dihedralType == '2':
                            # Harmonic torsion
                            k = float(params[6])
                            if k != 0:
                                if harmonicTorsion is None:
                                    harmonicTorsion = mm.CustomTorsionForce('0.5*k*(theta-theta0)^2')
                                    harmonicTorsion.addPerTorsionParameter('theta0')
                                    harmonicTorsion.addPerTorsionParameter('k')
                                    sys.addForce(harmonicTorsion)
                                harmonicTorsion.addTorsion(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], (float(params[5])*degToRad, k))
                        else:
                            # RB Torsion
                            c = [float(x) for x in params[5:11]]
                            if any(x != 0 for x in c):
                                if rb is None:
                                    rb = mm.RBTorsionForce()
                                    sys.addForce(rb)
                                rb.addTorsion(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], c[0], c[1], c[2], c[3], c[4], c[5])
                
                # Add CMAP terms.

                for fields in moleculeType.cmaps:
                    atoms = [int(x)-1 for x in fields[:5]]
                    types = tuple(atomTypes[i] for i in atoms)
                    if len(fields) >= 8 and len(fields) >= 8+int(fields[6])*int(fields[7]):
                        params = fields
                    elif types in self._cmapTypes:
                        params = self._cmapTypes[types]
                    elif types[::-1] in self._cmapTypes:
                        params = self._cmapTypes[types[::-1]]
                    else:
                        raise ValueError('No parameters specified for cmap: '+fields[0]+', '+fields[1]+', '+fields[2]+', '+fields[3]+', '+fields[4])
                    if cmap is None:
                        cmap = mm.CMAPTorsionForce()
                        sys.addForce(cmap)
                    mapSize = int(params[6])
                    if mapSize != int(params[7]):
                        raise ValueError('Non-square CMAPs are not supported')
                    map = []
                    for i in range(mapSize):
                        for j in range(mapSize):
                            map.append(float(params[8+mapSize*((j+mapSize/2)%mapSize)+((i+mapSize/2)%mapSize)]))
                    map = tuple(map)
                    if map not in mapIndices:
                        mapIndices[map] = cmap.addMap(mapSize, map)
                    cmap.addTorsion(mapIndices[map], baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3],
                                 baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], baseAtomIndex+atoms[4])

                # Set nonbonded parameters for particles.
                
                for fields in moleculeType.atoms:
                    params = self._atomTypes[fields[1]]
                    if len(fields) > 6:
                        q = float(fields[6])
                    else:
                        q = float(params[3])
                    nb.addParticle(q, float(params[5]), float(params[6]))
                    if implicitSolvent is OBC2:
                        if fields[1] not in self._implicitTypes:
                            raise ValueError('No implicit solvent parameters specified for atom type: '+fields[1])
                        gbparams = self._implicitTypes[fields[1]]
                        gb.addParticle(q, float(gbparams[4]), float(gbparams[5]))
                for fields in moleculeType.bonds:
                    atoms = [int(x)-1 for x in fields[:2]]
                    bondIndices.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1]))
                
                # Record nonbonded exceptions.
                
                for fields in moleculeType.pairs:
                    atoms = [int(x)-1 for x in fields[:2]]
                    types = tuple(atomTypes[i] for i in atoms)
                    if len(fields) >= 5:
                        params = fields[3:5]
                    elif types in self._pairTypes:
                        params = self._pairTypes[types][3:5]
                    elif types[::-1] in self._pairTypes:
                        params = self._pairTypes[types[::-1]][3:5]
                    else:
                        continue # We'll use the automatically generated parameters
                    atom1params = nb.getParticleParameters(baseAtomIndex+atoms[0])
                    atom2params = nb.getParticleParameters(baseAtomIndex+atoms[1])
                    exceptions.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], atom1params[0]*atom2params[0]*fudgeQQ, params[0], params[1]))
        
        # Create nonbonded exceptions.
        
        nb.createExceptionsFromBonds(bondIndices, fudgeQQ, float(self._defaults[3]))
        for exception in exceptions:
            nb.addException(exception[0], exception[1], exception[2], float(exception[3]), float(exception[4]), True)

        # Finish configuring the NonbondedForce.
        
        methodMap = {ff.NoCutoff:mm.NonbondedForce.NoCutoff,
                     ff.CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     ff.CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic,
                     ff.Ewald:mm.NonbondedForce.Ewald,
                     ff.PME:mm.NonbondedForce.PME}
        nb.setNonbondedMethod(methodMap[nonbondedMethod])
        nb.setCutoffDistance(nonbondedCutoff)
        nb.setEwaldErrorTolerance(ewaldErrorTolerance)

        # Add a CMMotionRemover.

        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())
        return sys
