"""
gromacstopfile.py: Used for loading Gromacs top files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2024 Stanford University and the Authors.
Authors: Peter Eastman
Contributors: Jason Swails

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

from openmm.app import Topology
from openmm.app import PDBFile
from . import forcefield as ff
from . import element as elem
from . import amberprmtopfile as prmtop
import openmm.unit as unit
import openmm as mm
import math
import os
import re
import shutil
from collections import OrderedDict, defaultdict
from itertools import combinations, combinations_with_replacement
from copy import deepcopy

HBonds = ff.HBonds
AllBonds = ff.AllBonds
HAngles = ff.HAngles

OBC2 = prmtop.OBC2

novarcharre = re.compile(r'\W')

def _find_all_instances_in_string(string, substr):
    """ Find indices of all instances of substr in string """
    indices = []
    idx = string.find(substr, 0)
    while idx > -1:
        indices.append(idx)
        idx = string.find(substr, idx+1)
    return indices

def _replace_defines(line, defines):
    """ Replaces defined tokens in a given line """
    if not defines: return line
    for define in reversed(defines):
        value = defines[define]
        indices = _find_all_instances_in_string(line, define)
        if not indices: continue
        # Check to see if it's inside of quotes
        inside = ''
        idx = 0
        n_to_skip = 0
        new_line = []
        for i, char in enumerate(line):
            if n_to_skip:
                n_to_skip -= 1
                continue
            if char in ('\'"'):
                if not inside:
                    inside = char
                else:
                    if inside == char:
                        inside = ''
            if idx < len(indices) and i == indices[idx]:
                if inside:
                    new_line.append(char)
                    idx += 1
                    continue
                if i == 0 or novarcharre.match(line[i-1]):
                    endidx = indices[idx] + len(define)
                    if endidx >= len(line) or novarcharre.match(line[endidx]):
                        new_line.extend(list(value))
                        n_to_skip = len(define) - 1
                        idx += 1
                        continue
                idx += 1
            new_line.append(char)
        line = ''.join(new_line)

    return line

class GromacsTopFile(object):
    """GromacsTopFile parses a Gromacs top file and constructs a Topology and (optionally) an OpenMM System from it."""

    class _MoleculeType(object):
        """Inner class to store information about a molecule type."""
        def __init__(self, name, nrexcl):
            self.name = name
            self.nrexcl = nrexcl
            self.atoms = []
            self.bonds = []
            self.angles = []
            self.dihedrals = []
            self.exclusions = []
            self.pairs = []
            self.constraints = []
            self.cmaps = []
            self.vsites2 = []
            self.vsites3 = []
            self.has_virtual_sites = False
            self.has_nbfix_terms = False

        def findExclusionsFromBonds(self, genpairs):
            """Find exclusions between atoms separated by up to nrexcl bonds if genpairs is false,
               or up to 2 bonds if genpairs is true.
            """
            bondedTo = [set() for i in range(len(self.atoms))]
            for fields in self.bonds:
                i = int(fields[0])-1
                j = int(fields[1])-1
                bondedTo[i].add(j)
                bondedTo[j].add(i)

            # Identify all neighbors of each atom with each separation.

            bondedWithSeparation = [bondedTo]
            maxBonds = self.nrexcl
            if genpairs:
                maxBonds = min(maxBonds, 2)
            for i in range(maxBonds-1):
                lastBonds = bondedWithSeparation[-1]
                newBonds = deepcopy(lastBonds)
                for atom in range(len(self.atoms)):
                    for a1 in lastBonds[atom]:
                        for a2 in bondedTo[a1]:
                            newBonds[atom].add(a2)
                bondedWithSeparation.append(newBonds)

            # Build the list of pairs.

            pairs = []
            for atom in range(len(self.atoms)):
                for otherAtom in bondedWithSeparation[-1][atom]:
                    if otherAtom > atom:
                        pairs.append((atom, otherAtom))
            return pairs

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
            if len(self._ifStack) != len(self._elseStack):
                raise RuntimeError('#if/#else stack out of sync')

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
                valueStart = stripped.find(name, len(command))+len(name)+1
                value = line[valueStart:].strip()
                value = value or '1' # Default define is 1
                self._defines[name] = value
            elif command == '#ifdef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                name = fields[1]
                self._ifStack.append(name in self._defines)
                self._elseStack.append(False)
            elif command == '#undef':
                # Un-define a variable
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                if fields[1] in self._defines:
                    self._defines.pop(fields[1])
            elif command == '#ifndef':
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError('Illegal line in .top file: '+line)
                name = fields[1]
                self._ifStack.append(name not in self._defines)
                self._elseStack.append(False)
            elif command == '#endif':
                # Pop an entry off the if stack.
                if len(self._ifStack) == 0:
                    raise ValueError('Unexpected line in .top file: '+line)
                del(self._ifStack[-1])
                del(self._elseStack[-1])
            elif command == '#else':
                # Reverse the last entry on the if stack
                if len(self._ifStack) == 0:
                    raise ValueError('Unexpected line in .top file: '+line)
                if self._elseStack[-1]:
                    raise ValueError('Unexpected line in .top file: '
                                     '#else has already been used ' + line)
                self._ifStack[-1] = (not self._ifStack[-1])
                self._elseStack[-1] = True

        elif not ignore:
            # Gromacs occasionally uses #define's to introduce specific
            # parameters for individual terms (for instance, this is how
            # ff99SB-ILDN is implemented). So make sure we do the appropriate
            # pre-processor replacements necessary
            line = _replace_defines(line, self._defines)
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
            elif self._currentCategory == 'constraints':
                self._processConstraint(line)
            elif self._currentCategory == 'settles':
                self._processSettles(line)
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
            elif self._currentCategory == 'pairtypes':
                self._processPairType(line)
            elif self._currentCategory == 'cmaptypes':
                self._processCmapType(line)
            elif self._currentCategory == 'nonbond_params':
                self._processNonbondType(line)
            elif self._currentCategory == 'virtual_sites2' or self._currentCategory == 'dummies2':
                self._processVirtualSites2(line)
            elif self._currentCategory == 'virtual_sites3' or self._currentCategory == 'dummies3':
                self._processVirtualSites3(line)
            elif self._currentCategory.startswith('virtual_sites') or self._currentCategory.startswith('dummies'):
                if self._currentMoleculeType is None:
                    raise ValueError('Found %s before [ moleculetype ]' %
                                     self._currentCategory)
                self._currentMoleculeType.has_virtual_sites = True

    def _processDefaults(self, line):
        """Process the [ defaults ] line."""
        fields = line.split()
        if len(fields) < 5:
            # fudgeLJ and fudgeQQ not specified, assumed 1.0 by default
            if len(fields) == 3:
                fields.append(1.0)
                fields.append(1.0)
            else:
                raise ValueError('Too few fields in [ defaults ] line: '+line)
        if fields[0] != '1':
            raise ValueError('Unsupported nonbonded type: '+fields[0])
        if not fields[1] in ('1', '2', '3'):
            raise ValueError('Unsupported combination rule: '+fields[1])
        if fields[2].lower() == 'no':
            self._genpairs = False
        self._defaults = fields

    def _processMoleculeType(self, line):
        """Process a line in the [ moleculetypes ] category."""
        fields = line.split()
        if len(fields) < 1:
            raise ValueError('Too few fields in [ moleculetypes ] line: '+line)
        type = GromacsTopFile._MoleculeType(fields[0], int(fields[1]))
        self._moleculeTypes[fields[0]] = type
        self._currentMoleculeType = type

    def _processMolecule(self, line):
        """Process a line in the [ molecules ] category."""
        fields = line.split()
        if len(fields) < 2:
            raise ValueError('Too few fields in [ molecules ] line: '+line)
        self._molecules.append((fields[0], int(fields[1])))

    def _processAtom(self, line):
        """Process a line in the [ atoms ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ atoms ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ atoms ] line: '+line)
        self._currentMoleculeType.atoms.append(fields)

    def _processBond(self, line):
        """Process a line in the [ bonds ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ bonds ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 3:
            raise ValueError('Too few fields in [ bonds ] line: '+line)
        if fields[2] not in ('1', '2'):
                raise ValueError('Unsupported function type in [ bonds ] line: '+line)
        self._currentMoleculeType.bonds.append(fields)

    def _processAngle(self, line):
        """Process a line in the [ angles ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ angles ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 4:
            raise ValueError('Too few fields in [ angles ] line: '+line)
        if fields[3] not in ('1', '2', '5'):
            raise ValueError('Unsupported function type in [ angles ] line: '+line)
        self._currentMoleculeType.angles.append(fields)

    def _processDihedral(self, line):
        """Process a line in the [ dihedrals ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ dihedrals ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ dihedrals ] line: '+line)
        if fields[4] not in ('1', '2', '3', '4', '5', '9'):
            raise ValueError('Unsupported function type in [ dihedrals ] line: '+line)
        self._currentMoleculeType.dihedrals.append(fields)

    def _processExclusion(self, line):
        """Process a line in the [ exclusions ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ exclusions ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 2:
            raise ValueError('Too few fields in [ exclusions ] line: '+line)
        self._currentMoleculeType.exclusions.append(fields)

    def _processPair(self, line):
        """Process a line in the [ pairs ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ pairs ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 3:
            raise ValueError('Too few fields in [ pairs ] line: '+line)
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ pairs ] line: '+line)
        self._currentMoleculeType.pairs.append(fields)

    def _processConstraint(self, line):
        """Process a line in the [ constraints ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ constraints ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 4:
            raise ValueError('Too few fields in [ constraints ] line: '+line)
        self._currentMoleculeType.constraints.append(fields)

    def _processSettles(self, line):
        """Process a line in the [ settles ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ settles ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 4:
            raise ValueError('Too few fields in [ settles ] line: '+line)
        atom = int(fields[0])
        self._currentMoleculeType.constraints.append([str(atom), str(atom+1), fields[1], fields[2]])
        self._currentMoleculeType.constraints.append([str(atom), str(atom+2), fields[1], fields[2]])
        self._currentMoleculeType.constraints.append([str(atom+1), str(atom+2), fields[1], fields[3]])

    def _processCmap(self, line):
        """Process a line in the [ cmaps ] category."""
        if self._currentMoleculeType is None:
            raise ValueError('Found [ cmap ] section before [ moleculetype ]')
        fields = line.split()
        if len(fields) < 6:
            raise ValueError('Too few fields in [ cmap ] line: '+line)
        self._currentMoleculeType.cmaps.append(fields)

    def _processAtomType(self, line):
        """Process a line in the [ atomtypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError('Too few fields in [ atomtypes ] line: '+line)
        if len(fields[3]) == 1:
            # Bonded type and atomic number are both missing.
            fields.insert(1, None)
            fields.insert(1, None)
        elif len(fields[4]) == 1 and fields[4].isalpha():
            if fields[1][0].isalpha():
                # Atomic number is missing.
                fields.insert(2, None)
            else:
                # Bonded type is missing.
                fields.insert(1, None)
        self._atomTypes[fields[0]] = fields

    def _processBondType(self, line):
        """Process a line in the [ bondtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ bondtypes ] line: '+line)
        if fields[2] not in ('1', '2'):
            raise ValueError('Unsupported function type in [ bondtypes ] line: '+line)
        self._bondTypes[tuple(fields[:3])] = fields

    def _processAngleType(self, line):
        """Process a line in the [ angletypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            raise ValueError('Too few fields in [ angletypes ] line: '+line)
        if fields[3] not in ('1', '2', '5'):
            raise ValueError('Unsupported function type in [ angletypes ] line: '+line)
        self._angleTypes[tuple(fields[:3])] = fields

    def _processDihedralType(self, line):
        """Process a line in the [ dihedraltypes ] category."""
        fields = line.split()
        if len(fields[2]) == 1 and fields[2].isdigit():
            # The third field contains the function type, meaning only two atom types are specified.
            # Interpret them as the two inner ones.
            fields = ['X', fields[0], fields[1], 'X']+fields[2:]
        if len(fields) < 7:
            raise ValueError('Too few fields in [ dihedraltypes ] line: '+line)
        if fields[4] not in ('1', '2', '3', '4', '5', '9'):
            raise ValueError('Unsupported function type in [ dihedraltypes ] line: '+line)
        key = tuple(fields[:5])
        if fields[4] == '9' and key in self._dihedralTypes:
            # There are multiple dihedrals defined for these atom types.
            self._dihedralTypes[key].append(fields)
        else:
            self._dihedralTypes[key] = [fields]

    def _processPairType(self, line):
        """Process a line in the [ pairtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ pairtypes] line: '+line)
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ pairtypes ] line: '+line)
        self._pairTypes[tuple(fields[:2])] = fields

    def _processCmapType(self, line):
        """Process a line in the [ cmaptypes ] category."""
        fields = line.split()
        if len(fields) < 8 or len(fields) < 8+int(fields[6])*int(fields[7]):
            raise ValueError('Too few fields in [ cmaptypes ] line: '+line)
        if fields[5] != '1':
            raise ValueError('Unsupported function type in [ cmaptypes ] line: '+line)
        self._cmapTypes[tuple(fields[:5])] = fields

    def _processNonbondType(self, line):
        """Process a line in the [ nonbond_params ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ nonbond_params ] line: '+line)
        if fields[2] != '1':
            raise ValueError('Unsupported function type in [ nonbond_params ] line: '+line)
        self._nonbondTypes[tuple(sorted(fields[:2]))] = fields

    def _processVirtualSites2(self, line):
        """Process a line in the [ virtual_sites2 ] category."""
        fields = line.split()
        if len(fields) < 5:
            raise ValueError('Too few fields in [ virtual_sites2 ] line: ' + line)
        if fields[3] != '1':
            raise ValueError('Unsupported function type in [ virtual_sites2 ] line: '+line)
        self._currentMoleculeType.vsites2.append(fields[:5])

    def _processVirtualSites3(self, line):
        """Process a line in the [ virtual_sites3 ] category."""
        fields = line.split()
        if len(fields) < 7:
            raise ValueError('Too few fields in [ virtual_sites3 ] line: ' + line)
        if fields[4] not in ('1', '4'):
            raise ValueError('Unsupported function type in [ virtual_sites3 ] line: '+line)
        self._currentMoleculeType.vsites3.append(fields)

    def __init__(self, file, periodicBoxVectors=None, unitCellDimensions=None, includeDir=None, defines=None):
        """Load a top file.

        Parameters
        ----------
        file : str
            the name of the file to load
        periodicBoxVectors : tuple of Vec3=None
            the vectors defining the periodic box
        unitCellDimensions : Vec3=None
            the dimensions of the crystallographic unit cell.  For
            non-rectangular unit cells, specify periodicBoxVectors instead.
        includeDir : string=None
            A directory in which to look for other files included from the
            top file. If not specified, we will attempt to locate a gromacs
            installation on your system. When gromacs is installed in
            /usr/local, this will resolve to /usr/local/gromacs/share/gromacs/top
        defines : dict={}
            preprocessor definitions that should be predefined when parsing the file
         """
        if includeDir is None:
            includeDir = _defaultGromacsIncludeDir()
        self._includeDirs = (os.path.dirname(file), includeDir)
        # Most of the gromacs water itp files for different forcefields,
        # unless the preprocessor #define FLEXIBLE is given, don't define
        # bonds between the water hydrogen and oxygens, but only give the
        # constraint distances and exclusions.
        self._defines = OrderedDict()
        self._defines['FLEXIBLE'] = True
        self._genpairs = True
        if defines is not None:
            for define, value in defines.items():
                self._defines[define] = value

        # Parse the file.

        self._currentCategory = None
        self._ifStack = []
        self._elseStack = []
        self._moleculeTypes = {}
        self._molecules = []
        self._currentMoleculeType = None
        self._atomTypes = {}
        self._bondTypes= {}
        self._angleTypes = {}
        self._dihedralTypes = {}
        self._pairTypes = {}
        self._cmapTypes = {}
        self._nonbondTypes = {}
        self._processFile(file)

        # Create the Topology from it.

        top = Topology()
        ## The Topology read from the prmtop file
        self.topology = top
        if periodicBoxVectors is not None:
            if unitCellDimensions is not None:
                raise ValueError("specify either periodicBoxVectors or unitCellDimensions, but not both")
            top.setPeriodicBoxVectors(periodicBoxVectors)
        else:
            top.setUnitCellDimensions(unitCellDimensions)
        PDBFile._loadNameReplacementTables()
        for moleculeName, moleculeCount in self._molecules:
            if moleculeName not in self._moleculeTypes:
                raise ValueError("Unknown molecule type: "+moleculeName)
            moleculeType = self._moleculeTypes[moleculeName]
            if moleculeCount > 0 and moleculeType.has_virtual_sites:
                raise ValueError('Virtual sites not yet supported by Gromacs parsers')

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

                    # Try to determine the element.

                    atomicNumber = self._atomTypes[fields[1]][2]
                    if atomicNumber is None:
                        # Try to guess the element from the name.
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
                    elif atomicNumber == '0':
                        element = None
                    else:
                        element = elem.Element.getByAtomicNumber(int(atomicNumber))
                    atoms.append(top.addAtom(atomName, element, r))

                # Add bonds to the topology

                for fields in moleculeType.bonds:
                    top.addBond(atoms[int(fields[0])-1], atoms[int(fields[1])-1])

    def createSystem(self, nonbondedMethod=ff.NoCutoff, nonbondedCutoff=1.0*unit.nanometer, constraints=None,
                     rigidWater=True, ewaldErrorTolerance=0.0005, removeCMMotion=True, hydrogenMass=None, switchDistance=None):
        """Construct an OpenMM System representing the topology described by this
        top file.

        Parameters
        ----------
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, or LJPME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with
            constraints. Allowed values are None, HBonds, AllBonds, or HAngles.
            Regardless of this value, constraints that are explicitly specified
            in the top file will always be included.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        ewaldErrorTolerance : float=0.0005
            The error tolerance to use if nonbondedMethod is Ewald, PME or LJPME.
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep their
            total mass the same.  If rigidWater is used to make water molecules
            rigid, then water hydrogens are not altered.
        switchDistance : float=None
            The distance at which the potential energy switching function is turned on for
            Lennard-Jones interactions. If this is None, no switching function will be used.

        Returns
        -------
        System
             the newly created System
        """

        # Build a list of atom types for NBFIX

        atom_types = []
        for moleculeName, moleculeCount in self._molecules:
            moleculeType = self._moleculeTypes[moleculeName]
            for _ in range(moleculeCount):
                for atom in moleculeType.atoms:
                    atom_types.append(atom[1])
        has_nbfix_terms = any([pair in self._nonbondTypes for pair in combinations_with_replacement(sorted(set(atom_types)), 2)])

        # Create the System.

        sys = mm.System()
        boxVectors = self.topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(*boxVectors)
        elif nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME, ff.LJPME):
            raise ValueError('Illegal nonbonded method for a non-periodic system')
        nb = mm.NonbondedForce()
        sys.addForce(nb)
        lj = None
        if has_nbfix_terms:
            lj = mm.CustomNonbondedForce('(a/r6)^2-b/r6; r6=r^6; a=acoef(type1, type2); b=bcoef(type1, type2)')
            lj.addPerParticleParameter('type')
            sys.addForce(lj)
        elif self._defaults[1] in ('1', '3'):
            lj = mm.CustomNonbondedForce('A1*A2/r^12-C1*C2/r^6')
            lj.addPerParticleParameter('C')
            lj.addPerParticleParameter('A')
            sys.addForce(lj)
        bonds = {}
        angles = {}
        periodic = None
        rb = None
        harmonicTorsion = None
        cmap = None
        mapIndices = {}
        bondIndices = []
        topologyAtoms = list(self.topology.atoms())
        exclusions = []
        pairs = []
        fudgeQQ = float(self._defaults[4])
        fudgeLJ = float(self._defaults[3])

        # Build a lookup table to let us process dihedrals more quickly.

        dihedralTypeTable = {}
        for key in self._dihedralTypes:
            if key[1] != 'X' and key[2] != 'X':
                if (key[1], key[2]) not in dihedralTypeTable:
                    dihedralTypeTable[(key[1], key[2])] = []
                dihedralTypeTable[(key[1], key[2])].append(key)
                if (key[2], key[1]) not in dihedralTypeTable:
                    dihedralTypeTable[(key[2], key[1])] = []
                dihedralTypeTable[(key[2], key[1])].append(key)
        wildcardDihedralTypes = []
        for key in self._dihedralTypes:
            if key[1] == 'X' or key[2] == 'X':
                wildcardDihedralTypes.append(key)
                for types in dihedralTypeTable.values():
                    types.append(key)

        if has_nbfix_terms:
            # Build a lookup table and angle/dihedral indices list to
            # let us handle exclusion manually.
            angleIndices = []
            torsionIndices = []
            atom_partners = defaultdict(lambda : defaultdict(set))
            atom_charges = []

        # Loop over molecules and create the specified number of each type.

        for moleculeName, moleculeCount in self._molecules:
            moleculeType = self._moleculeTypes[moleculeName]
            exclusionsFromBonds = moleculeType.findExclusionsFromBonds(self._genpairs)
            for i in range(moleculeCount):

                # Record the types of all atoms.

                baseAtomIndex = sys.getNumParticles()
                atomTypes = [atom[1] for atom in moleculeType.atoms]
                try:
                    bondedTypes = [self._atomTypes[t][1] for t in atomTypes]
                except KeyError as e:
                    raise ValueError('Unknown atom type: ' + e.message)
                bondedTypes = [b if b is not None else a for a, b in zip(atomTypes, bondedTypes)]

                # Add atoms.

                for fields in moleculeType.atoms:
                    if len(fields) >= 8:
                        mass = float(fields[7])
                    else:
                        mass = float(self._atomTypes[fields[1]][3])
                    sys.addParticle(mass)

                # Add bonds.

                atomBonds = [{} for x in range(len(moleculeType.atoms))]
                for fields in moleculeType.bonds:
                    atoms = [int(x)-1 for x in fields[:2]]
                    types = tuple(bondedTypes[i] for i in atoms)
                    bondType = fields[2]
                    reversedTypes = types[::-1]+(bondType,)
                    types = types+(bondType,)
                    if len(fields) >= 5:
                        params = fields[3:5]
                    elif types in self._bondTypes:
                        params = self._bondTypes[types][3:5]
                    elif reversedTypes in self._bondTypes:
                        params = self._bondTypes[reversedTypes][3:5]
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
                    elif bondType == '1':
                        if bondType not in bonds:
                            bonds[bondType] = mm.HarmonicBondForce()
                            sys.addForce(bonds[bondType])
                        bonds[bondType].addBond(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], length, float(params[1]))
                    elif bondType == '2':
                        if bondType not in bonds:
                            bonds[bondType] = mm.CustomBondForce('0.25*k*(r^2-r0^2)^2')
                            bonds[bondType].addPerBondParameter('r0')
                            bonds[bondType].addPerBondParameter('k')
                            bonds[bondType].setName('GROMOSBondForce')
                            sys.addForce(bonds[bondType])
                        bonds[bondType].addBond(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], (length, float(params[1])))
                    else:
                        raise ValueError('Internal error: bondType has unexpected value: '+bondType)
                    # Record information that will be needed for constraining angles.
                    atomBonds[atoms[0]][atoms[1]] = length
                    atomBonds[atoms[1]][atoms[0]] = length

                # Add angles.

                degToRad = math.pi/180
                for fields in moleculeType.angles:
                    atoms = [int(x)-1 for x in fields[:3]]
                    types = tuple(bondedTypes[i] for i in atoms)
                    angleType = fields[3]
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
                        if angleType in ('1', '5'):
                            if angleType not in angles:
                                angles[angleType] = mm.HarmonicAngleForce()
                                sys.addForce(angles[angleType])
                            angles[angleType].addAngle(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], theta, float(params[1]))
                            if angleType == '5':
                                # This is a Urey-Bradley term, so also add the bond.
                                if '1' not in bonds:
                                    bonds['1'] = mm.HarmonicBondForce()
                                    sys.addForce(bonds['1'])
                                k = float(params[3])
                                if k != 0:
                                    bonds['1'].addBond(baseAtomIndex + atoms[0], baseAtomIndex + atoms[2], float(params[2]), k)
                        elif angleType == '2':
                            if angleType not in angles:
                                angles[angleType] = mm.CustomAngleForce('0.5*k*(cos(theta)-cos(theta0))^2')
                                angles[angleType].addPerAngleParameter('theta0')
                                angles[angleType].addPerAngleParameter('k')
                                angles[angleType].setName('GROMOSAngleForce')
                                sys.addForce(angles[angleType])
                            angles[angleType].addAngle(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], (theta, float(params[1])))
                        else:
                            raise ValueError('Internal error: angleType has unexpected value: '+angleType)

                # Add torsions.

                for fields in moleculeType.dihedrals:
                    atoms = [int(x)-1 for x in fields[:4]]
                    types = tuple(bondedTypes[i] for i in atoms)
                    dihedralType = fields[4]
                    reversedTypes = types[::-1]+(dihedralType,)
                    types = types+(dihedralType,)
                    if (dihedralType in ('1', '4', '5', '9') and len(fields) > 7) or (dihedralType == '3' and len(fields) > 10) or (dihedralType == '2' and len(fields) > 6):
                        paramsList = [fields]
                    else:
                        # Look for a matching dihedral type.
                        paramsList = None
                        if (types[1], types[2]) in dihedralTypeTable:
                            dihedralTypes = dihedralTypeTable[(types[1], types[2])]
                        else:
                            dihedralTypes = wildcardDihedralTypes
                        for key in dihedralTypes:
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
                                periodic.addTorsion(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], int(float(params[7])), float(params[5])*degToRad, k)
                        elif dihedralType == '2':
                            # Harmonic torsion
                            k = float(params[6])
                            phi0 = float(params[5])
                            if k != 0:
                                if harmonicTorsion is None:
                                    harmonicTorsion = mm.CustomTorsionForce('0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi); pi = %.15g' % math.pi)
                                    harmonicTorsion.addPerTorsionParameter('theta0')
                                    harmonicTorsion.addPerTorsionParameter('k')
                                    harmonicTorsion.setName('HarmonicTorsionForce')
                                    sys.addForce(harmonicTorsion)
                                # map phi0 into correct space
                                phi0 = phi0 - 360 if phi0 > 180 else phi0
                                harmonicTorsion.addTorsion(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], (phi0*degToRad, k))
                        else:
                            # RB Torsion
                            c = [float(x) for x in params[5:11]]
                            if any(x != 0 for x in c):
                                if rb is None:
                                    rb = mm.RBTorsionForce()
                                    sys.addForce(rb)
                                if dihedralType == '5':
                                    # Convert Fourier coefficients to RB coefficients.
                                    c = [c[1]+0.5*(c[0]+c[2]), 0.5*(-c[0]+3*c[2]), -c[1]+4*c[3], -2*c[2], -4*c[3], 0]
                                rb.addTorsion(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], c[0], c[1], c[2], c[3], c[4], c[5])

                # Add CMAP terms.

                for fields in moleculeType.cmaps:
                    atoms = [int(x)-1 for x in fields[:5]]
                    types = tuple(bondedTypes[i] for i in atoms)
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
                            map.append(float(params[8+mapSize*((j+mapSize//2)%mapSize)+((i+mapSize//2)%mapSize)]))
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
                        q = float(params[4])

                    if has_nbfix_terms:
                        nb.addParticle(q, 1.0, 0.0)
                        atom_charges.append(q)
                        lj.addParticle([0])
                    else:
                        if self._defaults[1] == '1':
                            nb.addParticle(q, 1.0, 0.0)
                            lj.addParticle([math.sqrt(float(params[6])), math.sqrt(float(params[7]))])
                        elif self._defaults[1] == '2':
                            nb.addParticle(q, float(params[6]), float(params[7]))
                        elif self._defaults[1] == '3':
                            nb.addParticle(q, 1.0, 0.0)
                            sigma = float(params[6])
                            epsilon = float(params[7])
                            lj.addParticle([math.sqrt(4*epsilon*sigma**6), math.sqrt(4*epsilon*sigma**12)])

                for fields in moleculeType.bonds:
                    atoms = [int(x)-1 for x in fields[:2]]
                    bondIndices.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1]))
                for fields in moleculeType.constraints:
                    if fields[2] == '1':
                        atoms = [int(x)-1 for x in fields[:2]]
                        bondIndices.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1]))
                if has_nbfix_terms:
                    for fields in moleculeType.bonds:
                        atoms = [int(x)-1 for x in fields[:2]]
                        atom_partners[baseAtomIndex+atoms[0]]['bond'].add(baseAtomIndex+atoms[1])
                        atom_partners[baseAtomIndex+atoms[1]]['bond'].add(baseAtomIndex+atoms[0])
                    for fields in moleculeType.angles:
                        atoms = [int(x)-1 for x in fields[:3]]
                        angleIndices.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2]))
                        for pair in combinations(atoms, 2):
                            atom_partners[baseAtomIndex+pair[0]]['angle'].add(baseAtomIndex+pair[1])
                            atom_partners[baseAtomIndex+pair[1]]['angle'].add(baseAtomIndex+pair[0])
                    for fields in moleculeType.dihedrals:
                        atoms = [int(x)-1 for x in fields[:4]]
                        torsionIndices.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3]))
                        for pair in combinations(atoms, 2):
                            atom_partners[baseAtomIndex+pair[0]]['torsion'].add(baseAtomIndex+pair[1])
                            atom_partners[baseAtomIndex+pair[1]]['torsion'].add(baseAtomIndex+pair[0])

                # Record nonbonded exceptions.

                for fields in moleculeType.pairs:
                    atoms = [int(x)-1 for x in fields[:2]]
                    types = tuple(atomTypes[i] for i in atoms)
                    atom1params = nb.getParticleParameters(baseAtomIndex+atoms[0])
                    atom2params = nb.getParticleParameters(baseAtomIndex+atoms[1])
                    atom1params = [x.value_in_unit_system(unit.md_unit_system) for x in atom1params]
                    atom2params = [x.value_in_unit_system(unit.md_unit_system) for x in atom2params]
                    if len(fields) >= 5:
                        params = [float(x) for x in fields[3:5]]
                    elif types in self._pairTypes:
                        params = [float(x) for x in self._pairTypes[types][3:5]]
                    elif types[::-1] in self._pairTypes:
                        params = [float(x) for x in self._pairTypes[types[::-1]][3:5]]
                    elif not self._genpairs:
                        raise ValueError('No pair parameters defined for atom '
                                         'types %s and gen-pairs is "no"' % types)
                    elif has_nbfix_terms:
                        continue
                    else:
                        # Generate the parameters based on the atom parameters.
                        if self._defaults[1] == '2':
                            params = [0.5*(atom1params[1]+atom2params[1]), fudgeLJ*math.sqrt(atom1params[2]*atom2params[2])]
                        else:
                            atom1lj = lj.getParticleParameters(baseAtomIndex+atoms[0])
                            atom2lj = lj.getParticleParameters(baseAtomIndex+atoms[1])
                            params = [fudgeLJ*atom1lj[0]*atom2lj[0], fudgeLJ*atom1lj[1]*atom2lj[1]]
                    pairs.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], atom1params[0]*atom2params[0]*fudgeQQ, params[0], params[1]))
                for fields in moleculeType.exclusions:
                    atoms = [int(x)-1 for x in fields]
                    for atom in atoms[1:]:
                        exclusions.append((baseAtomIndex+atoms[0], baseAtomIndex+atom))
                for atoms in exclusionsFromBonds:
                    exclusions.append((baseAtomIndex+atoms[0], baseAtomIndex+atoms[1]))

                # Record virtual sites

                for fields in moleculeType.vsites2:
                    atoms = [int(x)-1 for x in fields[:3]]
                    c1 = float(fields[4])
                    vsite = mm.TwoParticleAverageSite(baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], (1-c1), c1)
                    sys.setVirtualSite(baseAtomIndex+atoms[0], vsite)
                for fields in moleculeType.vsites3:
                    atoms = [int(x)-1 for x in fields[:4]]
                    vsiteType = fields[4]
                    c1 = float(fields[5])
                    c2 = float(fields[6])
                    if vsiteType == '1':
                        vsite = mm.ThreeParticleAverageSite(baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], 1-c1-c2, c1, c2)
                    elif vsiteType == '4':
                        c3 = float(fields[7])
                        vsite = mm.OutOfPlaneSite(baseAtomIndex+atoms[1], baseAtomIndex+atoms[2], baseAtomIndex+atoms[3], c1, c2, c3)
                    else:
                        raise ValueError('Internal error: vsites3 has unexpected type: '+vsiteType)
                    sys.setVirtualSite(baseAtomIndex+atoms[0], vsite)

                # Add explicitly specified constraints.

                for fields in moleculeType.constraints:
                    atoms = [int(x)-1 for x in fields[:2]]
                    length = float(fields[3])
                    sys.addConstraint(baseAtomIndex+atoms[0], baseAtomIndex+atoms[1], length)

        # Create nonbonded exceptions.

        if not has_nbfix_terms:
            nb.createExceptionsFromBonds(bondIndices, fudgeQQ, fudgeLJ)
        else:
            excluded_atom_pairs = set() # save these pairs so we don't zero them out
            for tor in torsionIndices:
                # First check to see if atoms 1 and 4 are already excluded because
                # they are 1-2 or 1-3 pairs (would happen in 6-member rings or
                # fewer). Then check that they're not already added as exclusions
                if 'bond' in atom_partners[tor[3]] and tor[0] in atom_partners[tor[3]]['bond']: continue
                if 'angle' in atom_partners[tor[3]] and tor[0] in atom_partners[tor[3]]['angle']: continue
                key = min((tor[0], tor[3]),
                          (tor[3], tor[0]))
                if key in excluded_atom_pairs: continue # multiterm...
                params1 = self._atomTypes[atom_types[tor[0]]]
                params4 = self._atomTypes[atom_types[tor[3]]]
                q1 = atom_charges[tor[0]]
                rmin1 = float(params1[6])
                eps1 = float(params1[7])
                q4 = atom_charges[tor[3]]
                rmin4 = float(params4[6])
                eps4 = float(params4[7])

                charge_prod = fudgeQQ*q1*q4
                epsilon = math.sqrt(abs(eps1 * eps4))
                if self._defaults[1] == '2':
                   rmin14 = (rmin1 + rmin4) / 2
                else:
                   rmin14 = math.sqrt(rmin1 * rmin4)
                nb.addException(tor[0], tor[3], charge_prod, rmin14, epsilon)
                excluded_atom_pairs.add(key)

            # Add excluded atoms
            for atom_idx, atom in atom_partners.items():
                # Exclude all bonds and angles
                for atom2 in atom['bond']:
                    if atom2 > atom_idx:
                        nb.addException(atom_idx, atom2, 0.0, 1.0, 0.0)
                        excluded_atom_pairs.add((atom_idx, atom2))
                for atom2 in atom['angle']:
                    if ((atom_idx, atom2) in excluded_atom_pairs):
                        continue
                    if atom2 > atom_idx:
                        nb.addException(atom_idx, atom2, 0.0, 1.0, 0.0)
                        excluded_atom_pairs.add((atom_idx, atom2))
                for atom2 in atom['dihedral']:
                    if atom2 <= atom_idx: continue
                    if ((atom_idx, atom2) in excluded_atom_pairs):
                        continue
                    nb.addException(atom_idx, atom2, 0.0, 1.0, 0.0)

        for exclusion in exclusions:
            nb.addException(exclusion[0], exclusion[1], 0.0, 1.0, 0.0, True)

        if lj is not None:
            # We're using a CustomNonbondedForce for LJ interactions, so also create a CustomBondForce
            # to handle the exceptions.

            pair_bond = mm.CustomBondForce('-C/r^6+A/r^12')
            pair_bond.addPerBondParameter('C')
            pair_bond.addPerBondParameter('A')
            pair_bond.setName('LennardJonesExceptions')
            sys.addForce(pair_bond)
            for pair in pairs:
                nb.addException(pair[0], pair[1], pair[2], 1.0, 0.0, True)
                pair_bond.addBond(pair[0], pair[1], [pair[3], pair[4]])
            for i in range(nb.getNumExceptions()):
                ii, jj, q, eps, sig = nb.getExceptionParameters(i)
                lj.addExclusion(ii, jj)
        elif self._defaults[1] == '2':
            for pair in pairs:
                nb.addException(pair[0], pair[1], pair[2], pair[3], pair[4], True)

        # Finish configuring the NonbondedForce.

        methodMap = {ff.NoCutoff:mm.NonbondedForce.NoCutoff,
                     ff.CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     ff.CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic,
                     ff.Ewald:mm.NonbondedForce.Ewald,
                     ff.PME:mm.NonbondedForce.PME,
                     ff.LJPME:mm.NonbondedForce.LJPME}
        nb.setNonbondedMethod(methodMap[nonbondedMethod])
        nb.setCutoffDistance(nonbondedCutoff)
        nb.setEwaldErrorTolerance(ewaldErrorTolerance)
        if switchDistance is not None:
            nb.setUseSwitchingFunction(True)
            nb.setSwitchingDistance(switchDistance)
        if lj is not None:
            methodMap = {ff.NoCutoff:mm.CustomNonbondedForce.NoCutoff,
                         ff.CutoffNonPeriodic:mm.CustomNonbondedForce.CutoffNonPeriodic,
                         ff.CutoffPeriodic:mm.CustomNonbondedForce.CutoffPeriodic,
                         ff.Ewald:mm.CustomNonbondedForce.CutoffPeriodic,
                         ff.PME:mm.CustomNonbondedForce.CutoffPeriodic,
                         ff.LJPME:mm.CustomNonbondedForce.CutoffPeriodic}
            lj.setNonbondedMethod(methodMap[nonbondedMethod])
            lj.setCutoffDistance(nonbondedCutoff)
            if nonbondedMethod in (ff.PME, ff.LJPME, ff.Ewald, ff.CutoffPeriodic):
                lj.setUseLongRangeCorrection(True)
            if switchDistance is not None:
                lj.setUseSwitchingFunction(True)
                lj.setSwitchingDistance(switchDistance)
            lj.setName('LennardJonesForce')

        if has_nbfix_terms:
            atom_nbfix_types = set([])
            for pair in self._nonbondTypes:
                atom_nbfix_types.add(pair[0])
                atom_nbfix_types.add(pair[1])

            lj_idx_list = [0 for _ in atom_types]
            lj_radii, lj_depths = [], []
            atom_params = []
            num_lj_types = 0
            lj_type_list = []
            for i,atom_type in enumerate(atom_types):
                atom = self._atomTypes[atom_type]
                if lj_idx_list[i]: continue # already assigned
                ljtype = (float(atom[6]), float(atom[7]))
                atom_params.append(ljtype)
                num_lj_types += 1
                lj_idx_list[i] = num_lj_types
                lj_type_list.append(atom)
                for j in range(i+1, len(atom_types)):
                    atom_type2 = atom_types[j]
                    if lj_idx_list[j] > 0: continue # already assigned
                    atom2 = self._atomTypes[atom_type2]
                    ljtype2 = (float(atom2[6]), float(atom2[7]))
                    if atom2 is atom:
                        lj_idx_list[j] = num_lj_types
                    elif atom_type not in atom_nbfix_types:
                        # Only non-NBFIXed atom types can be compressed
                        if ljtype == ljtype2:
                            lj_idx_list[j] = num_lj_types

            # Now everything is assigned. Create the A-coefficient and
            # B-coefficient arrays
            acoef = [0 for i in range(num_lj_types*num_lj_types)]
            bcoef = acoef[:]
            for i in range(num_lj_types):
                namei = lj_type_list[i][0]
                for j in range(num_lj_types):
                    namej = lj_type_list[j][0]
                    try:
                        types = self._nonbondTypes[tuple(sorted((namei, namej)))]
                        params = (float(types[3]), float(types[4]))
                        if self._defaults[1] == '2':
                            c6 = 4 * params[1] * params[0]**6
                            c12 = 4 * params[1] * params[0]**12
                        else:
                            c6 = params[0]
                            c12 = params[1]
                    except KeyError:
                        params1 = atom_params[i]
                        params2 = atom_params[j]
                        if self._defaults[1] == '1':
                            c6 = math.sqrt(params1[0]*params2[0])
                            c12 = math.sqrt(params1[1]*params2[1])
                        else:
                            if self._defaults[1] == '2':
                                sigma = (params1[0] + params2[0]) / 2
                            else:
                                sigma = math.sqrt(params1[0] + params2[0])
                            epsilon = math.sqrt(params1[1] * params2[1])
                            c6 = 4 * epsilon * sigma**6
                            c12 = 4 * epsilon * sigma**12
                    acoef[i+num_lj_types*j] = math.sqrt(c12)
                    bcoef[i+num_lj_types*j] = c6
            lj.addTabulatedFunction('acoef', mm.Discrete2DFunction(num_lj_types, num_lj_types, acoef))
            lj.addTabulatedFunction('bcoef', mm.Discrete2DFunction(num_lj_types, num_lj_types, bcoef))
            for i, idx in enumerate(lj_idx_list):
                lj.setParticleParameters(i, [idx-1]) # adjust for indexing from 0

        # Adjust masses.

        if hydrogenMass is not None:
            for atom1, atom2 in self.topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if rigidWater and atom2.residue.name == 'HOH':
                    continue
                if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Add a CMMotionRemover.

        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())
        return sys

def _defaultGromacsIncludeDir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) for the gromacs binary
    'pdb2gmx' or 'gmx' in the PATH, or (3) just using the default gromacs
    install location, /usr/local/gromacs/share/gromacs/top """
    if 'GMXDATA' in os.environ:
        return os.path.join(os.environ['GMXDATA'], 'top')
    if 'GMXBIN' in os.environ:
        return os.path.abspath(os.path.join(os.environ['GMXBIN'], '..', 'share', 'gromacs', 'top'))

    pdb2gmx_path = shutil.which('pdb2gmx')
    if pdb2gmx_path is not None:
        return os.path.abspath(os.path.join(os.path.dirname(pdb2gmx_path), '..', 'share', 'gromacs', 'top'))
    else:
        gmx_path = shutil.which('gmx')
        if gmx_path is not None:
            return os.path.abspath(os.path.join(os.path.dirname(gmx_path), '..', 'share', 'gromacs', 'top'))

    return '/usr/local/gromacs/share/gromacs/top'
