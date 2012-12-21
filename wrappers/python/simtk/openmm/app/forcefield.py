"""
forcefield.py: Constructs OpenMM System objects based on a Topology and an XML force field description

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman, Mark Friedrichs
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

import os
import itertools
import xml.etree.ElementTree as etree
from math import sqrt, cos
import simtk.openmm as mm
import simtk.unit as unit
import element as elem
from simtk.openmm.app import Topology

# Enumerated values for nonbonded method

NoCutoff = object()
CutoffNonPeriodic = object()
CutoffPeriodic = object()
Ewald = object()
PME = object()

# Enumerated values for constraint type

HBonds = object()
AllBonds = object()
HAngles = object()

# A map of functions to parse elements of the XML file.

parsers = {}

class ForceField(object):
    """A ForceField constructs OpenMM System objects based on a Topology."""

    def __init__(self, *files):
        """Load one or more XML files and create a ForceField object based on them.
        
        Parameters:
         - files A list of XML files defining the force field.  Each entry may be an absolute file path, a path relative to the
           current working directory, or a path relative to this module's data subdirectory (for built in force fields).
        """
        self._atomTypes = {}
        self._templates = {}
        self._templateSignatures = {None:[]}
        self._atomClasses = {'':set()}
        self._forces = []
        self._scripts = []
        for file in files:
            try:
                tree = etree.parse(file)
            except IOError:
                tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', file))
            root = tree.getroot()
            
            # Load the atom types.
            
            if tree.getroot().find('AtomTypes') is not None:
                for type in tree.getroot().find('AtomTypes').findall('Type'):
                    element = None
                    if 'element' in type.attrib:
                        element = elem.get_by_symbol(type.attrib['element'])
                    self._atomTypes[type.attrib['name']] = (type.attrib['class'], float(type.attrib['mass']), element)
            
            # Load the residue templates.
            
            if tree.getroot().find('Residues') is not None:
                for residue in root.find('Residues').findall('Residue'):
                    resName = residue.attrib['name']
                    template = ForceField._TemplateData(resName)
                    self._templates[resName] = template
                    for atom in residue.findall('Atom'):
                        template.atoms.append(ForceField._TemplateAtomData(atom.attrib['name'], atom.attrib['type'], self._atomTypes[atom.attrib['type']][2]))
                    for site in residue.findall('VirtualSite'):
                        template.virtualSites.append(ForceField._VirtualSiteData(site))
                    for bond in residue.findall('Bond'):
                        b = (int(bond.attrib['from']), int(bond.attrib['to']))
                        template.bonds.append(b)
                        template.atoms[b[0]].bondedTo.append(b[1])
                        template.atoms[b[1]].bondedTo.append(b[0])
                    for bond in residue.findall('ExternalBond'):
                        b = int(bond.attrib['from'])
                        template.externalBonds.append(b)
                        template.atoms[b].externalBonds += 1
            for template in self._templates.values():
                template.signature = _createResidueSignature([atom.element for atom in template.atoms])
                if template.signature is None:
                    sigString = None
                else:
                    sigString = _signatureToString(template.signature)
                if sigString in self._templateSignatures:
                    self._templateSignatures[sigString].append(template)
                else:
                    self._templateSignatures[sigString] = [template]
            
            # Build sets of every atom type belonging to each class
            
            for type in self._atomTypes:
                atomClass = self._atomTypes[type][0]
                if atomClass in self._atomClasses:
                    typeSet = self._atomClasses[atomClass]
                else:
                    typeSet = set()
                    self._atomClasses[atomClass] = typeSet
                typeSet.add(type)
                self._atomClasses[''].add(type)
            
            # Load force definitions
            
            for child in root:
                if child.tag in parsers:
                    parsers[child.tag](child, self)
            
            # Load scripts
            
            for node in tree.getroot().findall('Script'):
                self._scripts.append(node.text)

    def _findAtomTypes(self, node, num):
        """Parse the attributes on an XML tag to find the set of atom types for each atom it involves."""
        types = []
        attrib = node.attrib
        for i in range(num):
            if num == 1:
                suffix = ''
            else:
                suffix = str(i+1)
            classAttrib = 'class'+suffix
            typeAttrib = 'type'+suffix
            if classAttrib in attrib:
                if typeAttrib in attrib:
                    raise ValueError('Tag specifies both a type and a class for the same atom: '+etree.tostring(node))
                if attrib[classAttrib] not in self._atomClasses:
                    return None # Unknown atom class
                types.append(self._atomClasses[attrib[classAttrib]])
            else:
                if typeAttrib not in attrib or attrib[typeAttrib] not in self._atomTypes:
                    return None # Unknown atom type
                types.append([attrib[typeAttrib]])
        return types

    def _parseTorsion(self, node):
        """Parse the node defining a torsion."""
        types = self._findAtomTypes(node, 4)
        if types is None:
            return None
        torsion = PeriodicTorsion(types)
        attrib = node.attrib
        index = 1
        while 'phase%d'%index in attrib:
            torsion.periodicity.append(int(attrib['periodicity%d'%index]))
            torsion.phase.append(float(attrib['phase%d'%index]))
            torsion.k.append(float(attrib['k%d'%index]))
            index += 1
        return torsion            
        
    class _SystemData:
        """Inner class used to encapsulate data about the system being created."""
        def __init__(self):
            self.atomType = {}
            self.atoms = []
            self.virtualSites = {}
            self.bonds = []
            self.angles = []
            self.propers = []
            self.impropers = []
            self.atomBonds = []
            self.isAngleConstrained = []

    class _TemplateData:
        """Inner class used to encapsulate data about a residue template definition."""
        def __init__(self, name):
            self.name = name
            self.atoms = []
            self.virtualSites = []
            self.bonds = []
            self.externalBonds = []

    class _TemplateAtomData:
        """Inner class used to encapsulate data about an atom in a residue template definition."""
        def __init__(self, name, type, element):
            self.name = name
            self.type = type
            self.element = element
            self.bondedTo = []
            self.externalBonds = 0

    class _BondData:
        """Inner class used to encapsulate data about a bond."""
        def __init__(self, atom1, atom2):
            self.atom1 = atom1
            self.atom2 = atom2
            self.isConstrained = False
            self.length = 0.0
    
    class _VirtualSiteData:
        """Inner class used to encapsulate data about a virtual site."""
        def __init__(self, node):
            attrib = node.attrib
            self.index = int(attrib['index'])
            self.type = attrib['type']
            if self.type == 'average2':
                self.atoms = [int(attrib['atom1']), int(attrib['atom2'])]
                self.weights = [float(attrib['weight1']), float(attrib['weight2'])]
            elif self.type == 'average3':
                self.atoms = [int(attrib['atom1']), int(attrib['atom2']), int(attrib['atom3'])]
                self.weights = [float(attrib['weight1']), float(attrib['weight2']), float(attrib['weight3'])]
            elif self.type == 'outOfPlane':
                self.atoms = [int(attrib['atom1']), int(attrib['atom2']), int(attrib['atom3'])]
                self.weights = [float(attrib['weight12']), float(attrib['weight13']), float(attrib['weightCross'])]
            else:
                raise ValueError('Unknown virtual site type: %s' % self.type)
            self.atoms = [x-self.index for x in self.atoms]

    def createSystem(self, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=True, removeCMMotion=True, **args):
        """Construct an OpenMM System representing a Topology with this force field.
        
        Parameters:
         - topology (Topology) The Topology for which to create a System
         - nonbondedMethod (object=NoCutoff) The method to use for nonbonded interactions.  Allowed values are
           NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
         - nonbondedCutoff (distance=1*nanometer) The cutoff distance to use for nonbonded interactions
         - constraints (object=None) Specifies which bonds and angles should be implemented with constraints.
           Allowed values are None, HBonds, AllBonds, or HAngles.
         - rigidWater (boolean=True) If true, water molecules will be fully rigid regardless of the value passed for the constraints argument
         - removeCMMotion (boolean=True) If true, a CMMotionRemover will be added to the System
         - args Arbitrary additional keyword arguments may also be specified.  This allows extra parameters to be specified that are specific to
           particular force fields.
        Returns: the newly created System
        """
        
        # Record atom indices
        
        data = ForceField._SystemData()
        atomIndices = {}
        for index, atom in enumerate(topology.atoms()):
            data.atoms.append(atom)
            atomIndices[atom] = index

        # Make a list of all bonds
        
        for bond in topology.bonds():
            if bond[0] in atomIndices and bond[1] in atomIndices:
                data.bonds.append(ForceField._BondData(atomIndices[bond[0]], atomIndices[bond[1]]))

        # Record which atoms are bonded to each other atom
        
        bondedToAtom = []
        for i in range(len(data.atoms)):
            bondedToAtom.append(set())
            data.atomBonds.append([])
        for i in range(len(data.bonds)):
            bond = data.bonds[i]
            bondedToAtom[bond.atom1].add(bond.atom2)
            bondedToAtom[bond.atom2].add(bond.atom1)
            data.atomBonds[bond.atom1].append(i)
            data.atomBonds[bond.atom2].append(i)

        # Find the template matching each residue and assign atom types.
        
        for chain in topology.chains():
            for res in chain.residues():
                template = None
                matches = None
                sig = _createResidueSignature([atom.element for atom in res.atoms()])
                if sig is not None:
                    signature = _signatureToString(sig)
                    if signature in self._templateSignatures:
                        for t in self._templateSignatures[signature]:
                            matches = _matchResidue(res, t, bondedToAtom, atomIndices)
                            if matches is not None:
                                template = t
                                break
                if matches is None:
                    # Check templates involving virtual sites
                    for t in self._templateSignatures[None]:
                        matches = _matchResidue(res, t, bondedToAtom, atomIndices)
                        if matches is not None:
                            template = t
                            break
                if matches is None:
                    raise ValueError('No template found for residue %d (%s)' % (res.index+1, res.name))
                for atom, match in zip(res.atoms(), matches):
                    data.atomType[atom] = template.atoms[match].type
                    for site in template.virtualSites:
                        if match == site.index:
                            data.virtualSites[atom] = site

        # Create the System and add atoms
        
        sys = mm.System()
        for atom in topology.atoms():
            sys.addParticle(self._atomTypes[data.atomType[atom]][1])
        
        # Set periodic boundary conditions.
        
        boxSize = topology.getUnitCellDimensions()
        if boxSize is not None:
            sys.setDefaultPeriodicBoxVectors((boxSize[0], 0, 0), (0, boxSize[1], 0), (0, 0, boxSize[2]))
        elif nonbondedMethod is not NoCutoff and nonbondedMethod is not CutoffNonPeriodic:
            raise ValueError('Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')

        # Make a list of all unique angles
        
        uniqueAngles = set()
        for bond in data.bonds:
            for atom in bondedToAtom[bond.atom1]:
                if atom != bond.atom2:
                    if atom < bond.atom2:
                        uniqueAngles.add((atom, bond.atom1, bond.atom2))
                    else:
                        uniqueAngles.add((bond.atom2, bond.atom1, atom))
            for atom in bondedToAtom[bond.atom2]:
                if atom != bond.atom1:
                    if atom > bond.atom1:
                        uniqueAngles.add((bond.atom1, bond.atom2, atom))
                    else:
                        uniqueAngles.add((atom, bond.atom2, bond.atom1))
        data.angles = sorted(list(uniqueAngles))
        
        # Make a list of all unique proper torsions
        
        uniquePropers = set()
        for angle in data.angles:
            for atom in bondedToAtom[angle[0]]:
                if atom != angle[1]:
                    if atom < angle[2]:
                        uniquePropers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniquePropers.add((angle[2], angle[1], angle[0], atom))
            for atom in bondedToAtom[angle[2]]:
                if atom != angle[1]:
                    if atom > angle[0]:
                        uniquePropers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniquePropers.add((atom, angle[2], angle[1], angle[0]))
        data.propers = sorted(list(uniquePropers))
        
        # Make a list of all unique improper torsions
        
        for atom in range(len(bondedToAtom)):
            bondedTo = bondedToAtom[atom]
            if len(bondedTo) > 2:
                for subset in itertools.combinations(bondedTo, 3):
                    data.impropers.append((atom, subset[0], subset[1], subset[2]))
        
        # Identify bonds that should be implemented with constraints
        
        if constraints == AllBonds or constraints == HAngles:
            for bond in data.bonds:
                bond.isConstrained = True
        elif constraints == HBonds:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                bond.isConstrained = atom1.name.startswith('H') or atom2.name.startswith('H')
        if rigidWater:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH':
                    bond.isConstrained = True
        
        # Identify angles that should be implemented with constraints
        
        if constraints == HAngles:
            for angle in data.angles:
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                numH = 0
                if atom1.name.startswith('H'):
                    numH += 1
                if atom3.name.startswith('H'):
                    numH += 1
                data.isAngleConstrained.append(numH == 2 or (numH == 1 and atom2.name.startswith('O')))
        else:
            data.isAngleConstrained = len(data.angles)*[False]
        if rigidWater:
            for i in range(len(data.angles)):
                angle = data.angles[i]
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                    data.isAngleConstrained[i] = True
        
        # Add virtual sites
        
        for atom in data.virtualSites:
            site = data.virtualSites[atom]
            index = atomIndices[atom]
            if site.type == 'average2':
                sys.setVirtualSite(index, mm.TwoParticleAverageSite(index+site.atoms[0], index+site.atoms[1], site.weights[0], site.weights[1]))
            elif site.type == 'average3':
                sys.setVirtualSite(index, mm.ThreeParticleAverageSite(index+site.atoms[0], index+site.atoms[1], index+site.atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'outOfPlane':
                sys.setVirtualSite(index, mm.OutOfPlaneSite(index+site.atoms[0], index+site.atoms[1], index+site.atoms[2], site.weights[0], site.weights[1], site.weights[2]))
        
        # Add forces to the System
        
        for force in self._forces:
            force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())
        
        # Execute scripts found in the XML files.
        
        for script in self._scripts:
            exec script
        return sys


def _createResidueSignature(elements):
    """Create a signature for a residue based on the elements of the atoms it contains."""
    counts = {}
    for element in elements:
        if element is None:
            return None # This residue contains "atoms" (probably virtual sites) that should match any element
        if element in counts:
            counts[element] += 1
        else:
            counts[element] = 1
    sig = []
    for c in counts:
        sig.append((c, counts[c]))
    sig.sort(key=lambda x: -x[0].mass)
    return sig


def _signatureToString(signature):
    """Convert the signature returned by _createResidueSignature() to a string."""
    s = ''
    for element, count in signature:
        s += element.symbol+str(count)
    return s


def _matchResidue(res, template, bondedToAtom, atomIndices):
    """Determine whether a residue matches a template and return a list of corresponding atoms.
    
    Parameters:
     - res (Residue) The residue to check
     - template (_TemplateData) The template to compare it to
     - bondedToAtom (list) Enumerates which other atoms each atom is bonded to
     - atomIndices (map) Maps from atoms to their indices in the System
    Returns: a list specifying which atom of the template each atom of the residue corresponds to,
    or None if it does not match the template
    """
    atoms = list(res.atoms())
    if len(atoms) != len(template.atoms):
        return None
    residueAtomBonds = []
    templateAtomBonds = []
    matches = len(atoms)*[0]
    hasMatch = len(atoms)*[False]
    
    # Translate from global to local atom indices, and record the bonds for each atom.
    
    renumberAtoms = {}
    for i in range(len(atoms)):
        renumberAtoms[atomIndices[atoms[i]]] = i
    bondedTo = []
    externalBonds = []
    for atom in atoms:
        bonds = [renumberAtoms[x] for x in bondedToAtom[atomIndices[atom]] if x in renumberAtoms]
        bondedTo.append(bonds)
        externalBonds.append(len([x for x in bondedToAtom[atomIndices[atom]] if x not in renumberAtoms]))
    if _findAtomMatches(atoms, template, bondedTo, externalBonds, matches, hasMatch, 0):
        return matches
    return None


def _findAtomMatches(atoms, template, bondedTo, externalBonds, matches, hasMatch, position):
    """This is called recursively from inside _matchResidue() to identify matching atoms."""
    if position == len(atoms):
        return True
    elem = atoms[position].element
    for i in range(len(atoms)):
        atom = template.atoms[i]
        if (atom.element == elem or atom.element is None) and not hasMatch[i] and len(atom.bondedTo) == len(bondedTo[position]) and atom.externalBonds == externalBonds[position]:
            # See if the bonds for this identification are consistent
            
            allBondsMatch = all((bonded > position or matches[bonded] in atom.bondedTo for bonded in bondedTo[position]))
            if allBondsMatch:
                # This is a possible match, so trying matching the rest of the residue.
                
                matches[position] = i
                hasMatch[i] = True
                if _findAtomMatches(atoms, template, bondedTo, externalBonds, matches, hasMatch, position+1):
                    return True
                hasMatch[i] = False
    return False


# The following classes are generators that know how to create Force subclasses and add them to a System that is being
# created.  Each generator class must define two methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding generator object; 2) a createForce() method that constructs the Force object and adds it
# to the System.  The static method should be added to the parsers map.

## @private
class HarmonicBondGenerator:
    """A HarmonicBondGenerator constructs a HarmonicBondForce."""
    
    def __init__(self):
        self.types1 = []
        self.types2 = []
        self.length = []
        self.k = []
    
    @staticmethod
    def parseElement(element, ff):
        generator = HarmonicBondGenerator()
        ff._forces.append(generator)
        for bond in element.findall('Bond'):
            types = ff._findAtomTypes(bond, 2)
            if types is not None:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.length.append(float(bond.attrib['length']))
                generator.k.append(float(bond.attrib['k']))
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.HarmonicBondForce]
        if len(existing) == 0:
            force = mm.HarmonicBondForce()
            sys.addForce(force)
        else:
            force = existing[0]
        for bond in data.bonds:
            type1 = data.atomType[data.atoms[bond.atom1]]
            type2 = data.atomType[data.atoms[bond.atom2]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
                    bond.length = self.length[i]
                    if bond.isConstrained:
                        sys.addConstraint(bond.atom1, bond.atom2, self.length[i])
                    elif self.k[i] != 0:
                        force.addBond(bond.atom1, bond.atom2, self.length[i], self.k[i])
                    break

parsers["HarmonicBondForce"] = HarmonicBondGenerator.parseElement


## @private
class HarmonicAngleGenerator:
    """A HarmonicAngleGenerator constructs a HarmonicAngleForce."""
    
    def __init__(self):
        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.angle = []
        self.k = []
    
    @staticmethod
    def parseElement(element, ff):
        generator = HarmonicAngleGenerator()
        ff._forces.append(generator)
        for angle in element.findall('Angle'):
            types = ff._findAtomTypes(angle, 3)
            if types is not None:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])
                generator.angle.append(float(angle.attrib['angle']))
                generator.k.append(float(angle.attrib['k']))
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.HarmonicAngleForce]
        if len(existing) == 0:
            force = mm.HarmonicAngleForce()
            sys.addForce(force)
        else:
            force = existing[0]
        for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):
            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                    if isConstrained:
                        # Find the two bonds that make this angle.
                        
                        bond1 = None
                        bond2 = None
                        for bond in data.atomBonds[angle[1]]:
                            atom1 = data.bonds[bond].atom1
                            atom2 = data.bonds[bond].atom2
                            if atom1 == angle[0] or atom2 == angle[0]:
                                bond1 = bond
                            elif atom1 == angle[2] or atom2 == angle[2]:
                                bond2 = bond
                        
                        # Compute the distance between atoms and add a constraint
                        
                        if bond1 is not None and bond2 is not None:
                            l1 = data.bonds[bond1].length
                            l2 = data.bonds[bond2].length
                            if l1 is not None and l2 is not None:
                                length = sqrt(l1*l1 + l2*l2 - 2*l1*l2*cos(self.angle[i]))
                                sys.addConstraint(angle[0], angle[2], length)
                    elif self.k[i] != 0:
                        force.addAngle(angle[0], angle[1], angle[2], self.angle[i], self.k[i])
                    break

parsers["HarmonicAngleForce"] = HarmonicAngleGenerator.parseElement


## @private
class PeriodicTorsion:
    """A PeriodicTorsion records the information for a periodic torsion definition."""

    def __init__(self, types):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.periodicity = []
        self.phase = []
        self.k = []

## @private
class PeriodicTorsionGenerator:
    """A PeriodicTorsionGenerator constructs a PeriodicTorsionForce."""
    
    def __init__(self):
        self.proper = []
        self.improper = []
    
    @staticmethod
    def parseElement(element, ff):
        generator = PeriodicTorsionGenerator()
        generator.ff = ff
        ff._forces.append(generator)
        for torsion in element.findall('Proper'):
            torsion = ff._parseTorsion(torsion)
            if torsion is not None:
                generator.proper.append(torsion)
        for torsion in element.findall('Improper'):
            torsion = ff._parseTorsion(torsion)
            if torsion is not None:
                generator.improper.append(torsion)
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.PeriodicTorsionForce]
        if len(existing) == 0:
            force = mm.PeriodicTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]
        wildcard = self.ff._atomClasses['']
        for torsion in data.propers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            match = None
            for tordef in self.proper:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                if (type2 in types2 and type3 in types3 and type4 in types4 and type1 in types1) or (type2 in types3 and type3 in types2 and type4 in types1 and type1 in types4):
                    hasWildcard = (wildcard in (types1, types2, types3, types4))
                    if match is None or not hasWildcard: # Prefer specific definitions over ones with wildcards
                        match = tordef
                    if not hasWildcard:
                        break
            if match is not None:
                for i in range(len(match.phase)):
                    if match.k[i] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], match.periodicity[i], match.phase[i], match.k[i])
        for torsion in data.impropers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            done = False
            for tordef in self.improper:
                if done:
                    break
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                if type1 in types1:
                    for (t2, t3, t4) in itertools.permutations(((type2, 1), (type3, 2), (type4, 3))):
                        if t2[0] in types2 and t3[0] in types3 and t4[0] in types4:
                            # Workaround to be more consistent with AMBER.  It uses wildcards to define most of its
                            # impropers, which leaves the ordering ambigous.  It then follows some bizarre rules
                            # to pick the order.
                            a1 = torsion[t2[1]]
                            a2 = torsion[t3[1]]
                            e1 = data.atoms[a1].element
                            e2 = data.atoms[a2].element
                            if e1 == e2 and a1 > a2:
                                (a1, a2) = (a2, a1)
                            elif e1 != elem.carbon and (e2 == elem.carbon or e1.mass < e2.mass):
                                (a1, a2) = (a2, a1)
                            for i in range(len(tordef.phase)):
                                if tordef.k[i] != 0:
                                    force.addTorsion(a1, a2, torsion[0], torsion[t4[1]], tordef.periodicity[i], tordef.phase[i], tordef.k[i])
                            done = True
                            break

parsers["PeriodicTorsionForce"] = PeriodicTorsionGenerator.parseElement


## @private
class RBTorsion:
    """An RBTorsion records the information for a Ryckaert-Bellemans torsion definition."""

    def __init__(self, types, c):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.c = c

## @private
class RBTorsionGenerator:
    """An RBTorsionGenerator constructs an RBTorsionForce."""
    
    def __init__(self):
        self.proper = []
        self.improper = []
    
    @staticmethod
    def parseElement(element, ff):
        generator = RBTorsionGenerator()
        generator.ff = ff
        ff._forces.append(generator)
        for torsion in element.findall('Proper'):
            types = ff._findAtomTypes(torsion, 4)
            if types is not None:
                generator.proper.append(RBTorsion(types, [float(torsion.attrib['c'+str(i)]) for i in range(6)]))
        for torsion in element.findall('Improper'):
            types = ff._findAtomTypes(torsion, 4)
            if types is not None:
                generator.improper.append(RBTorsion(types, [float(torsion.attrib['c'+str(i)]) for i in range(6)]))
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.RBTorsionForce]
        if len(existing) == 0:
            force = mm.RBTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]
        wildcard = self.ff._atomClasses['']
        for torsion in data.propers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            match = None
            for tordef in self.proper:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                if (type2 in types2 and type3 in types3 and type4 in types4 and type1 in types1) or (type2 in types3 and type3 in types2 and type4 in types1 and type1 in types4):
                    hasWildcard = (wildcard in (types1, types2, types3, types4))
                    if match is None or not hasWildcard: # Prefer specific definitions over ones with wildcards
                        match = tordef
                    if not hasWildcard:
                        break
            if match is not None:
                force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], match.c[0], match.c[1], match.c[2], match.c[3], match.c[4], match.c[5])
        for torsion in data.impropers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            done = False
            for tordef in self.improper:
                if done:
                    break
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                if type1 in types1:
                    for (t2, t3, t4) in itertools.permutations(((type2, 1), (type3, 2), (type4, 3))):
                        if t2[0] in types2 and t3[0] in types3 and t4[0] in types4:
                            # Workaround to be more consistent with AMBER.  It uses wildcards to define most of its
                            # impropers, which leaves the ordering ambigous.  It then follows some bizarre rules
                            # to pick the order.
                            a1 = torsion[t2[1]]
                            a2 = torsion[t3[1]]
                            e1 = data.atoms[a1].element
                            e2 = data.atoms[a2].element
                            if e1 == e2 and a1 > a2:
                                (a1, a2) = (a2, a1)
                            elif e1 != elem.carbon and (e2 == elem.carbon or e1.mass < e2.mass):
                                (a1, a2) = (a2, a1)
                            force.addTorsion(a1, a2, torsion[0], torsion[t4[1]], tordef.c[0], tordef.c[1], tordef.c[2], tordef.c[3], tordef.c[4], tordef.c[5])
                            done = True
                            break

parsers["RBTorsionForce"] = RBTorsionGenerator.parseElement


## @private
class CMAPTorsion:
    """A CMAPTorsion records the information for a CMAP torsion definition."""

    def __init__(self, types, map):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.types5 = types[4]
        self.map = map

## @private
class CMAPTorsionGenerator:
    """A CMAPTorsionGenerator constructs a CMAPTorsionForce."""
    
    def __init__(self):
        self.torsions = []
        self.maps = []
    
    @staticmethod
    def parseElement(element, ff):
        generator = CMAPTorsionGenerator()
        generator.ff = ff
        ff._forces.append(generator)
        for map in element.findall('Map'):
            values = [float(x) for x in map.text.split()]
            size = sqrt(len(values))
            if size*size != len(values):
                raise ValueError('CMAP must have the same number of elements along each dimension')
            generator.maps.append(values)
        for torsion in element.findall('Torsion'):
            types = ff._findAtomTypes(torsion, 5)
            if types is not None:
                generator.torsions.append(CMAPTorsion(types, int(torsion.attrib['map'])))
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.CMAPTorsionForce]
        if len(existing) == 0:
            force = mm.CMAPTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]
        for map in self.maps:
            force.addMap(int(sqrt(len(map))), map)
        
        # Find all chains of length 5
        
        uniqueTorsions = set()
        for torsion in data.propers:
            for bond in (data.bonds[x] for x in data.atomBonds[torsion[0]]):
                if bond.atom1 == torsion[0]:
                    atom = bond.atom2
                else:
                    atom = bond.atom1
                if atom != torsion[1]:
                    uniqueTorsions.add((atom, torsion[0], torsion[1], torsion[2], torsion[3]))
            for bond in (data.bonds[x] for x in data.atomBonds[torsion[3]]):
                if bond.atom1 == torsion[3]:
                    atom = bond.atom2
                else:
                    atom = bond.atom1
                if atom != torsion[2]:
                    uniqueTorsions.add((torsion[0], torsion[1], torsion[2], torsion[3], atom))
        torsions = sorted(list(uniqueTorsions))
        wildcard = self.ff._atomClasses['']
        for torsion in torsions:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            type5 = data.atomType[data.atoms[torsion[4]]]
            match = None
            for tordef in self.torsions:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                types5 = tordef.types5
                if (type1 in types1 and type2 in types2 and type3 in types3 and type4 in types4 and type5 in types5) or (type1 in types5 and type2 in types4 and type3 in types3 and type4 in types2 and type5 in types1):
                    hasWildcard = (wildcard in (types1, types2, types3, types4, types5))
                    if match is None or not hasWildcard: # Prefer specific definitions over ones with wildcards
                        match = tordef
                    if not hasWildcard:
                        break
            if match is not None:
                force.addTorsion(match.map, torsion[0], torsion[1], torsion[2], torsion[3], torsion[1], torsion[2], torsion[3], torsion[4])

parsers["CMAPTorsionForce"] = CMAPTorsionGenerator.parseElement


## @private
class NonbondedGenerator:
    """A NonbondedGenerator constructs a NonbondedForce."""
    
    def __init__(self, coulomb14scale, lj14scale):
        self.coulomb14scale = coulomb14scale
        self.lj14scale = lj14scale
        self.typeMap = {}

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, NonbondedGenerator)]
        if len(existing) == 0:
            generator = NonbondedGenerator(float(element.attrib['coulomb14scale']), float(element.attrib['lj14scale']))
            ff._forces.append(generator)
        else:
            # Multiple <NonbondedForce> tags were found, probably in different files.  Simply add more types to the existing one.
            generator = existing[0]
            if generator.coulomb14scale != float(element.attrib['coulomb14scale']) or generator.lj14scale != float(element.attrib['lj14scale']):
                raise ValueError('Found multiple NonbondedForce tags with different 1-4 scales') 
        for atom in element.findall('Atom'):
            types = ff._findAtomTypes(atom, 1)
            if types is not None:
                values = (float(atom.attrib['charge']), float(atom.attrib['sigma']), float(atom.attrib['epsilon']))
                for t in types[0]:
                    generator.typeMap[t] = values
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic,
                     Ewald:mm.NonbondedForce.Ewald,
                     PME:mm.NonbondedForce.PME}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for NonbondedForce') 
        force = mm.NonbondedForce()
        for atom in data.atoms:
            t = data.atomType[atom]
            if t in self.typeMap:
                values = self.typeMap[t]
                force.addParticle(values[0], values[1], values[2])
            else:
                raise ValueError('No nonbonded parameters defined for atom type '+t)
        # Create exceptions based on bonds and virtual sites.
        bondIndices = []
        for bond in data.bonds:
            bondIndices.append((bond.atom1, bond.atom2))
        for i in range(sys.getNumParticles()):
            if sys.isVirtualSite(i):
                site = sys.getVirtualSite(i)
                for j in range(site.getNumParticles()):
                    bondIndices.append((i, site.getParticle(j)))
        force.createExceptionsFromBonds(bondIndices, self.coulomb14scale, self.lj14scale)
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        if 'ewaldErrorTolerance' in args:
            force.setEwaldErrorTolerance(args['ewaldErrorTolerance'])
        sys.addForce(force)

parsers["NonbondedForce"] = NonbondedGenerator.parseElement


## @private
class GBSAOBCGenerator:
    """A GBSAOBCGenerator constructs a GBSAOBCForce."""
    
    def __init__(self):
        self.typeMap = {}

    @staticmethod
    def parseElement(element, ff):
        generator = GBSAOBCGenerator()
        ff._forces.append(generator)
        for atom in element.findall('Atom'):
            types = ff._findAtomTypes(atom, 1)
            if types is not None:
                values = (float(atom.attrib['charge']), float(atom.attrib['radius']), float(atom.attrib['scale']))
                for t in types[0]:
                    generator.typeMap[t] = values
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for GBSAOBCForce') 
        force = mm.GBSAOBCForce()
        for atom in data.atoms:
            t = data.atomType[atom]
            if t in self.typeMap:
                values = self.typeMap[t]
                force.addParticle(values[0], values[1], values[2])
            else:
                raise ValueError('No GBSAOBC parameters defined for atom type '+t) 
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        if 'soluteDielectric' in args:
            force.setSoluteDielectric(float(args['soluteDielectric']))
        if 'solventDielectric' in args:
            force.setSolventDielectric(float(args['solventDielectric']))
        sys.addForce(force)

parsers["GBSAOBCForce"] = GBSAOBCGenerator.parseElement


## @private
class GBVIGenerator:

    """A GBVIGenerator constructs a GBVIForce."""
    
    def __init__(self,ff):

        self.ff = ff
        self.fixedParameters = {}
        self.fixedParameters['soluteDielectric'] = 1.0
        self.fixedParameters['solventDielectric'] = 78.3
        self.fixedParameters['scalingMethod'] = 1
        self.fixedParameters['quinticUpperBornRadiusLimit'] = 5.0
        self.fixedParameters['quinticLowerLimitFactor'] = 0.8

        self.typeMap = {}

    @staticmethod
    def parseElement(element, ff):
        generator = GBVIGenerator(ff)
        for key in generator.fixedParameters.iterkeys():
            if (key in element.attrib):
                generator.fixedParameters[key] = float(element.attrib[key])

        ff._forces.append(generator)
        for atom in element.findall('Atom'):
            types = ff._findAtomTypes(atom, 1)
            if types is not None:
                values = (float(atom.attrib['charge']), float(atom.attrib['radius']), float(atom.attrib['gamma']))
                for t in types[0]:
                    generator.typeMap[t] = values
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        methodMap = {NoCutoff:mm.NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic}

        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for GB/VI Force') 

        # add particles

        force = mm.GBVIForce()
        for atom in data.atoms:
            t = data.atomType[atom]
            if t in self.typeMap:
                values = self.typeMap[t]
                force.addParticle(values[0], values[1], values[2])
            else:
                raise ValueError('No GB/VI parameters defined for atom type '+t) 

        # get HarmonicBond generator -- exit if not found

        hbGenerator = 0
        for generator in self.ff._forces:
            if (generator.__class__.__name__ == 'HarmonicBondGenerator'): 
               hbGenerator = generator
               break

        if (hbGenerator == 0):
            raise ValueError('HarmonicBondGenerator not found.') 

        # add bonds

        for bond in data.bonds:
            type1 = data.atomType[data.atoms[bond.atom1]]
            type2 = data.atomType[data.atoms[bond.atom2]]
            hit = 0
            for i in range(len(hbGenerator.types1)):
                types1 = hbGenerator.types1[i]
                types2 = hbGenerator.types2[i]
                if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
                    #bond.length = hbGenerator.length[i]
                    force.addBond(bond.atom1, bond.atom2, hbGenerator.length[i])

        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        force.setSolventDielectric(self.fixedParameters['solventDielectric'])
        force.setSoluteDielectric(self.fixedParameters['soluteDielectric'])
        force.setBornRadiusScalingMethod(self.fixedParameters['scalingMethod'])
        force.setQuinticLowerLimitFactor(self.fixedParameters['quinticLowerLimitFactor'])
        force.setQuinticUpperBornRadiusLimit(self.fixedParameters['quinticUpperBornRadiusLimit'])
        
        sys.addForce(force)

parsers["GBVIForce"] = GBVIGenerator.parseElement

## @private
class CustomBondGenerator:
    """A CustomBondGenerator constructs a CustomBondForce."""
    
    def __init__(self):
        self.types1 = []
        self.types2 = []
        self.globalParams = {}
        self.perBondParams = []
        self.paramValues = []
    
    @staticmethod
    def parseElement(element, ff):
        generator = CustomBondGenerator()
        ff._forces.append(generator)
        generator.energy = element.attrib['energy']
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerBondParameter'):
            generator.perBondParams.append(param.attrib['name'])
        for bond in element.findall('Bond'):
            types = ff._findAtomTypes(bond, 2)
            if types is not None:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.paramValues.append([float(bond.attrib[param]) for param in generator.perBondParams])
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = mm.CustomBondForce(self.energy)
        sys.addForce(force)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perBondParams:
            force.addPerBondParameter(param)
        for bond in data.bonds:
            type1 = data.atomType[data.atoms[bond.atom1]]
            type2 = data.atomType[data.atoms[bond.atom2]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
                    force.addBond(bond.atom1, bond.atom2, self.paramValues[i])
                    break

parsers["CustomBondForce"] = CustomBondGenerator.parseElement


## @private
class CustomAngleGenerator:
    """A CustomAngleGenerator constructs a CustomAngleForce."""
    
    def __init__(self):
        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.globalParams = {}
        self.perAngleParams = []
        self.paramValues = []
    
    @staticmethod
    def parseElement(element, ff):
        generator = CustomAngleGenerator()
        ff._forces.append(generator)
        generator.energy = element.attrib['energy']
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerAngleParameter'):
            generator.perAngleParams.append(param.attrib['name'])
        for angle in element.findall('Angle'):
            types = ff._findAtomTypes(angle, 3)
            if types is not None:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])
                generator.paramValues.append([float(angle.attrib[param]) for param in generator.perAngleParams])
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = mm.CustomAngleForce(self.energy)
        sys.addForce(force)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perAngleParams:
            force.addPerAngleParameter(param)
        for angle in data.angles:
            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                    force.addAngle(angle[0], angle[1], angle[2], self.paramValues[i])
                    break

parsers["CustomAngleForce"] = CustomAngleGenerator.parseElement


## @private
class CustomTorsion:
    """A CustomTorsion records the information for a custom torsion definition."""

    def __init__(self, types, paramValues):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.paramValues = paramValues

## @private
class CustomTorsionGenerator:
    """A CustomTorsionGenerator constructs a CustomTorsionForce."""
    
    def __init__(self):
        self.proper = []
        self.improper = []
        self.globalParams = {}
        self.perTorsionParams = []
    
    @staticmethod
    def parseElement(element, ff):
        generator = CustomTorsionGenerator()
        generator.ff = ff
        ff._forces.append(generator)
        generator.energy = element.attrib['energy']
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerTorsionParameter'):
            generator.perTorsionParams.append(param.attrib['name'])
        for torsion in element.findall('Proper'):
            types = ff._findAtomTypes(torsion, 4)
            if types is not None:
                generator.proper.append(CustomTorsion(types, [float(torsion.attrib[param]) for param in generator.perTorsionParams]))
        for torsion in element.findall('Improper'):
            types = ff._findAtomTypes(torsion, 4)
            if types is not None:
                generator.improper.append(CustomTorsion(types, [float(torsion.attrib[param]) for param in generator.perTorsionParams]))
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = mm.CustomTorsionForce(self.energy)
        sys.addForce(force)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perTorsionParams:
            force.addPerTorsionParameter(param)
        wildcard = self.ff._atomClasses['']
        for torsion in data.propers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            match = None
            for tordef in self.proper:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                if (type2 in types2 and type3 in types3 and type4 in types4 and type1 in types1) or (type2 in types3 and type3 in types2 and type4 in types1 and type1 in types4):
                    hasWildcard = (wildcard in (types1, types2, types3, types4))
                    if match is None or not hasWildcard: # Prefer specific definitions over ones with wildcards
                        match = tordef
                    if not hasWildcard:
                        break
            if match is not None:
                force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], match.paramValues)
        for torsion in data.impropers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            done = False
            for tordef in self.improper:
                if done:
                    break
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                if type1 in types1:
                    for (t2, t3, t4) in itertools.permutations(((type2, 1), (type3, 2), (type4, 3))):
                        if t2[0] in types2 and t3[0] in types3 and t4[0] in types4:
                            # Workaround to be more consistent with AMBER.  It uses wildcards to define most of its
                            # impropers, which leaves the ordering ambigous.  It then follows some bizarre rules
                            # to pick the order.
                            a1 = torsion[t2[1]]
                            a2 = torsion[t3[1]]
                            e1 = data.atoms[a1].element
                            e2 = data.atoms[a2].element
                            if e1 == e2 and a1 > a2:
                                (a1, a2) = (a2, a1)
                            elif e1 != elem.carbon and (e2 == elem.carbon or e1.mass < e2.mass):
                                (a1, a2) = (a2, a1)
                            force.addTorsion(a1, a2, torsion[0], torsion[t4[1]], tordef.paramValues)
                            done = True
                            break

parsers["CustomTorsionForce"] = CustomTorsionGenerator.parseElement


## @private
class CustomGBGenerator:
    """A CustomGBGenerator constructs a CustomGBForce."""
    
    def __init__(self):
        self.typeMap = {}
        self.globalParams = {}
        self.perParticleParams = []
        self.paramValues = []
        self.computedValues = []
        self.energyTerms = []
        self.functions = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomGBGenerator()
        ff._forces.append(generator)
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerParticleParameter'):
            generator.perParticleParams.append(param.attrib['name'])
        for atom in element.findall('Atom'):
            types = ff._findAtomTypes(atom, 1)
            if types is not None:
                values = [float(atom.attrib[param]) for param in generator.perParticleParams]
                for t in types[0]:
                    generator.typeMap[t] = values
        computationMap = {"SingleParticle" : mm.CustomGBForce.SingleParticle,
                          "ParticlePair" : mm.CustomGBForce.ParticlePair,
                          "ParticlePairNoExclusions" : mm.CustomGBForce.ParticlePairNoExclusions}
        for value in element.findall('ComputedValue'):
            generator.computedValues.append((value.attrib['name'], value.text, computationMap[value.attrib['type']]))
        for term in element.findall('EnergyTerm'):
            generator.energyTerms.append((term.text, computationMap[term.attrib['type']]))
        for function in element.findall("Function"):
            values = [float(x) for x in function.text.split()]
            generator.functions.append((function.attrib['name'], values, float(function.attrib['min']), float(function.attrib['max'])))
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.CustomGBForce.NoCutoff,
                     CutoffNonPeriodic:mm.CustomGBForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.CustomGBForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for CustomGBForce') 
        force = mm.CustomGBForce()
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perParticleParams:
            force.addPerParticleParameter(param)
        for value in self.computedValues:
            force.addComputedValue(value[0], value[1], value[2])
        for term in self.energyTerms:
            force.addEnergyTerm(term[0], term[1])
        for function in self.functions:
            force.addFunction(function[0], function[1], function[2], function[3])
        for atom in data.atoms:
            t = data.atomType[atom]
            if t in self.typeMap:
                values = self.typeMap[t]
                force.addParticle(self.typeMap[t])
            else:
                raise ValueError('No CustomGB parameters defined for atom type '+t) 
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        sys.addForce(force)

parsers["CustomGBForce"] = CustomGBGenerator.parseElement

def getAtomPrint(data, atomIndex):

    if (atomIndex < len(data.atoms)):
        atom = data.atoms[atomIndex]
        returnString = "%4s %4s %5d" % (atom.name, atom.residue.name, atom.residue.index)
    else:
        returnString = "NA"

    return returnString

#=============================================================================================

def countConstraint(data):

    bondCount = 0
    angleCount = 0
    for bond in data.bonds:
        if bond.isConstrained:
            bondCount += 1

    angleCount = 0
    for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):
        if (isConstrained):
            angleCount += 1
 
    print "Constraints bond=%d angle=%d  total=%d" % (bondCount, angleCount, (bondCount+angleCount))

## @private
class AmoebaBondGenerator:

    #=============================================================================================

    """An AmoebaBondGenerator constructs a AmoebaBondForce."""

    #=============================================================================================
    
    def __init__(self, cubic, quartic):

        self.cubic = cubic
        self.quartic = quartic
        self.types1 = []
        self.types2 = []
        self.length = []
        self.k = []
    
    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        # <AmoebaBondForce bond-cubic="-25.5" bond-quartic="379.3125">
        # <Bond class1="1" class2="2" length="0.1437" k="156900.0"/>
    
        generator = AmoebaBondGenerator(float(element.attrib['bond-cubic']), float(element.attrib['bond-quartic']))
        forceField._forces.append(generator)
        for bond in element.findall('Bond'):
            types = forceField._findAtomTypes(bond, 2)
            if types is not None:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.length.append(float(bond.attrib['length']))
                generator.k.append(float(bond.attrib['k']))
            else:
                outputString = "AmoebaBondGenerator: error getting types: %s %s" % (
                                    bond.attrib['class1'],
                                    bond.attrib['class2'])
                raise ValueError(outputString) 
    
    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        #countConstraint(data)

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaBondForce]
        if len(existing) == 0:
            force = mm.AmoebaBondForce()
            sys.addForce(force)
        else:
            force = existing[0]

        force.setAmoebaGlobalBondCubic(self.cubic)
        force.setAmoebaGlobalBondQuartic(self.quartic)

        for bond in data.bonds:
            type1 = data.atomType[data.atoms[bond.atom1]]
            type2 = data.atomType[data.atoms[bond.atom2]]
            hit = 0
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
                    bond.length = self.length[i]
                    hit = 1
                    if bond.isConstrained:
                        sys.addConstraint(bond.atom1, bond.atom2, self.length[i])
                    elif self.k[i] != 0:
                        force.addBond(bond.atom1, bond.atom2, self.length[i], self.k[i])
                    break

            if (hit == 0): 
                outputString = "AmoebaBondGenerator missing: types=[%5s %5s] atoms=[%6d %6d] " % (type1, type2, bond.atom1, bond.atom2)
                raise ValueError(outputString) 

parsers["AmoebaBondForce"] = AmoebaBondGenerator.parseElement

#=============================================================================================
# Add angle constraint
#=============================================================================================
    
def addAngleConstraint(angle, idealAngle, data, sys):

    # Find the two bonds that make this angle.
                    
    bond1 = None
    bond2 = None
    for bond in data.atomBonds[angle[1]]:
        atom1 = data.bonds[bond].atom1
        atom2 = data.bonds[bond].atom2
        if atom1 == angle[0] or atom2 == angle[0]:
            bond1 = bond
        elif atom1 == angle[2] or atom2 == angle[2]:
            bond2 = bond
                    
        # Compute the distance between atoms and add a constraint
                    
        if bond1 is not None and bond2 is not None:
            l1 = data.bonds[bond1].length
            l2 = data.bonds[bond2].length
            if l1 is not None and l2 is not None:
                length = sqrt(l1*l1 + l2*l2 - 2*l1*l2*cos(idealAngle))
                sys.addConstraint(angle[0], angle[2], length)
                return

#=============================================================================================
## @private
class AmoebaAngleGenerator:

    #=============================================================================================
    """An AmoebaAngleGenerator constructs a AmoebaAngleForce."""
    #=============================================================================================
    
    def __init__(self, forceField, cubic, quartic, pentic, sextic):

        self.forceField = forceField
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic

        self.types1 = []
        self.types2 = []
        self.types3 = []

        self.angle = []
        self.k = []
    
    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        # <AmoebaAngleForce angle-cubic="-0.014" angle-quartic="5.6e-05" angle-pentic="-7e-07" angle-sextic="2.2e-08">
        #   <Angle class1="2" class2="1" class3="3" k="0.0637259642196" angle1="122.00"  />

        generator = AmoebaAngleGenerator(forceField, float(element.attrib['angle-cubic']), float(element.attrib['angle-quartic']),  float(element.attrib['angle-pentic']), float(element.attrib['angle-sextic']))
        forceField._forces.append(generator)
        for angle in element.findall('Angle'):
            types = forceField._findAtomTypes(angle, 3)
            if types is not None:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

                angleList = []
                angleList.append(float(angle.attrib['angle1']))

                try:
                    angleList.append(float(angle.attrib['angle2']))
                    try:
                        angleList.append(float(angle.attrib['angle3']))
                    except:
                        pass
                except:
                    pass
                generator.angle.append(angleList)
                generator.k.append(float(angle.attrib['k']))
            else:
                outputString = "AmoebaAngleGenerator: error getting types: %s %s %s" % (
                                    angle.attrib['class1'],
                                    angle.attrib['class2'],
                                    angle.attrib['class3'])
                raise ValueError(outputString) 
    
    #=============================================================================================
    # createForce is bypassed here since the AmoebaOutOfPlaneBendForce generator must first execute
    # and partition angles into in-plane and non-in-plane angles
    #=============================================================================================
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        pass

    #=============================================================================================
    # createForcePostOpBendAngle is called by AmoebaOutOfPlaneBendForce with the list of
    # non-in-plane angles
    #=============================================================================================
    
    def createForcePostOpBendAngle(self, sys, data, nonbondedMethod, nonbondedCutoff, angleList, args):

        # get force

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaAngleForce]

        if len(existing) == 0:
            force = mm.AmoebaAngleForce()
            sys.addForce(force)
        else:
            force = existing[0]

        # set scalars

        force.setAmoebaGlobalAngleCubic(self.cubic)
        force.setAmoebaGlobalAngleQuartic(self.quartic)
        force.setAmoebaGlobalAnglePentic(self.pentic)
        force.setAmoebaGlobalAngleSextic(self.sextic)

        for angleDict in angleList:
            angle = angleDict['angle']
            isConstrained = angleDict['isConstrained']

            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]
            hit = 0
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                    hit = 1
                    if isConstrained and self.k[i] != 0.0:
                        angleDict['idealAngle'] = self.angle[i][0]
                        addAngleConstraint(angle, self.angle[i][0], data, sys)
                    elif self.k[i] != 0:
                        lenAngle = len(self.angle[i])
                        if (lenAngle > 1):
                            # get k-index by counting number of non-angle hydrogens on the central atom
                            # based on kangle.f
                            numberOfHydrogens = 0
                            for bond in data.atomBonds[angle[1]]:
                                atom1 = data.bonds[bond].atom1
                                atom2 = data.bonds[bond].atom2
                                if (atom1 == angle[1] and atom2 != angle[0] and atom2 != angle[2] and (sys.getParticleMass(atom2)/unit.dalton) < 1.90):
                                    numberOfHydrogens += 1
                                if (atom2 == angle[1] and atom1 != angle[0] and atom1 != angle[2] and (sys.getParticleMass(atom1)/unit.dalton) < 1.90):
                                    numberOfHydrogens += 1
                            if (numberOfHydrogens < lenAngle):
                                angleValue =  self.angle[i][numberOfHydrogens]
                            else:
                                outputString = "AmoebaAngleGenerator angle index=%d is out of range: [0, %5d] " % (numberOfHydrogens, lenAngle)
                                raise ValueError(outputString) 
                        else:
                            angleValue =  self.angle[i][0]
               
                        angleDict['idealAngle'] = angleValue
                        force.addAngle(angle[0], angle[1], angle[2], angleValue, self.k[i])
                    break
            if (hit == 0): 
                outputString = "AmoebaAngleGenerator missing types: [%s %s %s] for atoms: " % (type1, type2, type3)
                outputString += getAtomPrint( data, angle[0] ) + ' '
                outputString += getAtomPrint( data, angle[1] ) + ' '
                outputString += getAtomPrint( data, angle[2] )
                outputString += " indices: [%6d %6d %6d]" % (angle[0], angle[1], angle[2])
                raise ValueError(outputString) 

    #=============================================================================================
    # createForcePostOpBendInPlaneAngle is called by AmoebaOutOfPlaneBendForce with the list of
    # in-plane angles
    #=============================================================================================
    
    def createForcePostOpBendInPlaneAngle(self, sys, data, nonbondedMethod, nonbondedCutoff, angleList, args):

        # get force

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaInPlaneAngleForce]

        if len(existing) == 0:
            force = mm.AmoebaInPlaneAngleForce()
            sys.addForce(force)
        else:
            force = existing[0]

        # scalars

        force.setAmoebaGlobalInPlaneAngleCubic(self.cubic)
        force.setAmoebaGlobalInPlaneAngleQuartic(self.quartic)
        force.setAmoebaGlobalInPlaneAnglePentic(self.pentic)
        force.setAmoebaGlobalInPlaneAngleSextic(self.sextic)

        for angleDict in angleList:
 
            angle = angleDict['angle']
            isConstrained = angleDict['isConstrained']

            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]

            hit = 0
            for i in range(len(self.types1)):

                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]

                if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                    hit = 1
                    angleDict['idealAngle'] = self.angle[i][0]
                    if (isConstrained and self.k[i] != 0.0):
                        addAngleConstraint(angle, self.angle[i][0], data, sys)
                    else:
                        force.addAngle(angle[0], angle[1], angle[2], angle[3], self.angle[i][0], self.k[i])
                    break

            if (hit == 0): 
                outputString = "AmoebaInPlaneAngleGenerator missing types: [%s %s %s] atoms: " % (type1, type2, type3)
                outputString += getAtomPrint( data, angle[0] ) + ' '
                outputString += getAtomPrint( data, angle[1] ) + ' '
                outputString += getAtomPrint( data, angle[2] )
                outputString += " indices: [%6d %6d %6d]" % (angle[0], angle[1], angle[2])
                raise ValueError(outputString) 

parsers["AmoebaAngleForce"] = AmoebaAngleGenerator.parseElement

#=============================================================================================
# Generator for the AmoebaOutOfPlaneBend covalent force; also calls methods in the
# AmoebaAngleGenerator to generate the AmoebaAngleForce and
# AmoebaInPlaneAngleForce
#=============================================================================================

## @private
class AmoebaOutOfPlaneBendGenerator:

    #=============================================================================================

    """An AmoebaOutOfPlaneBendGenerator constructs a AmoebaOutOfPlaneBendForce."""
    
    #=============================================================================================

    def __init__(self, forceField, type, cubic, quartic, pentic, sextic):

        self.forceField = forceField
        self.type = type
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.types4 = []

        self.ks = []

    #=============================================================================================
    # Local version of findAtomTypes needed since class indices are 0 (i.e., not recognized)
    # for types3 and 4
    #=============================================================================================
    
    def findAtomTypes(self, forceField, node, num):
        """Parse the attributes on an XML tag to find the set of atom types for each atom it involves."""
        types = []
        attrib = node.attrib
        for i in range(num):
            if num == 1:
                suffix = ''
            else:
                suffix = str(i+1)
            classAttrib = 'class'+suffix
            if classAttrib in attrib:
                if attrib[classAttrib] in forceField._atomClasses:
                    types.append(forceField._atomClasses[attrib[classAttrib]])
                else:
                    types.append(set())
        return types

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaOutOfPlaneBendForce type="ALLINGER" opbend-cubic="-0.014" opbend-quartic="5.6e-05" opbend-pentic="-7e-07" opbend-sextic="2.2e-08">
        #   <Angle class1="2" class2="1" class3="0" class4="0" k="0.0531474541591"/>
        #   <Angle class1="3" class2="1" class3="0" class4="0" k="0.0898536095496"/>
         
        # get global scalar parameters

        generator = AmoebaOutOfPlaneBendGenerator(forceField, element.attrib['type'],
                                                   float(element.attrib['opbend-cubic']),
                                                   float(element.attrib['opbend-quartic']),
                                                   float(element.attrib['opbend-pentic']),
                                                   float(element.attrib['opbend-sextic']))

        forceField._forces.append(generator)

        for angle in element.findall('Angle'):
            types = generator.findAtomTypes(forceField, angle, 4)
            if types is not None:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])
                generator.types4.append(types[3])

                generator.ks.append(float(angle.attrib['k']))

            else:
                outputString = "AmoebaOutOfPlaneBendGenerator error getting types: %s %s %s %s." % (
                               angle.attrib['class1'], angle.attrib['class2'], angle.attrib['class3'], angle.attrib['class4'])
                raise ValueError(outputString) 
    
    #=============================================================================================
    # Get middle atom in a angle
    # return index of middle atom or -1 if no middle is found
    # This method appears not to be needed since the angle[1] entry appears to always
    # be the middle atom. However, was unsure if this is guaranteed
    #=============================================================================================
    
    def getMiddleAtom(self, angle, data):

        # find atom shared by both bonds making up the angle

        middleAtom = -1
        for atomIndex in angle: 
            isMiddle = 0
            for bond in data.atomBonds[atomIndex]:
                atom1 = data.bonds[bond].atom1
                atom2 = data.bonds[bond].atom2
                if (atom1 != atomIndex):
                    partner = atom1
                else:
                    partner = atom2
                if (partner == angle[0] or partner == angle[1] or partner == angle[2]): 
                    isMiddle += 1

            if (isMiddle == 2):
                return atomIndex
        return -1

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        # get force

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaOutOfPlaneBendForce]
        if len(existing) == 0:
            force = mm.AmoebaOutOfPlaneBendForce()
            sys.addForce(force)
        else:
            force = existing[0]

        # set scalars

        force.setAmoebaGlobalOutOfPlaneBendCubic(  self.cubic)
        force.setAmoebaGlobalOutOfPlaneBendQuartic(self.quartic)
        force.setAmoebaGlobalOutOfPlaneBendPentic( self.pentic)
        force.setAmoebaGlobalOutOfPlaneBendSextic( self.sextic)

        # this hash is used to insure the out-of-plane-bend bonds
        # are only added once

        skipAtoms = dict()

        # these lists are used in the partitioning of the angles into
        # angle and inPlane angles

        inPlaneAngles = []
        nonInPlaneAngles = []
        nonInPlaneAnglesConstrained = []
        idealAngles = []*len(data.angles)

        for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):

            middleAtom = self.getMiddleAtom(angle, data)
            if (middleAtom > -1):
                middleType = data.atomType[data.atoms[middleAtom]]
                middleCovalency = len(data.atomBonds[middleAtom])
            else:
                middleType = -1
                middleCovalency = -1

            # if middle atom has covalency of 3 and 
            # the types of the middle atom and the partner atom (atom bonded to
            # middle atom, but not in angle) match types1 and types2, then
            # three out-of-plane bend angles are generated. Three in-plane angle 
            # are also generated. If the conditions are not satisfied, the angle is marked as 'generic' angle (not a in-plane angle)

            if (middleAtom > -1 and middleCovalency == 3 and middleAtom not in skipAtoms):

                partners = []
                partnerSet = set()
                partnerTypes = []
                partnerK = []

                for bond in data.atomBonds[middleAtom]:
                    atom1 = data.bonds[bond].atom1
                    atom2 = data.bonds[bond].atom2
                    if (atom1 != middleAtom):
                        partner = atom1
                    else:
                        partner = atom2

                    partnerType = data.atomType[data.atoms[partner]]
                    for i in range(len(self.types1)):
                        types1 = self.types1[i]
                        types2 = self.types2[i]
                        if (middleType in types2 and partnerType in types1):
                            partners.append(partner)
                            partnerSet.add(partner)
                            partnerTypes.append(partnerType)
                            partnerK.append(self.ks[i])
             
                if (len(partners) == 3):

                    force.addOutOfPlaneBend(partners[0], middleAtom, partners[1], partners[2], partnerK[2])
                    force.addOutOfPlaneBend(partners[0], middleAtom, partners[2], partners[1], partnerK[1])
                    force.addOutOfPlaneBend(partners[1], middleAtom, partners[2], partners[0], partnerK[0])

                    # skipAtoms is used to insure angles are only included once

                    skipAtoms[middleAtom] = set()
                    skipAtoms[middleAtom].add(partners[0])
                    skipAtoms[middleAtom].add(partners[1])
                    skipAtoms[middleAtom].add(partners[2])

                    # in-plane angle

                    angleDict = {}
                    angleList = []
                    angleList.append(angle[0])
                    angleList.append(angle[1])
                    angleList.append(angle[2])
                    angleDict['angle'] = angleList

                    angleDict['isConstrained'] = 0

                    angleSet = set()
                    angleSet.add(angle[0])
                    angleSet.add(angle[1])
                    angleSet.add(angle[2])

                    for atomIndex in partnerSet:
                        if (atomIndex not in angleSet):
                            angleList.append(atomIndex)

                    inPlaneAngles.append(angleDict)

                else:
                    angleDict = {}
                    angleDict['angle'] = angle
                    angleDict['isConstrained'] = isConstrained
                    nonInPlaneAngles.append(angleDict)
            else:
                if (middleAtom > -1 and middleCovalency == 3 and middleAtom in skipAtoms):

                    partnerSet = skipAtoms[middleAtom]
                  
                    angleDict = {}

                    angleList = []
                    angleList.append(angle[0])
                    angleList.append(angle[1])
                    angleList.append(angle[2])
                    angleDict['angle'] = angleList

                    angleDict['isConstrained'] = isConstrained

                    angleSet = set()
                    angleSet.add(angle[0])
                    angleSet.add(angle[1])
                    angleSet.add(angle[2])

                    for atomIndex in partnerSet:
                        if (atomIndex not in angleSet):
                            angleList.append(atomIndex)

                    inPlaneAngles.append(angleDict)

                else:
                    angleDict = {}
                    angleDict['angle'] = angle
                    angleDict['isConstrained'] = isConstrained
                    nonInPlaneAngles.append(angleDict)

        # get AmoebaAngleGenerator and add AmoebaAngle and AmoebaInPlaneAngle forces

        for force in self.forceField._forces:
            if (force.__class__.__name__ == 'AmoebaAngleGenerator'): 
                force.createForcePostOpBendAngle(sys, data, nonbondedMethod, nonbondedCutoff, nonInPlaneAngles, args)
                force.createForcePostOpBendInPlaneAngle(sys, data, nonbondedMethod, nonbondedCutoff, inPlaneAngles, args)

        for force in self.forceField._forces:
            if (force.__class__.__name__ == 'AmoebaStretchBendGenerator'): 
                for angleDict in inPlaneAngles:
                    nonInPlaneAngles.append(angleDict)
                force.createForcePostAmoebaBondForce(sys, data, nonbondedMethod, nonbondedCutoff, nonInPlaneAngles, args)

parsers["AmoebaOutOfPlaneBendForce"] = AmoebaOutOfPlaneBendGenerator.parseElement

#=============================================================================================

## @private
class AmoebaTorsionGenerator:

    #=============================================================================================
    """An AmoebaTorsionGenerator constructs a AmoebaTorsionForce."""
    #=============================================================================================

    def __init__(self, torsionUnit):

        self.torsionUnit = torsionUnit

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.types4 = []

        self.t1 = []
        self.t2 = []
        self.t3 = []
    
    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaTorsionForce torsionUnit="0.5">
        #   <Torsion class1="3" class2="1" class3="2" class4="3"   amp1="0.0" angle1="0.0"   amp2="0.0" angle2="3.14159265359"   amp3="0.0" angle3="0.0" />
        #   <Torsion class1="3" class2="1" class3="2" class4="6"   amp1="0.0" angle1="0.0"   amp2="0.0" angle2="3.14159265359"   amp3="-0.263592" angle3="0.0" />
         
        generator = AmoebaTorsionGenerator(float(element.attrib['torsionUnit']))
        forceField._forces.append(generator)

        # collect particle classes and t1,t2,t3,
        # where ti=[amplitude_i,angle_i]

        for torsion in element.findall('Torsion'):
            types = forceField._findAtomTypes(torsion, 4)
            if types is not None:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])
                generator.types4.append(types[3])

                for ii in range(1,4):
                    tInfo = []
                    suffix = str(ii)
                    ampName = 'amp' + suffix
                    tInfo.append(float(torsion.attrib[ampName]))

                    angName = 'angle' + suffix
                    tInfo.append(float(torsion.attrib[angName]))

                    if (ii == 1):
                        generator.t1.append(tInfo)
                    elif (ii == 2):
                        generator.t2.append(tInfo)
                    elif (ii == 3):
                        generator.t3.append(tInfo)

            else:
                outputString = "AmoebaTorsionGenerator: error getting types: %s %s %s %s" % (
                                    stretchBend.attrib['class1'],
                                    stretchBend.attrib['class2'],
                                    stretchBend.attrib['class3'],
                                    stretchBend.attrib['class4'])
                raise ValueError(outputString) 
    
    #=============================================================================================

    def createForce(self, sys, data, nontorsionedMethod, nontorsionedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.PeriodicTorsionForce]
        if len(existing) == 0:
            force = mm.PeriodicTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for torsion in data.propers:

            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]

            hit = 0
            for i in range(len(self.types1)):

                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                types4 = self.types4[i]

                # match types in forward or reverse direction

                if (type1 in types1 and type2 in types2 and type3 in types3 and type4 in types4) or (type4 in types1 and type3 in types2 and type2 in types3 and type1 in types4):
                    hit = 1
                    if self.t1[i][0] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], 1, self.t1[i][1], self.t1[i][0])
                    if self.t2[i][0] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], 2, self.t2[i][1], self.t2[i][0])
                    if self.t3[i][0] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], 3, self.t3[i][1], self.t3[i][0])
                    break

            if (hit == 0): 
                outputString = "AmoebaTorsionGenerator missing type: [%s %s %s %s] atoms: " % (type1, type2, type3, type4)
                outputString += getAtomPrint( data, torsion[0] ) + ' '
                outputString += getAtomPrint( data, torsion[1] ) + ' '
                outputString += getAtomPrint( data, torsion[2] ) + ' '
                outputString += getAtomPrint( data, torsion[3] )
                outputString += " indices: [%6d %6d %6d %6d]" % (torsion[0], torsion[1], torsion[2], torsion[3])
                raise ValueError(outputString) 

parsers["AmoebaTorsionForce"] = AmoebaTorsionGenerator.parseElement

#=============================================================================================

## @private
class AmoebaPiTorsionGenerator:

    #=============================================================================================

    """An AmoebaPiTorsionGenerator constructs a AmoebaPiTorsionForce."""

    #=============================================================================================
    
    def __init__(self, piTorsionUnit):
        self.piTorsionUnit = piTorsionUnit 
        self.types1 = []
        self.types2 = []
        self.k = []
    
    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaPiTorsionForce piTorsionUnit="1.0">
        #   <PiTorsion class1="1" class2="3" k="28.6604" />

        generator = AmoebaPiTorsionGenerator(float(element.attrib['piTorsionUnit']))
        forceField._forces.append(generator)

        for piTorsion in element.findall('PiTorsion'):
            types = forceField._findAtomTypes(piTorsion, 2)
            if types is not None:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.k.append(float(piTorsion.attrib['k']))
            else:
                outputString = "AmoebaPiTorsionGenerator: error getting types: %s %s " % (
                                    piTorsion.attrib['class1'],
                                    piTorsion.attrib['class2'])
                raise ValueError(outputString) 
    
    #=============================================================================================

    def createForce(self, sys, data, nonpiTorsionedMethod, nonpiTorsionedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaPiTorsionForce]

        if len(existing) == 0:
            force = mm.AmoebaPiTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for bond in data.bonds:

            # search for bonds with both atoms in bond having covalency == 3
 
            atom1 = bond.atom1
            atom2 = bond.atom2
 
            if (len(data.atomBonds[atom1]) == 3 and len(data.atomBonds[atom1]) == 3):

                type1 = data.atomType[data.atoms[atom1]]
                type2 = data.atomType[data.atoms[atom2]]

                for i in range(len(self.types1)):

                   types1 = self.types1[i]
                   types2 = self.types2[i]

                   if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):

                       # piTorsionAtom1, piTorsionAtom2 are the atoms bonded to atom1, excluding atom2 
                       # piTorsionAtom5, piTorsionAtom6 are the atoms bonded to atom2, excluding atom1 

                       piTorsionAtom1 = -1
                       piTorsionAtom2 = -1
                       piTorsionAtom3 = atom1

                       piTorsionAtom4 = atom2
                       piTorsionAtom5 = -1
                       piTorsionAtom6 = -1

                       for bond in data.atomBonds[atom1]:
                           bondedAtom1 = data.bonds[bond].atom1
                           bondedAtom2 = data.bonds[bond].atom2
                           if (bondedAtom1 != atom1):
                               b1 = bondedAtom1
                           else:
                               b1 = bondedAtom2
                           if (b1 != atom2):
                               if (piTorsionAtom1 == -1):
                                   piTorsionAtom1 = b1 
                               else:
                                   piTorsionAtom2 = b1

                       for bond in data.atomBonds[atom2]:
                           bondedAtom1 = data.bonds[bond].atom1
                           bondedAtom2 = data.bonds[bond].atom2
                           if (bondedAtom1 != atom2):
                               b1 = bondedAtom1
                           else:
                               b1 = bondedAtom2

                           if (b1 != atom1):
                               if (piTorsionAtom5 == -1):
                                   piTorsionAtom5 = b1 
                               else:
                                   piTorsionAtom6 = b1
    
                       force.addPiTorsion(piTorsionAtom1, piTorsionAtom2, piTorsionAtom3, piTorsionAtom4, piTorsionAtom5, piTorsionAtom6, self.k[i])

parsers["AmoebaPiTorsionForce"] = AmoebaPiTorsionGenerator.parseElement

#=============================================================================================

## @private
class AmoebaTorsionTorsionGenerator:

    #=============================================================================================
    """An AmoebaTorsionTorsionGenerator constructs a AmoebaTorsionTorsionForce."""
    #=============================================================================================

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.types4 = []
        self.types5 = []

        self.gridIndex = []

        self.grids = []
    
    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        generator = AmoebaTorsionTorsionGenerator()
        forceField._forces.append(generator)
        maxGridIndex = -1

        # <AmoebaTorsionTorsionForce >
        # <TorsionTorsion class1="3" class2="1" class3="2" class4="3" class5="1" grid="0" nx="25" ny="25" />

        for torsionTorsion in element.findall('TorsionTorsion'):
            types = forceField._findAtomTypes(torsionTorsion, 5)
            if types is not None:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])
                generator.types4.append(types[3])
                generator.types5.append(types[4])

                gridIndex = int(torsionTorsion.attrib['grid'])
                if (gridIndex > maxGridIndex):
                    maxGridIndex = gridIndex

                generator.gridIndex.append(gridIndex)
            else:
                outputString = "AmoebaTorsionTorsionGenerator: error getting types: %s %s %s %s %s" % (
                                    torsionTorsion.attrib['class1'],
                                    torsionTorsion.attrib['class2'],
                                    torsionTorsion.attrib['class3'],
                                    torsionTorsion.attrib['class4'],
                                    torsionTorsion.attrib['class5'] )
                raise ValueError(outputString) 
    
        # load grid

        # xml source

        # <TorsionTorsionGrid grid="0" nx="25" ny="25" >
        # <Grid angle1="-180.00" angle2="-180.00" f="0.0" fx="2.31064374824e-05" fy="0.0" fxy="-0.0052801799672" />
        # <Grid angle1="-165.00" angle2="-180.00" f="-0.66600912" fx="-0.06983370052" fy="-0.075058725744" fxy="-0.0044462732032" />

        # output grid:

        #     grid[x][y][0] = x value 
        #     grid[x][y][1] = y value 
        #     grid[x][y][2] = function value 
        #     grid[x][y][3] = dfdx value 
        #     grid[x][y][4] = dfdy value 
        #     grid[x][y][5] = dfd(xy) value 

        maxGridIndex    += 1
        generator.grids = maxGridIndex*[]
        for torsionTorsionGrid in element.findall('TorsionTorsionGrid'):

            gridIndex = int(torsionTorsionGrid.attrib[ "grid"])
            nx = int(torsionTorsionGrid.attrib[ "nx"])
            ny = int(torsionTorsionGrid.attrib[ "ny"])

            grid = []
            gridCol = []

            gridColIndex = 0

            for gridEntry in torsionTorsionGrid.findall('Grid'):

                gridRow = []
                gridRow.append(float(gridEntry.attrib['angle1']))
                gridRow.append(float(gridEntry.attrib['angle2']))
                gridRow.append(float(gridEntry.attrib['f']))
                gridRow.append(float(gridEntry.attrib['fx']))
                gridRow.append(float(gridEntry.attrib['fy']))
                gridRow.append(float(gridEntry.attrib['fxy']))
                gridCol.append(gridRow)

                gridColIndex  += 1
                if (gridColIndex == nx):
                    grid.append(gridCol)
                    gridCol = []
                    gridColIndex = 0

            
            if (gridIndex == len(generator.grids)):
                generator.grids.append(grid)
            else:
                while(len(generator.grids) < gridIndex):
                    generator.grids.append([])
                generator.grids[gridIndex] = grid

    #=============================================================================================

    def getChiralAtomIndex(self, data, sys, atomB, atomC, atomD):

        chiralAtomIndex = -1

        # if atomC has four bonds, find the
        # two bonds that do not include atomB and atomD
        # set chiralAtomIndex to one of these, if they are
        # not the same atom(type/mass)

        if (len(data.atomBonds[atomC]) == 4):
            atomE = -1
            atomF = -1
            for bond in data.atomBonds[atomC]:
                bondedAtom1 = data.bonds[bond].atom1
                bondedAtom2 = data.bonds[bond].atom2
                hit = -1
                if (  bondedAtom1 == atomC and bondedAtom2 != atomB and bondedAtom2 != atomD):
                    hit = bondedAtom2
                elif (bondedAtom2 == atomC and bondedAtom1 != atomB and bondedAtom1 != atomD):
                    hit = bondedAtom1

                if (hit > -1):
                    if (atomE == -1):
                        atomE = hit
                    else:
                        atomF = hit
       
            # raise error if atoms E or F not found

            if (atomE == -1 or atomF == -1):
                outputString = "getChiralAtomIndex: error getting bonded partners of atomC=%s %d %s" % (atomC.name, atomC.resiude.index, atomC.resiude.name,)
                raise ValueError(outputString) 

            # check for different type/mass between atoms E & F

            typeE = int(data.atomType[data.atoms[atomE]])
            typeF = int(data.atomType[data.atoms[atomF]])
            if (typeE > typeF):
                chiralAtomIndex = atomE 
            if (typeF > typeE):
                chiralAtomIndex = atomF 

            massE = sys.getParticleMass(atomE)/unit.dalton
            massF = sys.getParticleMass(atomE)/unit.dalton
            if (massE > massF):
                chiralAtomIndex = massE 
            if (massF > massE):
                chiralAtomIndex = massF 

        return chiralAtomIndex

    #=============================================================================================
 
    def createForce(self, sys, data, nonpiTorsionedMethod, nonpiTorsionedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaTorsionTorsionForce]

        if len(existing) == 0:
            force = mm.AmoebaTorsionTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]

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
                    for bondIndex in data.atomBonds[id]:
                        bondedAtom1 = data.bonds[bondIndex].atom1
                        bondedAtom2 = data.bonds[bondIndex].atom2
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

                            for i in range(len(self.types1)):

                                types1 = self.types1[i]
                                types2 = self.types2[i]
                                types3 = self.types3[i]
                                types4 = self.types4[i]
                                types5 = self.types5[i]

                                # match in order

                                if (type1 in types1 and type2 in types2 and type3 in types3 and type4 in types4 and type5 in types5):
                                    chiralAtomIndex = self.getChiralAtomIndex(data, sys, ib, ic, id)
                                    force.addTorsionTorsion(ia, ib, ic, id, ie, chiralAtomIndex, self.gridIndex[i])

                                # match in reverse order

                                if (type5 in types1 and type4 in types2 and type3 in types3 and type2 in types4 and type1 in types5):
                                    chiralAtomIndex = self.getChiralAtomIndex(data, sys, ib, ic, id)
                                    force.addTorsionTorsion(ie, id, ic, ib, ia, chiralAtomIndex, self.gridIndex[i])

        # set grids

        for (index, grid) in enumerate(self.grids):
            force.setTorsionTorsionGrid(index, grid)
 
parsers["AmoebaTorsionTorsionForce"] = AmoebaTorsionTorsionGenerator.parseElement

#=============================================================================================

## @private
class AmoebaStretchBendGenerator:

    #=============================================================================================
    """An AmoebaStretchBendGenerator constructs a AmoebaStretchBendForce."""
    #=============================================================================================

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

        self.k1 = []
        self.k2 = []
    
    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):
        generator = AmoebaStretchBendGenerator()
        forceField._forces.append(generator)

        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />

        for stretchBend in element.findall('StretchBend'):
            types = forceField._findAtomTypes(stretchBend, 3)
            if types is not None:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

                generator.k1.append(float(stretchBend.attrib['k1']))
                generator.k2.append(float(stretchBend.attrib['k2']))

            else:
                outputString = "AmoebaStretchBendGenerator : error getting types: %s %s %s" % (
                                    stretchBend.attrib['class1'],
                                    stretchBend.attrib['class2'],
                                    stretchBend.attrib['class3'])
                raise ValueError(outputString) 
    
    #=============================================================================================

    # The setup of this force is dependent on AmoebaBondForce and AmoebaAngleForce 
    # having been called since the ideal bond lengths and angle are needed here.
    # As a conseqeunce, createForce() is not implemented since it is not guaranteed that the generator for
    # AmoebaBondForce and AmoebaAngleForce have been called prior to AmoebaStretchBendGenerator(). 
    # Instead, createForcePostAmoebaBondForce() is called 
    # after the generators for AmoebaBondForce and AmoebaAngleForce have been called

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        pass

    #=============================================================================================

    # Note: request for constrained bonds is ignored.

    #=============================================================================================

    def createForcePostAmoebaBondForce(self, sys, data, nonbondedMethod, nonbondedCutoff, angleList, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaStretchBendForce]
        if len(existing) == 0:
            force = mm.AmoebaStretchBendForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for angleDict in angleList:

            angle = angleDict['angle']
            if ('isConstrained' in angleDict):
                isConstrained = angleDict['isConstrained']
            else:
                isConstrained = 0

            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]

            radian = 57.2957795130
            for i in range(len(self.types1)):

                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]

                # match types
                # get ideal bond lengths, bondAB, bondCB
                # get ideal angle

                if (type2 in types2 and ((type1 in types1 and type3 in types3) or (type3 in types1 and type1 in types3))): 
                    bondAB = -1.0
                    bondCB = -1.0
                    swap = 0
                    for bond in data.atomBonds[angle[1]]:
                        atom1 = data.bonds[bond].atom1
                        atom2 = data.bonds[bond].atom2
                        length = data.bonds[bond].length
                        if (atom1 == angle[0]):
                            bondAB = length
                        if (atom1 == angle[2]):
                            bondCB = length
                        if (atom2 == angle[2]):
                            bondCB = length
                        if (atom2 == angle[0]):
                            bondAB = length
                    
                    # check that ideal angle and bonds are set

                    if ('idealAngle' not in angleDict):

                       outputString = "AmoebaStretchBendGenerator: ideal angle is not set for following entry:\n"
                       outputString += "   types: %5s %5s %5s atoms: " % (type1, type2, type3)
                       outputString += getAtomPrint( data, angle[0] ) + ' '
                       outputString += getAtomPrint( data, angle[1] ) + ' '
                       outputString += getAtomPrint( data, angle[2] )
                       raise ValueError(outputString) 

                    elif (bondAB < 0 or bondCB < 0):

                       outputString = "AmoebaStretchBendGenerator: bonds not set: %15.7e %15.7e. for following entry:" % (bondAB, bondCB)
                       outputString += "     types: [%5s %5s %5s] atoms: " % (type1, type2, type3)
                       outputString += getAtomPrint( data, angle[0] ) + ' '
                       outputString += getAtomPrint( data, angle[1] ) + ' '
                       outputString += getAtomPrint( data, angle[2] )
                       raise ValueError(outputString) 

                    else:
                        force.addStretchBend(angle[0], angle[1], angle[2], bondAB, bondCB, angleDict['idealAngle']/radian, self.k1[i])

                    break

parsers["AmoebaStretchBendForce"] = AmoebaStretchBendGenerator.parseElement

#=============================================================================================

## @private
class AmoebaVdwGenerator:

    """A AmoebaVdwGenerator constructs a AmoebaVdwForce."""
    
    #=============================================================================================

    def __init__(self, type, radiusrule, radiustype, radiussize, epsilonrule, vdw13Scale, vdw14Scale, vdw15Scale):

        self.type = type 

        self.radiusrule = radiusrule
        self.radiustype = radiustype
        self.radiussize = radiussize

        self.epsilonrule = epsilonrule

        self.vdw13Scale = vdw13Scale
        self.vdw14Scale = vdw14Scale
        self.vdw15Scale = vdw15Scale

        self.typeMap = {}

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        # <AmoebaVdwForce type="BUFFERED-14-7" radiusrule="CUBIC-MEAN" radiustype="R-MIN" radiussize="DIAMETER" epsilonrule="HHG" vdw-13-scale="0.0" vdw-14-scale="1.0" vdw-15-scale="1.0" >
        #   <Vdw class="1" sigma="0.371" epsilon="0.46024" reduction="1.0" /> 
        #   <Vdw class="2" sigma="0.382" epsilon="0.422584" reduction="1.0" /> 
         
        generator = AmoebaVdwGenerator(element.attrib['type'], element.attrib['radiusrule'], element.attrib['radiustype'], element.attrib['radiussize'], element.attrib['epsilonrule'], 
                                        float(element.attrib['vdw-13-scale']), float(element.attrib['vdw-14-scale']), float(element.attrib['vdw-15-scale'])) 
        forceField._forces.append(generator)
        two_six = 1.122462048309372

        # types[] = [ sigma, epsilon, reductionFactor, class ]
        # sigma is modified based on radiustype and radiussize

        for atom in element.findall('Vdw'):
            types = forceField._findAtomTypes(atom, 1)
            if types is not None:

                values = [float(atom.attrib['sigma']), float(atom.attrib['epsilon']), float(atom.attrib['reduction'])]


                if (generator.radiustype == 'SIGMA'):
                    values[0] *= two_six
      
                if (generator.radiussize == 'DIAMETER'):
                    values[0] *= 0.5

                for t in types[0]:
                    generator.typeMap[t] = values
    
            else:
                outputString = "AmoebaVdwGenerator: error getting type: %s" % (atom.attrib['class'])
                raise ValueError(outputString) 
    
    #=============================================================================================

    # Return a set containing the indices of particles bonded to particle with index=particleIndex

    #=============================================================================================

    @staticmethod
    def getBondedParticleSet(particleIndex, data):

        bondedParticleSet = set()

        for bond in data.atomBonds[particleIndex]:
            atom1 = data.bonds[bond].atom1
            atom2 = data.bonds[bond].atom2
            if (atom1 != particleIndex):
                bondedParticleSet.add(atom1)
            else:
                bondedParticleSet.add(atom2)

        return bondedParticleSet
          
    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        sigmaMap = {'ARITHMETIC':1, 'GEOMETRIC':1, 'CUBIC-MEAN':1}
        epsilonMap = {'ARITHMETIC':1, 'GEOMETRIC':1, 'HARMONIC':1, 'HHG':1}

        # get or create force depending on whether it has already been added to the system

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaVdwForce]
        if len(existing) == 0:
            force = mm.AmoebaVdwForce()
            sys.addForce(force)

            # sigma and epsilon combining rules

            if ('sigmaCombiningRule' in args):
                sigmaRule = args['sigmaCombiningRule'].upper()
                if (sigmaRule.upper() in sigmaMap):
                    force.setSigmaCombiningRule(sigmaRule.upper())
                else:
                    stringList = ' ' . join(str(x) for x in sigmaMap.keys())
                    raise ValueError( "AmoebaVdwGenerator: sigma combining rule %s not recognized; valid values are %s; using default." % (sigmaRule, stringList) )
            else:
                force.setSigmaCombiningRule(self.radiusrule)

            if ('epsilonCombiningRule' in args):
                epsilonRule = args['epsilonCombiningRule'].upper()
                if (epsilonRule.upper() in epsilonMap):
                    force.setEpsilonCombiningRule(epsilonRule.upper())
                else:
                    stringList = ' ' . join(str(x) for x in epsilonMap.keys())
                    raise ValueError( "AmoebaVdwGenerator: epsilon combining rule %s not recognized; valid values are %s; using default." % (epsilonRule, stringList) )  
            else:
                force.setEpsilonCombiningRule(self.epsilonrule)
               
            # cutoff

            if ('vdwCutoff' in args):
                force.setCutoff(args['vdwCutoff'])
            else:
                force.setCutoff(nonbondedCutoff)

            # dispersion correction

            if ('useDispersionCorrection' in args):
                force.setUseDispersionCorrection(int(args['useDispersionCorrection']))

            if (nonbondedMethod == PME):
                force.setNonbondedMethod(mm.AmoebaVdwForce.CutoffPeriodic)
           
        else:
            force = existing[0]

        # add particles to force
        # throw error if particle type not available 

        for (i, atom) in enumerate(data.atoms):
            t = data.atomType[atom]
            if t in self.typeMap:

                values = self.typeMap[t]
   
                # ivIndex = index of bonded partner for hydrogens; otherwise ivIndex = particle index

                ivIndex = i
                mass = sys.getParticleMass(i)/unit.dalton
                if (mass < 1.9 and len(data.atomBonds[i]) == 1):
                    bondIndex = data.atomBonds[i][0]
                    if (data.bonds[bondIndex].atom1 == i):
                        ivIndex = data.bonds[bondIndex].atom2
                    else:
                        ivIndex = data.bonds[bondIndex].atom1

                force.addParticle(ivIndex, values[0], values[1], values[2])
            else:
                raise ValueError('No vdw type for atom %s' % (atom.name))

        # set combining rules

        # set particle exclusions: self, 1-2 and 1-3 bonds
        # (1) collect in bondedParticleSets[i], 1-2 indices for all bonded partners of particle i
        # (2) add 1-2,1-3 and self to exclusion set

        bondedParticleSets = []
        for i in range(len(data.atoms)):
            bondedParticleSets.append(AmoebaVdwGenerator.getBondedParticleSet(i, data))

        for (i,atom) in enumerate(data.atoms):
 
            # 1-2 partners

            exclusionSet = bondedParticleSets[i].copy()

            # 1-3 partners

            if (self.vdw13Scale == 0.0):
                for bondedParticle in bondedParticleSets[i]:
                    exclusionSet = exclusionSet.union(bondedParticleSets[bondedParticle])

            # self

            exclusionSet.add(i)

            force.setParticleExclusions(i, exclusionSet)

parsers["AmoebaVdwForce"] = AmoebaVdwGenerator.parseElement

#=============================================================================================

## @private
class AmoebaMultipoleGenerator:

    #=============================================================================================

    """A AmoebaMultipoleGenerator constructs a AmoebaMultipoleForce."""
    
    #=============================================================================================

    def __init__(self, forceField,
                       direct11Scale, direct12Scale, direct13Scale, direct14Scale,
                       mpole12Scale,  mpole13Scale,  mpole14Scale,  mpole15Scale,
                       mutual11Scale, mutual12Scale, mutual13Scale, mutual14Scale,
                       polar12Scale,  polar13Scale,  polar14Scale,  polar15Scale):

        self.forceField = forceField

        self.direct11Scale = direct11Scale 
        self.direct12Scale = direct12Scale 
        self.direct13Scale = direct13Scale 
        self.direct14Scale = direct14Scale 

        self.mpole12Scale = mpole12Scale 
        self.mpole13Scale = mpole13Scale 
        self.mpole14Scale = mpole14Scale 
        self.mpole15Scale = mpole15Scale 

        self.mutual11Scale = mutual11Scale 
        self.mutual12Scale = mutual12Scale 
        self.mutual13Scale = mutual13Scale 
        self.mutual14Scale = mutual14Scale 

        self.polar12Scale = polar12Scale 
        self.polar13Scale = polar13Scale 
        self.polar14Scale = polar14Scale 
        self.polar15Scale = polar15Scale 

        self.typeMap = {}

    #=============================================================================================
    # Set axis type
    #=============================================================================================

    @staticmethod
    def setAxisType(kIndices):

                # set axis type

                kIndicesLen = len(kIndices)
                if (kIndicesLen > 3):
                    ky = kIndices[3]
                else:
                    ky = 0
   
                if (kIndicesLen > 2):
                    kx = kIndices[2]
                else:
                    kx = 0
   
                if (kIndicesLen > 1):
                    kz = kIndices[1]
                else:
                    kz = 0

                while(len(kIndices) < 4):
                    kIndices.append(0)

                axisType = mm.AmoebaMultipoleForce.ZThenX
                if (kz == 0):
                    axisType = mm.AmoebaMultipoleForce.NoAxisType
                if (kz != 0 and kx == 0):
                    axisType = mm.AmoebaMultipoleForce.ZOnly
                if (kz < 0 or kx < 0):
                    axisType = mm.AmoebaMultipoleForce.Bisector
                if (kx < 0 and ky < 0):
                    axisType = mm.AmoebaMultipoleForce.ZBisect
                if (kz < 0 and kx < 0 and ky  < 0):
                    axisType = mm.AmoebaMultipoleForce.ThreeFold

                kIndices[1] = abs(kz)   
                kIndices[2] = abs(kx)   
                kIndices[3] = abs(ky)   

                return axisType

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #   <AmoebaMultipoleForce  direct11Scale="0.0"  direct12Scale="1.0"  direct13Scale="1.0"  direct14Scale="1.0"  mpole12Scale="0.0"  mpole13Scale="0.0"  mpole14Scale="0.4"  mpole15Scale="0.8"  mutual11Scale="1.0"  mutual12Scale="1.0"  mutual13Scale="1.0"  mutual14Scale="1.0"  polar12Scale="0.0"  polar13Scale="0.0"  polar14Intra="0.5"  polar14Scale="1.0"  polar15Scale="1.0"  > 
        # <Multipole class="1"    kz="2"    kx="4"    c0="-0.22620" d1="0.08214" d2="0.00000" d3="0.34883" q11="0.11775" q21="0.00000" q22="-1.02185" q31="-0.17555" q32="0.00000" q33="0.90410"  />
        # <Multipole class="2"    kz="1"    kx="3"    c0="-0.15245" d1="0.19517" d2="0.00000" d3="0.19687" q11="-0.20677" q21="0.00000" q22="-0.48084" q31="-0.01672" q32="0.00000" q33="0.68761"  />

        generator = AmoebaMultipoleGenerator(forceField,
                                              element.attrib['direct11Scale'],
                                              element.attrib['direct12Scale'],
                                              element.attrib['direct13Scale'],
                                              element.attrib['direct14Scale'],

                                              element.attrib['mpole12Scale'],
                                              element.attrib['mpole13Scale'],
                                              element.attrib['mpole14Scale'],
                                              element.attrib['mpole15Scale'],

                                              element.attrib['mutual11Scale'],
                                              element.attrib['mutual12Scale'],
                                              element.attrib['mutual13Scale'],
                                              element.attrib['mutual14Scale'],

                                              element.attrib['polar12Scale'],
                                              element.attrib['polar13Scale'],
                                              element.attrib['polar14Scale'],
                                              element.attrib['polar15Scale'])



        forceField._forces.append(generator)

        # set type map: [ kIndices, multipoles, AMOEBA/OpenMM axis type]

        for atom in element.findall('Multipole'):
            types = forceField._findAtomTypes(atom, 1)
            if types is not None:

                # k-indices not provided default to 0

                kIndices = [int(atom.attrib['type'])]

                kStrings = [ 'kz', 'kx', 'ky' ]
                for kString in kStrings:
                    try:
                        if (atom.attrib[kString]):
                             kIndices.append(int(atom.attrib[kString]))
                    except: 
                        pass

                # set axis type based on k-Indices 

                axisType = AmoebaMultipoleGenerator.setAxisType(kIndices)

                # set multipole

                charge = float(atom.attrib['c0'])
 
                conversion = 1.0
                dipole = [ conversion*float(atom.attrib['d1']), conversion*float(atom.attrib['d2']), conversion*float(atom.attrib['d3'])]

                quadrupole = []
                quadrupole.append(conversion*float(atom.attrib['q11']))
                quadrupole.append(conversion*float(atom.attrib['q21']))
                quadrupole.append(conversion*float(atom.attrib['q31']))
                quadrupole.append(conversion*float(atom.attrib['q21']))
                quadrupole.append(conversion*float(atom.attrib['q22']))
                quadrupole.append(conversion*float(atom.attrib['q32']))
                quadrupole.append(conversion*float(atom.attrib['q31']))
                quadrupole.append(conversion*float(atom.attrib['q32']))
                quadrupole.append(conversion*float(atom.attrib['q33']))

                for t in types[0]:
                    if (t not in generator.typeMap):
                        generator.typeMap[t] = []

                    valueMap = dict()
                    valueMap['classIndex'] = atom.attrib['type']
                    valueMap['kIndices'] = kIndices
                    valueMap['charge'] = charge 
                    valueMap['dipole'] = dipole
                    valueMap['quadrupole'] = quadrupole
                    valueMap['axisType'] = axisType
                    generator.typeMap[t].append(valueMap)
    
            else:
                outputString = "AmoebaMultipoleGenerator: error getting type for multipole: %s" % (atom.attrib['class'])
                raise ValueError(outputString) 
    
        # polarization parameters
 
        for atom in element.findall('Polarize'):
            types = forceField._findAtomTypes(atom, 1)
            if types is not None:

                classIndex = atom.attrib['type']
                polarizability = float(atom.attrib['polarizability'])
                thole = float(atom.attrib['thole'])
                if (thole == 0):
                    pdamp = 0
                else:
                    pdamp = pow(polarizability, 1.0/6.0)

                pgrpMap = dict()
                for index in range(1, 7):
                    pgrp = 'pgrp' + str(index)
                    if (pgrp in atom.attrib):
                        pgrpMap[int(atom.attrib[pgrp])] = -1

                for t in types[0]:
                    if (t not in generator.typeMap):
                        outputString = "AmoebaMultipoleGenerator: polarize type not present: %s" % (atom.attrib['type'])
                        raise ValueError(outputString) 
                    else:
                        typeMapList = generator.typeMap[t]
                        hit = 0
                        for (ii, typeMap) in enumerate(typeMapList):

                            if (typeMap['classIndex'] == classIndex):
                                typeMap['polarizability'] = polarizability
                                typeMap['thole'] = thole
                                typeMap['pdamp'] = pdamp
                                typeMap['pgrpMap'] = pgrpMap 
                                typeMapList[ii] = typeMap
                                hit = 1

                        if (hit == 0):
                            outputString = "AmoebaMultipoleGenerator: error getting type for polarize: class index=%s not in multipole list?" % (atom.attrib['class'])
                            raise ValueError(outputString) 
    
            else:
                outputString = "AmoebaMultipoleGenerator: error getting type for polarize: %s" % (atom.attrib['class'])
                raise ValueError(outputString) 
    
    #=============================================================================================

    def setPolarGroups(self, data, bonded12ParticleSets, force):

        for (atomIndex, atom) in enumerate(data.atoms):

            # assign multipole parameters via only 1-2 connected atoms

            multipoleDict = atom.multipoleDict
            pgrpMap = multipoleDict['pgrpMap']
            bondedAtomIndices = bonded12ParticleSets[atomIndex]
            atom.stage = -1
            atom.polarizationGroupSet = list()
            atom.polarizationGroups[atomIndex] = 1
            for bondedAtomIndex in bondedAtomIndices:
                bondedAtomType = int(data.atomType[data.atoms[bondedAtomIndex]])
                bondedAtom = data.atoms[bondedAtomIndex]
                if (bondedAtomType in pgrpMap):
                    atom.polarizationGroups[bondedAtomIndex] = 1
                    bondedAtom.polarizationGroups[atomIndex] = 1
                    
        # pgrp11

        for (atomIndex, atom) in enumerate(data.atoms):

            if (len( data.atoms[atomIndex].polarizationGroupSet) > 0):
                continue

            group = set()
            visited = set()
            notVisited = set()
            for pgrpAtomIndex in atom.polarizationGroups:
                group.add(pgrpAtomIndex)
                notVisited.add(pgrpAtomIndex)
            visited.add(atomIndex)
            while(len(notVisited) > 0):
                nextAtom = notVisited.pop()
                if (nextAtom not in visited):
                   visited.add(nextAtom)
                   for ii in data.atoms[nextAtom].polarizationGroups:
                       group.add(ii)
                       if (ii not in visited):
                           notVisited.add(ii)

            pGroup = group
            for pgrpAtomIndex in group:
                data.atoms[pgrpAtomIndex].polarizationGroupSet.append(pGroup)

        for (atomIndex, atom) in enumerate(data.atoms):
            atom.polarizationGroupSet[0] = sorted(atom.polarizationGroupSet[0])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent11, atom.polarizationGroupSet[0])

        # pgrp12

        for (atomIndex, atom) in enumerate(data.atoms):

            if (len( data.atoms[atomIndex].polarizationGroupSet) > 1):
                continue

            pgrp11 = set(atom.polarizationGroupSet[0])
            pgrp12 = set()
            for pgrpAtomIndex in pgrp11:
                for bonded12 in bonded12ParticleSets[pgrpAtomIndex]:
                    pgrp12 = pgrp12.union(data.atoms[bonded12].polarizationGroupSet[0])
            pgrp12 = pgrp12 - pgrp11
            for pgrpAtomIndex in pgrp11:
                data.atoms[pgrpAtomIndex].polarizationGroupSet.append(pgrp12)
                
        for (atomIndex, atom) in enumerate(data.atoms):
            atom.polarizationGroupSet[1] = sorted(atom.polarizationGroupSet[1])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent12, atom.polarizationGroupSet[1])

        # pgrp13

        for (atomIndex, atom) in enumerate(data.atoms):

            if (len(data.atoms[atomIndex].polarizationGroupSet) > 2):
                continue

            pgrp11 = set(atom.polarizationGroupSet[0])
            pgrp12 = set(atom.polarizationGroupSet[1])
            pgrp13 = set()
            for pgrpAtomIndex in pgrp12:
                for bonded12 in bonded12ParticleSets[pgrpAtomIndex]:
                    pgrp13 = pgrp13.union(data.atoms[bonded12].polarizationGroupSet[0])
            pgrp13 = pgrp13 - pgrp12
            pgrp13 = pgrp13 - set(pgrp11)
            for pgrpAtomIndex in pgrp11:
                data.atoms[pgrpAtomIndex].polarizationGroupSet.append(pgrp13)
                
        for (atomIndex, atom) in enumerate(data.atoms):
            atom.polarizationGroupSet[2] = sorted(atom.polarizationGroupSet[2])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent13, atom.polarizationGroupSet[2])

        # pgrp14

        for (atomIndex, atom) in enumerate(data.atoms):

            if (len(data.atoms[atomIndex].polarizationGroupSet) > 3):
                continue

            pgrp11 = set(atom.polarizationGroupSet[0])
            pgrp12 = set(atom.polarizationGroupSet[1])
            pgrp13 = set(atom.polarizationGroupSet[2])
            pgrp14 = set()
            for pgrpAtomIndex in pgrp13:
                for bonded12 in bonded12ParticleSets[pgrpAtomIndex]:
                    pgrp14 = pgrp14.union(data.atoms[bonded12].polarizationGroupSet[0])

            pgrp14 = pgrp14 - pgrp13
            pgrp14 = pgrp14 - pgrp12
            pgrp14 = pgrp14 - set(pgrp11)

            for pgrpAtomIndex in pgrp11:
                data.atoms[pgrpAtomIndex].polarizationGroupSet.append(pgrp14)
                
        for (atomIndex, atom) in enumerate(data.atoms):
            atom.polarizationGroupSet[3] = sorted(atom.polarizationGroupSet[3])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent14, atom.polarizationGroupSet[3])

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        methodMap = {NoCutoff:mm.AmoebaMultipoleForce.NoCutoff,
                     PME:mm.AmoebaMultipoleForce.PME}

        # get or create force depending on whether it has already been added to the system

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaMultipoleForce]
        if len(existing) == 0:
            force = mm.AmoebaMultipoleForce()
            sys.addForce(force)
            if (nonbondedMethod not in methodMap): 
                raise ValueError( "AmoebaMultipoleForce: input cutoff method not available." )
            else:
                force.setNonbondedMethod(methodMap[nonbondedMethod])
            force.setCutoffDistance(nonbondedCutoff)

            if ('ewaldErrorTolerance' in args):
                force.setEwaldErrorTolerance(float(args['ewaldErrorTolerance']))

            if ('polarization' in args):
                polarizationType = args['polarization']
                if (polarizationType.lower() == 'direct'):
                    force.setPolarizationType(mm.AmoebaMultipoleForce.Direct)
                else:
                    force.setPolarizationType(mm.AmoebaMultipoleForce.Mutual)

            if ('aEwald' in args):
                force.setAEwald(float(args['aEwald']))

            if ('pmeGridDimensions' in args):
                force.setPmeGridDimensions(args['pmeGridDimensions'])

            if ('mutualInducedMaxIterations' in args):
                force.setMutualInducedMaxIterations(int(args['mutualInducedMaxIterations']))

            if ('mutualInducedTargetEpsilon' in args):
                force.setMutualInducedTargetEpsilon(float(args['mutualInducedTargetEpsilon']))

        else:
            force = existing[0]

        # add particles to force
        # throw error if particle type not available 

        # get 1-2, 1-3, 1-4, 1-5 bonded sets

        # 1-2

        bonded12ParticleSets = []
        for i in range(len(data.atoms)):
            bonded12ParticleSet = AmoebaVdwGenerator.getBondedParticleSet(i, data)
            bonded12ParticleSet = set(sorted(bonded12ParticleSet))
            bonded12ParticleSets.append(bonded12ParticleSet)

        # 1-3

        bonded13ParticleSets = []
        for i in range(len(data.atoms)):
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
        for i in range(len(data.atoms)):
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
        for i in range(len(data.atoms)):
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

        for (atomIndex, atom) in enumerate(data.atoms):
            t = data.atomType[atom]
            if t in self.typeMap:

                multipoleList = self.typeMap[t]
                hit = 0
                savedMultipoleDict = 0

                # assign multipole parameters via only 1-2 connected atoms

                for multipoleDict in multipoleList:

                    if (hit != 0):
                        break

                    kIndices = multipoleDict['kIndices']
    
                    kz = kIndices[1]   
                    kx = kIndices[2]
                    ky = kIndices[3]

                    # assign multipole parameters
                    #    (1) get bonded partners
                    #    (2) match parameter types
    
                    bondedAtomIndices = bonded12ParticleSets[atomIndex]
                    zaxis = -1
                    xaxis = -1
                    yaxis = -1
                    for bondedAtomZIndex in bondedAtomIndices:

                       if (hit != 0):
                           break

                       bondedAtomZType = int(data.atomType[data.atoms[bondedAtomZIndex]])
                       bondedAtomZ = data.atoms[bondedAtomZIndex]
                       if (bondedAtomZType == kz):
                          for bondedAtomXIndex in bondedAtomIndices:
                              if (bondedAtomXIndex == bondedAtomZIndex or hit != 0):
                                  continue
                              bondedAtomXType = int(data.atomType[data.atoms[bondedAtomXIndex]])
                              if (bondedAtomXType == kx):
                                  if (ky == 0):
                                      zaxis = bondedAtomZIndex
                                      xaxis = bondedAtomXIndex
                                      if( bondedAtomXType == bondedAtomZType and xaxis < zaxis ):
                                          swapI = zaxis
                                          zaxis = xaxis
                                          xaxis = swapI
                                      else:
                                          for bondedAtomXIndex in bondedAtomIndices:
                                              bondedAtomX1Type = int(data.atomType[data.atoms[bondedAtomXIndex]])
                                              if( bondedAtomX1Type == kx and bondedAtomXIndex != bondedAtomZIndex and bondedAtomXIndex < xaxis ):
                                                  xaxis = bondedAtomXIndex

                                      savedMultipoleDict = multipoleDict
                                      hit = 1
                                  else:
                                      for bondedAtomYIndex in bondedAtomIndices:
                                          if (bondedAtomYIndex == bondedAtomZIndex or bondedAtomYIndex == bondedAtomXIndex or hit != 0):
                                              continue
                                          bondedAtomYType = int(data.atomType[data.atoms[bondedAtomYIndex]])
                                          if (bondedAtomYType == ky):
                                              zaxis = bondedAtomZIndex
                                              xaxis = bondedAtomXIndex
                                              yaxis = bondedAtomYIndex
                                              savedMultipoleDict = multipoleDict
                                              hit = 2
                                         
                # assign multipole parameters via 1-2 and 1-3 connected atoms

                for multipoleDict in multipoleList:

                    if (hit != 0):
                        break

                    kIndices = multipoleDict['kIndices']
    
                    kz = kIndices[1]   
                    kx = kIndices[2]
                    ky = kIndices[3]
    
                    # assign multipole parameters
                    #    (1) get bonded partners
                    #    (2) match parameter types
    
                    bondedAtom12Indices = bonded12ParticleSets[atomIndex]
                    bondedAtom13Indices = bonded13ParticleSets[atomIndex]

                    zaxis = -1
                    xaxis = -1
                    yaxis = -1

                    for bondedAtomZIndex in bondedAtom12Indices:

                       if (hit != 0):
                           break

                       bondedAtomZType = int(data.atomType[data.atoms[bondedAtomZIndex]])
                       bondedAtomZ = data.atoms[bondedAtomZIndex]

                       if (bondedAtomZType == kz):
                          for bondedAtomXIndex in bondedAtom13Indices:

                              if (bondedAtomXIndex == bondedAtomZIndex or hit != 0):
                                  continue
                              bondedAtomXType = int(data.atomType[data.atoms[bondedAtomXIndex]])
                              if (bondedAtomXType == kx and bondedAtomZIndex in bonded12ParticleSets[bondedAtomXIndex]):
                                  if (ky == 0):
                                      zaxis = bondedAtomZIndex
                                      xaxis = bondedAtomXIndex

                                      # select xaxis w/ smallest index

                                      for bondedAtomXIndex in bondedAtom13Indices:
                                          bondedAtomX1Type = int(data.atomType[data.atoms[bondedAtomXIndex]])
                                          if( bondedAtomX1Type == kx and bondedAtomXIndex != bondedAtomZIndex and bondedAtomZIndex in bonded12ParticleSets[bondedAtomXIndex] and bondedAtomXIndex < xaxis ):
                                              xaxis = bondedAtomXIndex

                                      savedMultipoleDict = multipoleDict
                                      hit = 3
                                  else:
                                      for bondedAtomYIndex in bondedAtom13Indices:
                                          if (bondedAtomYIndex == bondedAtomZIndex or bondedAtomYIndex == bondedAtomXIndex or hit != 0):
                                              continue
                                          bondedAtomYType = int(data.atomType[data.atoms[bondedAtomYIndex]])
                                          if (bondedAtomYType == ky and bondedAtomZIndex in bonded12ParticleSets[bondedAtomYIndex]):
                                              zaxis = bondedAtomZIndex
                                              xaxis = bondedAtomXIndex
                                              yaxis = bondedAtomYIndex
                                              savedMultipoleDict = multipoleDict
                                              hit = 4
                                         
                # assign multipole parameters via only a z-defining atom

                for multipoleDict in multipoleList:

                    if (hit != 0):
                        break

                    kIndices = multipoleDict['kIndices']
    
                    kz = kIndices[1]   
                    kx = kIndices[2]   
    
                    zaxis = -1
                    xaxis = -1
                    yaxis = -1

                    for bondedAtomZIndex in bondedAtom12Indices:

                        if (hit != 0):
                            break

                        bondedAtomZType = int(data.atomType[data.atoms[bondedAtomZIndex]])
                        bondedAtomZ = data.atoms[bondedAtomZIndex]

                        if (kx == 0 and kz == bondedAtomZType):
                            kz = bondedAtomZIndex
                            savedMultipoleDict = multipoleDict
                            hit = 5

                # assign multipole parameters via no connected atoms

                for multipoleDict in multipoleList:

                    if (hit != 0):
                        break

                    kIndices = multipoleDict['kIndices']
    
                    kz = kIndices[1]   
    
                    zaxis = -1
                    xaxis = -1
                    yaxis = -1

                    if (kz == 0):
                        savedMultipoleDict = multipoleDict
                        hit = 6
                
                # add particle if there was a hit

                if (hit != 0):

                    atom.multipoleDict = savedMultipoleDict
                    atom.polarizationGroups = dict()
                    newIndex = force.addMultipole(savedMultipoleDict['charge'], savedMultipoleDict['dipole'], savedMultipoleDict['quadrupole'], savedMultipoleDict['axisType'],
                                                                 zaxis, xaxis, yaxis, savedMultipoleDict['thole'], savedMultipoleDict['pdamp'], savedMultipoleDict['polarizability'])
                    if (atomIndex == newIndex):
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent12, bonded12ParticleSets[atomIndex])
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent13, bonded13ParticleSets[atomIndex])
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent14, bonded14ParticleSets[atomIndex])
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent15, bonded15ParticleSets[atomIndex])
                    else:
                        raise ValueError("Atom %s of %s %d is out of synch!." %(atom.name, atom.residue.name, atom.residue.index))
                else:
                    raise ValueError("Atom %s of %s %d was not assigned." %(atom.name, atom.residue.name, atom.residue.index))
            else:
                raise ValueError('No multipole type for atom %s %s %d' % (atom.name, atom.residue.name, atom.residue.index))

        # set polar groups

        self.setPolarGroups(data, bonded12ParticleSets, force)

parsers["AmoebaMultipoleForce"] = AmoebaMultipoleGenerator.parseElement

#=============================================================================================

## @private
class AmoebaWcaDispersionGenerator:

    """A AmoebaWcaDispersionGenerator constructs a AmoebaWcaDispersionForce."""
    
    #=========================================================================================

    def __init__(self, epso, epsh, rmino, rminh, awater, slevy, dispoff, shctd):

        self.epso = epso 
        self.epsh = epsh 
        self.rmino = rmino
        self.rminh = rminh
        self.awater = awater
        self.slevy = slevy
        self.dispoff = dispoff
        self.shctd = shctd 

        self.typeMap = {}

    #=========================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaWcaDispersionForce epso="0.46024" epsh="0.056484" rmino="0.17025" rminh="0.13275" awater="33.428" slevy="1.0"  dispoff="0.026" shctd="0.81" >
        #   <WcaDispersion class="1" radius="0.1855" epsilon="0.46024" />
        #   <WcaDispersion class="2" radius="0.191" epsilon="0.422584" />
      
        generator = AmoebaWcaDispersionGenerator(element.attrib['epso'],
                                                  element.attrib['epsh'],
                                                  element.attrib['rmino'],
                                                  element.attrib['rminh'],
                                                  element.attrib['awater'], 
                                                  element.attrib['slevy'],
                                                  element.attrib['dispoff'],
                                                  element.attrib['shctd']) 
        forceField._forces.append(generator)

        # typeMap[] = [ radius, epsilon ]

        for atom in element.findall('WcaDispersion'):
            types = forceField._findAtomTypes(atom, 1)
            if types is not None:

                values = [float(atom.attrib['radius']), float(atom.attrib['epsilon'])]
                for t in types[0]:
                    generator.typeMap[t] = values
            else:
                outputString = "AmoebaWcaDispersionGenerator: error getting type: %s" % (atom.attrib['class'])
                raise ValueError(outputString) 
    
    #=========================================================================================
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        # get or create force depending on whether it has already been added to the system

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaWcaDispersionForce]
        if len(existing) == 0:
            force = mm.AmoebaWcaDispersionForce()
            sys.addForce(force)
        else:
            force = existing[0]

        # add particles to force
        # throw error if particle type not available 

        force.setEpso(   float(self.epso   ))
        force.setEpsh(   float(self.epsh   ))
        force.setRmino(  float(self.rmino  ))
        force.setRminh(  float(self.rminh  ))
        force.setDispoff(float(self.dispoff))
        force.setSlevy(  float(self.slevy  ))
        force.setAwater( float(self.awater ))
        force.setShctd(  float(self.shctd  ))

        for (i, atom) in enumerate(data.atoms):
            t = data.atomType[atom]
            if t in self.typeMap:

                values = self.typeMap[t]
                force.addParticle(values[0], values[1])
            else:
                raise ValueError('No WcaDispersion type for atom %s of %s %d' % (atom.name, atom.residue.name, atom.residue.index))

parsers["AmoebaWcaDispersionForce"] = AmoebaWcaDispersionGenerator.parseElement

#=============================================================================================

## @private
class AmoebaGeneralizedKirkwoodGenerator:

    """A AmoebaGeneralizedKirkwoodGenerator constructs a AmoebaGeneralizedKirkwoodForce."""
    
    #=========================================================================================

    def __init__(self, forceField, solventDielectric, soluteDielectric, includeCavityTerm, probeRadius, surfaceAreaFactor):

        self.forceField = forceField
        self.solventDielectric = solventDielectric
        self.soluteDielectric = soluteDielectric
        self.includeCavityTerm = includeCavityTerm
        self.probeRadius = probeRadius
        self.surfaceAreaFactor = surfaceAreaFactor

        self.radiusTypeMap = {}
        self.radiusTypeMap['Bondi'] = {}
        bondiMap = self.radiusTypeMap['Bondi'] 
        rscale = 1.03

        bondiMap[0] = 0.00
        bondiMap[1] = 0.12*rscale
        bondiMap[2] = 0.14*rscale
        bondiMap[5] = 0.18*rscale

        bondiMap[6] = 0.170*rscale
        bondiMap[7] = 0.155*rscale
        bondiMap[8] = 0.152*rscale
        bondiMap[9] = 0.147*rscale

        bondiMap[10] = 0.154*rscale
        bondiMap[14] = 0.210*rscale
        bondiMap[15] = 0.180*rscale
        bondiMap[16] = 0.180*rscale

        bondiMap[17] = 0.175 *rscale
        bondiMap[18] = 0.188*rscale
        bondiMap[34] = 0.190*rscale
        bondiMap[35] = 0.185*rscale

        bondiMap[36] = 0.202*rscale
        bondiMap[53] = 0.198*rscale
        bondiMap[54] = 0.216*rscale

    #=========================================================================================

    def getObcShct(self, data, atomIndex):

        atom = data.atoms[atomIndex]
        atomicNumber = atom.element.atomic_number
        shct = -1.0

        # shct
 
        if (atomicNumber == 1):                 # H(1)
            shct = 0.85       
        elif (atomicNumber == 6):               # C(6)
            shct = 0.72         
        elif (atomicNumber == 7):               # N(7)
            shct = 0.79        
        elif (atomicNumber == 8):               # O(8)
            shct = 0.85       
        elif (atomicNumber == 9):               # F(9)
            shct = 0.88    
        elif (atomicNumber == 15):              # P(15)              
            shct = 0.86 
        elif (atomicNumber == 16):              # S(16)
            shct = 0.96
        elif (atomicNumber == 26):              # Fe(26)
            shct = 0.88

        if (shct < 0.0): 
            raise ValueError( "getObcShct: no GK overlap scale factor for atom %s of %s %d" % (atom.name, atom.residue.name, atom.residue.index) )
 
        return shct 

    #=========================================================================================

    def getAmoebaTypeRadius(self, data, bondedAtomIndices, atomIndex):

        atom = data.atoms[atomIndex]
        atomicNumber = atom.element.atomic_number
        radius = -1.0

        if (atomicNumber == 1):                  # H(1)
 
            radius = 0.132
 
            if (len(bondedAtomIndices) < 1):
                 raise ValueError( "AmoebaGeneralizedKirkwoodGenerator: error getting atom bonded to %s of %s %d " % (atom.name, atom.residue.name, atom.residue.index) )
 
            for bondedAtomIndex in bondedAtomIndices:
                bondedAtomAtomicNumber = data.atoms[bondedAtomIndex].element.atomic_number

            if (bondedAtomAtomicNumber == 7):
                radius = 0.11
            if (bondedAtomAtomicNumber == 8):
                radius = 0.105
 
        elif (atomicNumber == 3):               # Li(3)
            radius = 0.15
        elif (atomicNumber == 6):               # C(6)
            
            radius = 0.20
            if (len(bondedAtomIndices) == 3):
                radius = 0.205

            elif (len(bondedAtomIndices) == 4):
                for bondedAtomIndex in bondedAtomIndices:
                   bondedAtomAtomicNumber = data.atoms[bondedAtomIndex].element.atomic_number
                   if (bondedAtomAtomicNumber == 7 or bondedAtomAtomicNumber == 8):
                       radius = 0.175

        elif (atomicNumber == 7):               # N(7)
            radius = 0.16
        elif (atomicNumber == 8):               # O(8)
            radius = 0.155
            if (len(bondedAtomIndices) == 2):
                radius = 0.145
        elif (atomicNumber == 9):               # F(9)
            radius = 0.154
        elif (atomicNumber == 10):              
            radius = 0.146
        elif (atomicNumber == 11):              
            radius = 0.209
        elif (atomicNumber == 12):              
            radius = 0.179
        elif (atomicNumber == 14):              
            radius = 0.189
        elif (atomicNumber == 15):              # P(15)              
            radius = 0.196
        elif (atomicNumber == 16):              # S(16)
            radius = 0.186
        elif (atomicNumber == 17):              
            radius = 0.182
        elif (atomicNumber == 18):              
            radius = 0.179
        elif (atomicNumber == 19):              
            radius = 0.223
        elif (atomicNumber == 20):              
            radius = 0.191
        elif (atomicNumber == 35):         
            radius = 2.00
        elif (atomicNumber == 36):   
            radius = 0.190
        elif (atomicNumber == 37):              
            radius = 0.226
        elif (atomicNumber == 53):              
            radius = 0.237
        elif (atomicNumber == 54):              
            radius = 0.207
        elif (atomicNumber == 55):              
            radius = 0.263
        elif (atomicNumber == 56):         
            radius = 0.230

        if (radius < 0.0): 
            outputString = "No GK radius for atom %s of %s %d" % (atom.name, atom.residue.name, atom.residue.index)
            raise ValueError( outputString )
 
        return radius

    #=========================================================================================

    def getBondiTypeRadius(self, data, bondedAtomIndices, atomIndex):

        bondiMap = self.radiusTypeMap['Bondi'] 
        atom = data.atoms[atomIndex]
        atomicNumber = atom.element.atomic_number
        if (atomicNumber in bondiMap): 
            radius = bondiMap[atomicNumber]
        else:
            outputString = "Warning no Bondi radius for atom %s of %s %d using default value=%f" % (atom.name, atom.residue.name, atom.residue.index, radius)
            raise ValueError( outputString )
 
        return radius

    #=========================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaGeneralizedKirkwoodForce solventDielectric="78.3" soluteDielectric="1.0" includeCavityTerm="1" probeRadius="0.14" surfaceAreaFactor="-170.351730663">
        #   <GeneralizedKirkwood type="1" charge="-0.22620" shct="0.79"  />
        #   <GeneralizedKirkwood type="2" charge="-0.15245" shct="0.72"  />
        
        generator = AmoebaGeneralizedKirkwoodGenerator(forceField, element.attrib['solventDielectric'], element.attrib['soluteDielectric'],
                                                        element.attrib['includeCavityTerm'], 
                                                        element.attrib['probeRadius'], element.attrib['surfaceAreaFactor']) 
        forceField._forces.append(generator)

    #=========================================================================================
    
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        if( nonbondedMethod != NoCutoff ):
            raise ValueError( "Only the nonbondedMethod=NoCutoff option is available for implicit solvent simulations." )

        # check if AmoebaMultipoleForce exists since charges needed
        # if it has not been created, raise an error

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        amoebaMultipoleForceList = [f for f in existing if type(f) == mm.AmoebaMultipoleForce]
        if (len(amoebaMultipoleForceList) > 0):
            amoebaMultipoleForce = amoebaMultipoleForceList[0]
        else:
            # call AmoebaMultipoleForceGenerator.createForce() to ensure charges have been set

            for force in self.forceField._forces:
                if (force.__class__.__name__ == 'AmoebaMultipoleGenerator'): 
                    force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
          
        # get or create force depending on whether it has already been added to the system

        existing = [f for f in existing if type(f) == mm.AmoebaGeneralizedKirkwoodForce]
        if len(existing) == 0:

            force = mm.AmoebaGeneralizedKirkwoodForce()
            sys.addForce(force)
 
            if ('solventDielectric' in args):
                force.setSolventDielectric(float(args['solventDielectric']))
            else:
                force.setSolventDielectric(   float(self.solventDielectric))

            if ('soluteDielectric' in args):
                force.setSoluteDielectric(float(args['soluteDielectric']))
            else:
                force.setSoluteDielectric(    float(self.soluteDielectric))

            if ('includeCavityTerm' in args):
                force.setIncludeCavityTerm(int(args['includeCavityTerm']))
            else:
               force.setIncludeCavityTerm(   int(self.includeCavityTerm))

        else:
            force = existing[0]

        # add particles to force
        # throw error if particle type not available 

        force.setProbeRadius(         float(self.probeRadius))
        force.setSurfaceAreaFactor(   float(self.surfaceAreaFactor))

        # 1-2

        bonded12ParticleSets = []
        for i in range(len(data.atoms)):
            bonded12ParticleSet = AmoebaVdwGenerator.getBondedParticleSet(i, data)
            bonded12ParticleSet = set(sorted(bonded12ParticleSet))
            bonded12ParticleSets.append(bonded12ParticleSet)

        radiusType = 'Bondi'
        for atomIndex in range(0, amoebaMultipoleForce.getNumMultipoles()):
            multipoleParameters = amoebaMultipoleForce.getMultipoleParameters(atomIndex)
            if (radiusType == 'Amoeba'):
                radius = self.getAmoebaTypeRadius(data, bonded12ParticleSets[atomIndex], atomIndex)
            else:
                radius = self.getBondiTypeRadius(data, bonded12ParticleSets[atomIndex], atomIndex)
            #shct = self.getObcShct(data, atomIndex)
            shct = 0.69
            force.addParticle(multipoleParameters[0], radius, shct)

parsers["AmoebaGeneralizedKirkwoodForce"] = AmoebaGeneralizedKirkwoodGenerator.parseElement

#=============================================================================================

## @private
class AmoebaUreyBradleyGenerator:

    #=============================================================================================
    """An AmoebaUreyBradleyGenerator constructs a AmoebaUreyBradleyForce."""
    #=============================================================================================
    
    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

        self.length = []
        self.k = []

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaUreyBradleyForce>
        #   <UreyBradley class1="74" class2="73" class3="74" k="16003.8" d="0.15537" /> 

        generator = AmoebaUreyBradleyGenerator()
        forceField._forces.append(generator)
        for bond in element.findall('UreyBradley'):
            types = forceField._findAtomTypes(bond, 3)
            if types is not None:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

                generator.length.append(float(bond.attrib['d']))
                generator.k.append(float(bond.attrib['k']))

            else:
                outputString = "AmoebaUreyBradleyGenerator: error getting types: %s %s %s" % (
                                    bond.attrib['class1'], bond.attrib['class2'], bond.attrib['class3'])
                raise ValueError(outputString) 
    
    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.HarmonicBondForce]

        if len(existing) == 0:
            force = mm.HarmonicBondForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):
            if (isConstrained):
                continue
            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                if ((type1 in types1 and type2 in types2 and type3 in types3) or (type3 in types1 and type2 in types2 and type1 in types3)):
                    force.addBond(angle[0], angle[2], self.length[i], 2*self.k[i])
                    break

parsers["AmoebaUreyBradleyForce"] = AmoebaUreyBradleyGenerator.parseElement

#=============================================================================================
