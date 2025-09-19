"""
forcefield.py: Constructs OpenMM System objects based on a Topology and an XML force field description

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2025 Stanford University and the Authors.
Authors: Peter Eastman, Mark Friedrichs
Contributors: Evan Pretti

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
from __future__ import absolute_import, print_function

__author__ = "Peter Eastman"
__version__ = "1.0"

import os
import itertools
import xml.etree.ElementTree as etree
import math
import warnings
from math import sqrt, cos
from copy import deepcopy
from collections import Counter, defaultdict
from difflib import SequenceMatcher
import openmm as mm
import openmm.unit as unit
from . import element as elem
from openmm.app.internal.singleton import Singleton
from openmm.app.internal import compiled, amoebaforces
from openmm.app.internal.argtracker import ArgTracker

# Directories from which to load built in force fields.

_dataDirectories = None

def _getDataDirectories():
    global _dataDirectories
    if _dataDirectories is None:
        _dataDirectories = [os.path.join(os.path.dirname(__file__), 'data')]
        try:
            from importlib_metadata import entry_points
        except:
            try:
                from importlib.metadata import entry_points
            except:
                pass
        try:
            for entry in entry_points().select(group='openmm.forcefielddir'):
                _dataDirectories.append(entry.load()())
        except:
            pass # importlib_metadata is not installed
    return _dataDirectories

def _convertParameterToNumber(param):
    if unit.is_quantity(param):
        if param.unit.is_compatible(unit.bar):
            return param / unit.bar
        return param.value_in_unit_system(unit.md_unit_system)
    return float(param)

def _parseFunctions(element):
    """Parse the attributes on an XML tag to find any tabulated functions it defines."""
    functions = []
    for function in element.findall('Function'):
        values = [float(x) for x in function.text.split()]
        if 'type' in function.attrib:
            functionType = function.attrib['type']
        else:
            functionType = 'Continuous1D'
        params = {}
        for key in function.attrib:
            if key.endswith('size'):
                params[key] = int(function.attrib[key])
            elif key.endswith('min') or key.endswith('max'):
                params[key] = float(function.attrib[key])
        if functionType.startswith('Continuous'):
            periodicStr = function.attrib.get('periodic', 'false').lower()
            if periodicStr in ['true', 'false', 'yes', 'no', '1', '0']:
                params['periodic'] = periodicStr in ['true', 'yes', '1']
            else:
                raise ValueError('ForceField: non-boolean value for periodic attribute in tabulated function definition')
        functions.append((function.attrib['name'], functionType, values, params))
    return functions

def _createFunctions(force, functions):
    """Add TabulatedFunctions to a Force based on the information that was recorded by _parseFunctions()."""
    for (name, type, values, params) in functions:
        if type == 'Continuous1D':
            force.addTabulatedFunction(
                name,
                mm.Continuous1DFunction(values, params['min'], params['max'], params['periodic']),
            )
        elif type == 'Continuous2D':
            force.addTabulatedFunction(
                name,
                mm.Continuous2DFunction(
                    params['xsize'], params['ysize'],
                    values,
                    params['xmin'], params['xmax'],
                    params['ymin'], params['ymax'],
                    params['periodic'],
                ),
            )
        elif type == 'Continuous3D':
            force.addTabulatedFunction(
                name,
                mm.Continuous3DFunction(
                    params['xsize'], params['ysize'], params['zsize'],
                    values,
                    params['xmin'], params['xmax'],
                    params['ymin'], params['ymax'],
                    params['zmin'], params['zmax'],
                    params['periodic'],
                ),
            )
        elif type == 'Discrete1D':
            force.addTabulatedFunction(name, mm.Discrete1DFunction(values))
        elif type == 'Discrete2D':
            force.addTabulatedFunction(name, mm.Discrete2DFunction(params['xsize'], params['ysize'], values))
        elif type == 'Discrete3D':
            force.addTabulatedFunction(name, mm.Discrete3DFunction(params['xsize'], params['ysize'], params['zsize'], values))

# Enumerated values for nonbonded method

class NoCutoff(Singleton):
    def __repr__(self):
        return 'NoCutoff'
NoCutoff = NoCutoff()

class CutoffNonPeriodic(Singleton):
    def __repr__(self):
        return 'CutoffNonPeriodic'
CutoffNonPeriodic = CutoffNonPeriodic()

class CutoffPeriodic(Singleton):
    def __repr__(self):
        return 'CutoffPeriodic'
CutoffPeriodic = CutoffPeriodic()

class Ewald(Singleton):
    def __repr__(self):
        return 'Ewald'
Ewald = Ewald()

class PME(Singleton):
    def __repr__(self):
        return 'PME'
PME = PME()

class LJPME(Singleton):
    def __repr__(self):
        return 'LJPME'
LJPME = LJPME()

# Enumerated values for constraint type

class HBonds(Singleton):
    def __repr__(self):
        return 'HBonds'
HBonds = HBonds()

class AllBonds(Singleton):
    def __repr__(self):
        return 'AllBonds'
AllBonds = AllBonds()

class HAngles(Singleton):
    def __repr__(self):
        return 'HAngles'
HAngles = HAngles()

# A map of functions to parse elements of the XML file.

parsers = {}

class ForceField(object):
    """A ForceField constructs OpenMM System objects based on a Topology."""

    def __init__(self, *files):
        """Load one or more XML files and create a ForceField object based on them.

        Parameters
        ----------
        files : list
            A list of XML files defining the force field.  Each entry may
            be an absolute file path, a path relative to the current working
            directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a
            read() method from which the forcefield XML data can be loaded.
        """
        self._atomTypes = {}
        self._templates = {}
        self._patches = {}
        self._templatePatches = {}
        self._templateSignatures = {None:[]}
        self._atomClasses = {'':set()}
        self._forces = []
        self._scripts = []
        self._templateMatchers = []
        self._templateGenerators = []
        self.loadFile(files)

    def loadFile(self, files, resname_prefix=''):
        """Load an XML file and add the definitions from it to this ForceField.

        Parameters
        ----------
        files : string or file or tuple
            An XML file or tuple of XML files containing force field definitions.
            Each entry may be either an absolute file path, a path relative to the current working
            directory, a path relative to this module's data subdirectory (for
            built in force fields), or an open file-like object with a read()
            method from which the forcefield XML data can be loaded.
        prefix : string
            An optional string to be prepended to each residue name found in the
            loaded files.
        """

        if isinstance(files, tuple):
            files = list(files)
        else:
            files = [files]

        trees = []

        i = 0
        while i < len(files):
            file = files[i]
            # this handles either filenames or open file-like objects
            if isinstance(file, str) and not os.path.isfile(file):
                for dataDir in _getDataDirectories():
                    f = os.path.join(dataDir, file)
                    if os.path.isfile(f):
                        file = f
                        break
            try:
                tree = etree.parse(file)
            except FileNotFoundError:
                raise ValueError('Could not locate file "%s"' % file)
            except Exception as e:
                # Fail with an error message about which file could not be read.
                if hasattr(file, 'name'):
                    filename = file.name
                else:
                    filename = str(file)
                raise Exception('ForceField.loadFile() encountered an error reading file "%s": %s' % (filename, e))

            trees.append(tree)
            i += 1

            # Process includes in this file.

            if isinstance(file, str):
                parentDir = os.path.dirname(file)
            else:
                parentDir = ''
            for included in tree.getroot().findall('Include'):
                includeFile = included.attrib['file']
                joined = os.path.join(parentDir, includeFile)
                if os.path.isfile(joined):
                    includeFile = joined
                if includeFile not in files:
                    files.append(includeFile)

        # Load the atom types.

        for tree in trees:
            if tree.getroot().find('AtomTypes') is not None:
                for type in tree.getroot().find('AtomTypes').findall('Type'):
                    self.registerAtomType(type.attrib)

        # Load the residue templates.

        for tree in trees:
            if tree.getroot().find('Residues') is not None:
                for residue in tree.getroot().find('Residues').findall('Residue'):
                    resName = resname_prefix+residue.attrib['name']
                    template = ForceField._TemplateData(resName)
                    if 'override' in residue.attrib:
                        template.overrideLevel = int(residue.attrib['override'])
                    if 'rigidWater' in residue.attrib:
                        template.rigidWater = (residue.attrib['rigidWater'].lower() == 'true')
                    for key in residue.attrib:
                        template.attributes[key] = residue.attrib[key]
                    atomIndices = template.atomIndices
                    for ia, atom in enumerate(residue.findall('Atom')):
                        params = {}
                        for key in atom.attrib:
                            if key not in ('name', 'type'):
                                params[key] = _convertParameterToNumber(atom.attrib[key])
                        atomName = atom.attrib['name']
                        if atomName in atomIndices:
                            raise ValueError('Residue '+resName+' contains multiple atoms named '+atomName)
                        typeName = atom.attrib['type']
                        atomIndices[atomName] = ia
                        template.atoms.append(ForceField._TemplateAtomData(atomName, typeName, self._atomTypes[typeName].element, params))
                    for site in residue.findall('VirtualSite'):
                        template.virtualSites.append(ForceField._VirtualSiteData(site, atomIndices))
                    for bond in residue.findall('Bond'):
                        if 'atomName1' in bond.attrib:
                            template.addBondByName(bond.attrib['atomName1'], bond.attrib['atomName2'])
                        else:
                            template.addBond(int(bond.attrib['from']), int(bond.attrib['to']))
                    for bond in residue.findall('ExternalBond'):
                        if 'atomName' in bond.attrib:
                            template.addExternalBondByName(bond.attrib['atomName'])
                        else:
                            template.addExternalBond(int(bond.attrib['from']))
                    for patch in residue.findall('AllowPatch'):
                        patchName = patch.attrib['name']
                        if ':' in patchName:
                            colonIndex = patchName.find(':')
                            self.registerTemplatePatch(resName, patchName[:colonIndex], int(patchName[colonIndex+1:])-1)
                        else:
                            self.registerTemplatePatch(resName, patchName, 0)
                    self.registerResidueTemplate(template)

        # Load the patch definitions.

        for tree in trees:
            if tree.getroot().find('Patches') is not None:
                for patch in tree.getroot().find('Patches').findall('Patch'):
                    patchName = patch.attrib['name']
                    if 'residues' in patch.attrib:
                        numResidues = int(patch.attrib['residues'])
                    else:
                        numResidues = 1
                    patchData = ForceField._PatchData(patchName, numResidues)
                    for key in patch.attrib:
                        patchData.attributes[key] = patch.attrib[key]
                    for atom in patch.findall('AddAtom'):
                        params = {}
                        for key in atom.attrib:
                            if key not in ('name', 'type'):
                                params[key] = _convertParameterToNumber(atom.attrib[key])
                        atomName = atom.attrib['name']
                        if atomName in patchData.allAtomNames:
                            raise ValueError('Patch '+patchName+' contains multiple atoms named '+atomName)
                        patchData.allAtomNames.add(atomName)
                        atomDescription = ForceField._PatchAtomData(atomName)
                        typeName = atom.attrib['type']
                        patchData.addedAtoms[atomDescription.residue].append(ForceField._TemplateAtomData(atomDescription.name, typeName, self._atomTypes[typeName].element, params))
                    for atom in patch.findall('ChangeAtom'):
                        params = {}
                        for key in atom.attrib:
                            if key not in ('name', 'type'):
                                params[key] = _convertParameterToNumber(atom.attrib[key])
                        atomName = atom.attrib['name']
                        if atomName in patchData.allAtomNames:
                            raise ValueError('Patch '+patchName+' contains multiple atoms named '+atomName)
                        patchData.allAtomNames.add(atomName)
                        atomDescription = ForceField._PatchAtomData(atomName)
                        typeName = atom.attrib['type']
                        patchData.changedAtoms[atomDescription.residue].append(ForceField._TemplateAtomData(atomDescription.name, typeName, self._atomTypes[typeName].element, params))
                    for atom in patch.findall('RemoveAtom'):
                        atomName = atom.attrib['name']
                        if atomName in patchData.allAtomNames:
                            raise ValueError('Patch '+patchName+' contains multiple atoms named '+atomName)
                        patchData.allAtomNames.add(atomName)
                        atomDescription = ForceField._PatchAtomData(atomName)
                        patchData.deletedAtoms.append(atomDescription)
                    for bond in patch.findall('AddBond'):
                        atom1 = ForceField._PatchAtomData(bond.attrib['atomName1'])
                        atom2 = ForceField._PatchAtomData(bond.attrib['atomName2'])
                        patchData.addedBonds.append((atom1, atom2))
                    for bond in patch.findall('RemoveBond'):
                        atom1 = ForceField._PatchAtomData(bond.attrib['atomName1'])
                        atom2 = ForceField._PatchAtomData(bond.attrib['atomName2'])
                        patchData.deletedBonds.append((atom1, atom2))
                    for bond in patch.findall('AddExternalBond'):
                        atom = ForceField._PatchAtomData(bond.attrib['atomName'])
                        patchData.addedExternalBonds.append(atom)
                    for bond in patch.findall('RemoveExternalBond'):
                        atom = ForceField._PatchAtomData(bond.attrib['atomName'])
                        patchData.deletedExternalBonds.append(atom)
                    # The following three lines are only correct for single residue patches.  Multi-residue patches with
                    # virtual sites currently don't work correctly.  See issue #2848.
                    atomIndices = dict((atom.name, i) for i, atom in enumerate(patchData.addedAtoms[0]+patchData.changedAtoms[0]))
                    for site in patch.findall('VirtualSite'):
                        patchData.virtualSites[0].append(ForceField._VirtualSiteData(site, atomIndices))
                    for residue in patch.findall('ApplyToResidue'):
                        name = residue.attrib['name']
                        if ':' in name:
                            colonIndex = name.find(':')
                            self.registerTemplatePatch(name[colonIndex+1:], patchName, int(name[:colonIndex])-1)
                        else:
                            self.registerTemplatePatch(name, patchName, 0)
                    self.registerPatch(patchData)

        # Load force definitions

        for tree in trees:
            for child in tree.getroot():
                if child.tag in parsers:
                    parsers[child.tag](child, self)

        # Load scripts

        for tree in trees:
            for node in tree.getroot().findall('Script'):
                self.registerScript(node.text)

        # Execute initialization scripts.

        for tree in trees:
            for node in tree.getroot().findall('InitializationScript'):
                exec(node.text, locals())

    def getGenerators(self):
        """Get the list of all registered generators."""
        return self._forces

    def registerGenerator(self, generator):
        """Register a new generator."""
        self._forces.append(generator)

    def registerAtomType(self, parameters):
        """Register a new atom type."""
        name = parameters['name']
        if name in self._atomTypes:
            #  allow multiple registrations of the same atom type provided the definitions are identical
            existing = self._atomTypes[name]
            elementsMatch = ((existing.element is None and 'element' not in parameters) or (existing.element is not None and 'element' in parameters and existing.element.symbol == parameters['element']))
            if existing.atomClass == parameters['class'] and existing.mass == float(parameters['mass']) and elementsMatch:
                return
            raise ValueError('Found multiple definitions for atom type: '+name)
        atomClass = parameters['class']
        mass = _convertParameterToNumber(parameters['mass'])
        element = None
        if 'element' in parameters:
            element = parameters['element']
            if not isinstance(element, elem.Element):
                element = elem.get_by_symbol(element)
        self._atomTypes[name] = ForceField._AtomType(name, atomClass, mass, element)
        if atomClass in self._atomClasses:
            typeSet = self._atomClasses[atomClass]
        else:
            typeSet = set()
            self._atomClasses[atomClass] = typeSet
        typeSet.add(name)
        self._atomClasses[''].add(name)

    def registerResidueTemplate(self, template):
        """Register a new residue template."""
        if template.name in self._templates:
            # There is already a template with this name, so check the override levels.

            existingTemplate = self._templates[template.name]
            if template.overrideLevel < existingTemplate.overrideLevel:
                # The existing one takes precedence, so just return.
                return
            if template.overrideLevel > existingTemplate.overrideLevel:
                # We need to delete the existing template.
                del self._templates[template.name]
                existingSignature = _createResidueSignature([atom.element for atom in existingTemplate.atoms])
                self._templateSignatures[existingSignature].remove(existingTemplate)
            else:
                raise ValueError('Residue template %s with the same override level %d already exists.' % (template.name, template.overrideLevel))

        # Register the template.

        self._templates[template.name] = template
        signature = _createResidueSignature([atom.element for atom in template.atoms])
        if signature in self._templateSignatures:
            self._templateSignatures[signature].append(template)
        else:
            self._templateSignatures[signature] = [template]

    def registerPatch(self, patch):
        """Register a new patch that can be applied to templates."""
        patch.index = len(self._patches)
        self._patches[patch.name] = patch

    def registerTemplatePatch(self, residue, patch, patchResidueIndex):
        """Register that a particular patch can be used with a particular residue."""
        if residue not in self._templatePatches:
            self._templatePatches[residue] = set()
        self._templatePatches[residue].add((patch, patchResidueIndex))

    def registerScript(self, script):
        """Register a new script to be executed after building the System."""
        self._scripts.append(script)

    def registerTemplateMatcher(self, matcher):
        """Register an object that can override the default logic for matching templates to residues.

        A template matcher is a callable object that can be invoked as::

            template = f(forcefield, residue, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)

        where ``forcefield`` is the ForceField invoking it, ``residue`` is an openmm.app.Residue object,
        ``bondedToAtom[i]`` is the set of atoms bonded to atom index i, and ``ignoreExternalBonds`` and
        ``ignoreExtraParticles`` indicate whether external bonds and extra particules should be considered
        in matching templates.

        It should return a _TemplateData object that matches the residue.  Alternatively it may return
        None, in which case the standard logic will be used to find a template for the residue.

        .. CAUTION:: This method is experimental, and its API is subject to change.
        """
        self._templateMatchers.append(matcher)

    def registerTemplateGenerator(self, generator):
        """Register a residue template generator that can be used to parameterize residues that do not match existing forcefield templates.

        This functionality can be used to add handlers to parameterize small molecules or unnatural/modified residues.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        generator : function
            A function that will be called when a residue is encountered that does not match an existing forcefield template.

        When a residue without a template is encountered, the ``generator`` function is called with:

        ::
           success = generator(forcefield, residue)

        where ``forcefield`` is the calling ``ForceField`` object and ``residue`` is a openmm.app.topology.Residue object.

        ``generator`` must conform to the following API:

        ::
           generator API

           Parameters
           ----------
           forcefield : openmm.app.ForceField
               The ForceField object to which residue templates and/or parameters are to be added.
           residue : openmm.app.Topology.Residue
               The residue topology for which a template is to be generated.

           Returns
           -------
           success : bool
               If the generator is able to successfully parameterize the residue, `True` is returned.
               If the generator cannot parameterize the residue, it should return `False` and not modify `forcefield`.

           The generator should either register a residue template directly with `forcefield.registerResidueTemplate(template)`
           or it should call `forcefield.loadFile(file)` to load residue definitions from an ffxml file.

           It can also use the `ForceField` programmatic API to add additional atom types (via `forcefield.registerAtomType(parameters)`)
           or additional parameters.

        """
        self._templateGenerators.append(generator)

    def _findAtomTypes(self, attrib, num):
        """Parse the attributes on an XML tag to find the set of atom types for each atom it involves.

        Parameters
        ----------
        attrib : dict of attributes
            The dictionary of attributes for an XML parameter tag.
        num : int
            The number of atom specifiers (e.g. 'class1' through 'class4') to extract.

        Returns
        -------
        types : list
            A list of atom types that match.

        """
        types = []
        for i in range(num):
            if num == 1:
                suffix = ''
            else:
                suffix = str(i+1)
            classAttrib = 'class'+suffix
            typeAttrib = 'type'+suffix
            if classAttrib in attrib:
                if typeAttrib in attrib:
                    raise ValueError('Specified both a type and a class for the same atom: '+str(attrib))
                if attrib[classAttrib] not in self._atomClasses:
                    types.append(None) # Unknown atom class
                else:
                    types.append(self._atomClasses[attrib[classAttrib]])
            elif typeAttrib in attrib:
                if attrib[typeAttrib] == '':
                    types.append(self._atomClasses[''])
                elif attrib[typeAttrib] not in self._atomTypes:
                    types.append(None) # Unknown atom type
                else:
                    types.append([attrib[typeAttrib]])
            else:
                types.append(None) # Unknown atom type
        return types

    def _parseTorsion(self, attrib):
        """Parse the node defining a torsion."""
        types = self._findAtomTypes(attrib, 4)
        if None in types:
            return None
        torsion = PeriodicTorsion(types)
        index = 1
        while 'phase%d'%index in attrib:
            torsion.periodicity.append(int(attrib['periodicity%d'%index]))
            torsion.phase.append(_convertParameterToNumber(attrib['phase%d'%index]))
            torsion.k.append(_convertParameterToNumber(attrib['k%d'%index]))
            index += 1
        return torsion

    class _SystemData(object):
        """Inner class used to encapsulate data about the system being created."""
        def __init__(self, topology):
            self.atomType = {}
            self.atomParameters = {}
            self.atomTemplateIndexes = {}
            self.atoms = list(topology.atoms())
            self.excludeAtomWith = [[] for _ in self.atoms]
            self.virtualSites = {}
            self.bonds = [ForceField._BondData(bond[0].index, bond[1].index) for bond in topology.bonds()]
            self.angles = []
            self.propers = []
            self.impropers = []
            self.atomBonds = [[] for _ in self.atoms]
            self.isAngleConstrained = []
            self.constraints = {}
            self.bondedToAtom = [set() for _ in self.atoms]

            # Record which atoms are bonded to each other atom

            for i in range(len(self.bonds)):
                bond = self.bonds[i]
                self.bondedToAtom[bond.atom1].add(bond.atom2)
                self.bondedToAtom[bond.atom2].add(bond.atom1)
                self.atomBonds[bond.atom1].append(i)
                self.atomBonds[bond.atom2].append(i)
            self.bondedToAtom = [sorted(b) for b in self.bondedToAtom]

        def addConstraint(self, system, atom1, atom2, distance):
            """Add a constraint to the system, avoiding duplicate constraints."""
            key = (min(atom1, atom2), max(atom1, atom2))
            if key in self.constraints:
                if self.constraints[key] != distance:
                    raise ValueError('Two constraints were specified between atoms %d and %d with different distances' % (atom1, atom2))
            else:
                self.constraints[key] = distance
                system.addConstraint(atom1, atom2, distance)

        def recordMatchedAtomParameters(self, residue, template, matches):
            """Record parameters for atoms based on having matched a residue to a template."""
            matchAtoms = dict(zip(matches, residue.atoms()))
            for atom, match in zip(residue.atoms(), matches):
                self.atomType[atom] = template.atoms[match].type
                self.atomParameters[atom] = template.atoms[match].parameters
                self.atomTemplateIndexes[atom] = match
                for site in template.virtualSites:
                    if match == site.index:
                        self.virtualSites[atom] = (site, [matchAtoms[i].index for i in site.atoms], matchAtoms[site.excludeWith].index)

    class _TemplateData(object):
        """Inner class used to encapsulate data about a residue template definition."""
        def __init__(self, name):
            self.name = name
            self.atoms = []
            self.atomIndices = {}
            self.virtualSites = []
            self.bonds = []
            self.externalBonds = []
            self.overrideLevel = 0
            self.rigidWater = True
            self.attributes = {}

        def getAtomIndexByName(self, atom_name):
            """Look up an atom index by atom name, providing a helpful error message if not found."""
            index = self.atomIndices.get(atom_name, None)
            if index is not None:
                return index

            # Provide a helpful error message if atom name not found.
            msg =  "Atom name '%s' not found in residue template '%s'.\n" % (atom_name, self.name)
            msg += "Possible atom names are: %s" % str(list(map(lambda x: x.name, self.atoms)))
            raise ValueError(msg)

        def addAtom(self, atom):
            self.atoms.append(atom)
            self.atomIndices[atom.name] = len(self.atoms)-1

        def addBond(self, atom1, atom2):
            """Add a bond between two atoms in a template given their indices in the template."""
            self.bonds.append((atom1, atom2))
            self.atoms[atom1].bondedTo.append(atom2)
            self.atoms[atom2].bondedTo.append(atom1)

        def addBondByName(self, atom1_name, atom2_name):
            """Add a bond between two atoms in a template given their atom names."""
            atom1 = self.getAtomIndexByName(atom1_name)
            atom2 = self.getAtomIndexByName(atom2_name)
            self.addBond(atom1, atom2)

        def addExternalBond(self, atom_index):
            """Designate that an atom in a residue template has an external bond, using atom index within template."""
            self.externalBonds.append(atom_index)
            self.atoms[atom_index].externalBonds += 1

        def addExternalBondByName(self, atom_name):
            """Designate that an atom in a residue template has an external bond, using atom name within template."""
            atom = self.getAtomIndexByName(atom_name)
            self.addExternalBond(atom)

        def areParametersIdentical(self, template2, matchingAtoms, matchingAtoms2):
            """Get whether this template and another one both assign identical atom types and parameters to all atoms.

            Parameters
            ----------
            template2: _TemplateData
                the template to compare this one to
            matchingAtoms: list
                the indices of atoms in this template that atoms of the residue are matched to
            matchingAtoms2: list
                the indices of atoms in template2 that atoms of the residue are matched to
            """
            atoms1 = [self.atoms[m] for m in matchingAtoms]
            atoms2 = [template2.atoms[m] for m in matchingAtoms2]
            if any(a1.type != a2.type or a1.parameters != a2.parameters for a1,a2 in zip(atoms1, atoms2)):
                return False
            # Properly comparing virtual sites really needs a much more complicated analysis.  This simple check
            # could easily fail for templates containing vsites, even if they're actually identical.  Since we
            # currently have no force fields that include both patches and vsites, I'm not going to worry about it now.
            if self.virtualSites != template2.virtualSites:
                return False
            return True

    class _TemplateAtomData(object):
        """Inner class used to encapsulate data about an atom in a residue template definition."""
        def __init__(self, name, type, element, parameters={}):
            self.name = name
            self.type = type
            self.element = element
            self.parameters = parameters
            self.bondedTo = []
            self.externalBonds = 0

    class _BondData(object):
        """Inner class used to encapsulate data about a bond."""
        def __init__(self, atom1, atom2):
            self.atom1 = atom1
            self.atom2 = atom2
            self.isConstrained = False
            self.length = 0.0

    class _VirtualSiteData(object):
        """Inner class used to encapsulate data about a virtual site."""
        def __init__(self, node, atomIndices):
            attrib = node.attrib
            self.type = attrib['type']
            if self.type == 'average2':
                numAtoms = 2
                self.weights = [float(attrib['weight1']), float(attrib['weight2'])]
            elif self.type == 'average3':
                numAtoms = 3
                self.weights = [float(attrib['weight1']), float(attrib['weight2']), float(attrib['weight3'])]
            elif self.type == 'outOfPlane':
                numAtoms = 3
                self.weights = [float(attrib['weight12']), float(attrib['weight13']), float(attrib['weightCross'])]
            elif self.type == 'localCoords':
                numAtoms = 0
                self.originWeights = []
                self.xWeights = []
                self.yWeights = []
                while ('wo%d' % (numAtoms+1)) in attrib:
                    numAtoms += 1
                    self.originWeights.append(float(attrib['wo%d' % numAtoms]))
                    self.xWeights.append(float(attrib['wx%d' % numAtoms]))
                    self.yWeights.append(float(attrib['wy%d' % numAtoms]))
                self.localPos = [float(attrib['p1']), float(attrib['p2']), float(attrib['p3'])]
            else:
                raise ValueError('Unknown virtual site type: %s' % self.type)
            if 'siteName' in attrib:
                self.index = atomIndices[attrib['siteName']]
                self.atoms = [atomIndices[attrib['atomName%d'%(i+1)]] for i in range(numAtoms)]
            else:
                self.index = int(attrib['index'])
                self.atoms = [int(attrib['atom%d'%(i+1)]) for i in range(numAtoms)]
            if 'excludeWith' in attrib:
                self.excludeWith = int(attrib['excludeWith'])
            else:
                self.excludeWith = self.atoms[0]

        def __eq__(self, other):
            if not isinstance(other, ForceField._VirtualSiteData):
                return False
            if self.type != other.type or self.index != other.index or self.atoms != other.atoms or self.excludeWith != other.excludeWith:
                return False
            if self.type in ('average2', 'average3', 'outOfPlane'):
                return self.weights == other.weights
            elif self.type == 'localCoords':
                return self.originWeights == other.originWeights and self.xWeights == other.xWeights and self.yWeights == other.yWeights and self.localPos == other.localPos
            return False

    class _PatchData(object):
        """Inner class used to encapsulate data about a patch definition."""
        def __init__(self, name, numResidues):
            self.name = name
            self.numResidues = numResidues
            self.addedAtoms = [[] for i in range(numResidues)]
            self.changedAtoms = [[] for i in range(numResidues)]
            self.deletedAtoms = []
            self.addedBonds = []
            self.deletedBonds = []
            self.addedExternalBonds = []
            self.deletedExternalBonds = []
            self.allAtomNames = set()
            self.virtualSites = [[] for i in range(numResidues)]
            self.attributes = {}
            self.index = None

        def __lt__(self, other):
            return self.index < other.index

        def createPatchedTemplates(self, templates):
            """Apply this patch to a set of templates, creating new modified ones."""
            if len(templates) != self.numResidues:
                raise ValueError("Patch '%s' expected %d templates, received %d", (self.name, self.numResidues, len(templates)))

            # Construct a new version of each template.

            newTemplates = []
            for index, template in enumerate(templates):
                newTemplate = ForceField._TemplateData("%s-%s" % (template.name, self.name))
                newTemplates.append(newTemplate)

                # Build the list of atoms in it.

                for atom in template.atoms:
                    if not any(deleted.name == atom.name and deleted.residue == index for deleted in self.deletedAtoms):
                        newTemplate.addAtom(ForceField._TemplateAtomData(atom.name, atom.type, atom.element, atom.parameters))
                for atom in self.addedAtoms[index]:
                    if any(a.name == atom.name for a in newTemplate.atoms):
                        raise ValueError("Patch '%s' adds an atom with the same name as an existing atom: %s" % (self.name, atom.name))
                    newTemplate.addAtom(ForceField._TemplateAtomData(atom.name, atom.type, atom.element, atom.parameters))
                oldAtomIndex = dict([(atom.name, i) for i, atom in enumerate(template.atoms)])
                newAtomIndex = dict([(atom.name, i) for i, atom in enumerate(newTemplate.atoms)])
                for atom in self.changedAtoms[index]:
                    if atom.name not in newAtomIndex:
                        raise ValueError("Patch '%s' modifies nonexistent atom '%s' in template '%s'" % (self.name, atom.name, template.name))
                    newTemplate.atoms[newAtomIndex[atom.name]] = ForceField._TemplateAtomData(atom.name, atom.type, atom.element, atom.parameters)

                # Copy over the virtual sites, translating the atom indices.

                indexMap = dict([(oldAtomIndex[name], newAtomIndex[name]) for name in newAtomIndex if name in oldAtomIndex])
                for site in template.virtualSites:
                    if site.index in indexMap and all(i in indexMap for i in site.atoms):
                        newSite = deepcopy(site)
                        newSite.index = indexMap[site.index]
                        newSite.atoms = [indexMap[i] for i in site.atoms]
                        newSite.excludeWith = indexMap[site.excludeWith]
                        newTemplate.virtualSites.append(newSite)

                # Build the lists of bonds and external bonds.

                atomMap = dict([(template.atoms[i], indexMap[i]) for i in indexMap])
                deletedBonds = [(atom1.name, atom2.name) for atom1, atom2 in self.deletedBonds if atom1.residue == index and atom2.residue == index]
                for atom1, atom2 in template.bonds:
                    a1 = template.atoms[atom1]
                    a2 = template.atoms[atom2]
                    if a1 in atomMap and a2 in atomMap and (a1.name, a2.name) not in deletedBonds and (a2.name, a1.name) not in deletedBonds:
                        newTemplate.addBond(atomMap[a1], atomMap[a2])
                deletedExternalBonds = [atom.name for atom in self.deletedExternalBonds if atom.residue == index]
                for atom in template.externalBonds:
                    if template.atoms[atom].name not in deletedExternalBonds:
                        newTemplate.addExternalBond(indexMap[atom])
                for atom1, atom2 in self.addedBonds:
                    if atom1.residue == index and atom2.residue == index:
                        newTemplate.addBondByName(atom1.name, atom2.name)
                    elif atom1.residue == index:
                        newTemplate.addExternalBondByName(atom1.name)
                    elif atom2.residue == index:
                        newTemplate.addExternalBondByName(atom2.name)
                for atom in self.addedExternalBonds:
                    newTemplate.addExternalBondByName(atom.name)

                # Add new virtual sites.

                indexMap = dict((i, newAtomIndex[atom.name]) for i, atom in enumerate(self.addedAtoms[index]+self.changedAtoms[index]))
                for site in self.virtualSites[index]:
                    newSite = deepcopy(site)
                    newSite.index = indexMap[site.index]
                    newSite.atoms = [indexMap[i] for i in site.atoms]
                    newSite.excludeWith = indexMap[site.excludeWith]
                    newTemplate.virtualSites = [site for site in newTemplate.virtualSites if site.index != newSite.index]
                    newTemplate.virtualSites.append(newSite)
            return newTemplates

    class _PatchAtomData(object):
        """Inner class used to encapsulate data about an atom in a patch definition."""
        def __init__(self, description):
            if ':' in description:
                colonIndex = description.find(':')
                self.residue = int(description[:colonIndex])-1
                self.name = description[colonIndex+1:]
            else:
                self.residue = 0
                self.name = description

    class _AtomType(object):
        """Inner class used to record atom types and associated properties."""
        def __init__(self, name, atomClass, mass, element):
            self.name = name
            self.atomClass = atomClass
            self.mass = mass
            self.element = element

    class _AtomTypeParameters(object):
        """Inner class used to record parameter values for atom types."""
        def __init__(self, forcefield, forceName, atomTag, paramNames):
            self.ff = forcefield
            self.forceName = forceName
            self.atomTag = atomTag
            self.paramNames = paramNames
            self.paramsForType = {}
            self.extraParamsForType = {}

        def registerAtom(self, parameters, expectedParams=None):
            if expectedParams is None:
                expectedParams = self.paramNames
            types = self.ff._findAtomTypes(parameters, 1)
            if None not in types:
                values = {}
                extraValues = {}
                for key in parameters:
                    if key in expectedParams:
                        values[key] = _convertParameterToNumber(parameters[key])
                    else:
                        extraValues[key] = parameters[key]
                if len(values) < len(expectedParams):
                    for key in expectedParams:
                        if key not in values:
                            raise ValueError('%s: No value specified for "%s"' % (self.forceName, key))
                for t in types[0]:
                    self.paramsForType[t] = values
                    self.extraParamsForType[t] = extraValues

        def parseDefinitions(self, element):
            """"Load the definitions from an XML element."""
            expectedParams = list(self.paramNames)
            excludedParams = [node.attrib['name'] for node in element.findall('UseAttributeFromResidue')]
            for param in excludedParams:
                if param not in expectedParams:
                    raise ValueError('%s: <UseAttributeFromResidue> specified an invalid attribute: %s' % (self.forceName, param))
                expectedParams.remove(param)
            for atom in element.findall(self.atomTag):
                for param in excludedParams:
                    if param in atom.attrib:
                        raise ValueError('%s: The attribute "%s" appeared in both <%s> and <UseAttributeFromResidue> tags' % (self.forceName, param, self.atomTag))
                self.registerAtom(atom.attrib, expectedParams)

        def getAtomParameters(self, atom, data):
            """Get the parameter values for a particular atom."""
            t = data.atomType[atom]
            p = data.atomParameters[atom]
            if t in self.paramsForType:
                values = self.paramsForType[t]
                result = [None]*len(self.paramNames)
                for i, name in enumerate(self.paramNames):
                    if name in values:
                        result[i] = values[name]
                    elif name in p:
                        result[i] = p[name]
                    else:
                        raise ValueError('%s: No value specified for "%s"' % (self.forceName, name))
                return result
            else:
                raise ValueError('%s: No parameters defined for atom type %s' % (self.forceName, t))

        def getExtraParameters(self, atom, data):
            """Get extra parameter values for an atom that appeared in the <Atom> tag but were not included in paramNames."""
            t = data.atomType[atom]
            if t in self.paramsForType:
                return self.extraParamsForType[t]
            else:
                raise ValueError('%s: No parameters defined for atom type %s' % (self.forceName, t))


    def _getResidueTemplateMatches(self, res, bondedToAtom, templateSignatures=None, ignoreExternalBonds=False, ignoreExtraParticles=False):
        """Return the templates that match a residue, or None if none are found.

        Parameters
        ----------
        res : Topology.Residue
            The residue for which template matches are to be retrieved.
        bondedToAtom : list of set of int
            bondedToAtom[i] is the set of atoms bonded to atom index i

        Returns
        -------
        template : _TemplateData
            The matching forcefield residue template, or None if no matches are found.
        matches : list
            a list specifying which atom of the template each atom of the residue
            corresponds to, or None if it does not match the template

        """
        template = None
        matches = None
        for matcher in self._templateMatchers:
            template = matcher(self, res, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
            if template is not None:
                match = compiled.matchResidueToTemplate(res, template, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                if match is None:
                    raise ValueError('A custom template matcher returned a template for residue %d (%s), but it does not match the residue.' % (res.index, res.name))
                return [template, match]
        if templateSignatures is None:
            templateSignatures = self._templateSignatures
        signature = _createResidueSignature([atom.element for atom in res.atoms()])
        if signature in templateSignatures:
            allMatches = []
            for t in templateSignatures[signature]:
                match = compiled.matchResidueToTemplate(res, t, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                if match is not None:
                    allMatches.append((t, match))
            if len(allMatches) == 1:
                template = allMatches[0][0]
                matches = allMatches[0][1]
            elif len(allMatches) > 1:
                # We found multiple matches.  This is OK if and only if they assign identical types and parameters to all atoms.
                t1, m1 = allMatches[0]
                for t2, m2 in allMatches[1:]:
                    if not t1.areParametersIdentical(t2, m1, m2):
                        raise Exception('Multiple non-identical matching templates found for residue %d (%s): %s.' % (res.index, res.name, ', '.join(match[0].name for match in allMatches)))
                template = allMatches[0][0]
                matches = allMatches[0][1]
        return [template, matches]

    def _buildBondedToAtomList(self, topology):
        """Build a list of which atom indices are bonded to each atom.

        Parameters
        ----------
        topology : Topology
            The Topology whose bonds are to be indexed.

        Returns
        -------
        bondedToAtom : list of list of int
            bondedToAtom[index] is the list of atom indices bonded to atom `index`

        """
        bondedToAtom = [set() for _ in topology.atoms()]
        for (atom1, atom2) in topology.bonds():
            bondedToAtom[atom1.index].add(atom2.index)
            bondedToAtom[atom2.index].add(atom1.index)
        bondedToAtom = [sorted(b) for b in bondedToAtom]
        return bondedToAtom

    def getUnmatchedResidues(self, topology, residueTemplates=dict()):
        """Return a list of Residue objects from specified topology for which no forcefield templates are available.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        topology : Topology
            The Topology whose residues are to be checked against the forcefield residue templates.
        residueTemplates : dict=dict()
            Specifies which template to use for particular residues.  The keys should be Residue
            objects from the Topology, and the values should be the names of the templates to
            use for them.  This is useful when a ForceField contains multiple templates that
            can match the same residue (e.g Fe2+ and Fe3+ templates in the ForceField for a
            monoatomic iron ion in the Topology).

        Returns
        -------
        unmatched_residues : list of Residue
            List of Residue objects from `topology` for which no forcefield residue templates are available.
            Note that multiple instances of the same residue appearing at different points in the topology may be returned.

        This method may be of use in generating missing residue templates or diagnosing parameterization failures.
        """
        # Find the template matching each residue, compiling a list of residues for which no templates are available.
        bondedToAtom = self._buildBondedToAtomList(topology)
        unmatched_residues = list() # list of unmatched residues
        for res in topology.residues():
            if res in residueTemplates:
                # Make sure the specified template matches.
                template = self._templates[residueTemplates[res]]
                matches = compiled.matchResidueToTemplate(res, template, bondedToAtom, False, False)
            else:
                # Attempt to match one of the existing templates.
                [template, matches] = self._getResidueTemplateMatches(res, bondedToAtom)
            if matches is None:
                # No existing templates match.
                unmatched_residues.append(res)

        return unmatched_residues

    def getMatchingTemplates(self, topology, ignoreExternalBonds=False):
        """Return a list of forcefield residue templates matching residues in the specified topology.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        topology : Topology
            The Topology whose residues are to be checked against the forcefield residue templates.
        ignoreExternalBonds : bool=False
            If true, ignore external bonds when matching residues to templates.
        Returns
        -------
        templates : list of _TemplateData
            List of forcefield residue templates corresponding to residues in the topology.
            templates[index] is template corresponding to residue `index` in topology.residues()

        This method may be of use in debugging issues related to parameter assignment.
        """
        # Find the template matching each residue, compiling a list of residues for which no templates are available.
        bondedToAtom = self._buildBondedToAtomList(topology)
        templates = list() # list of templates matching the corresponding residues
        for residue in topology.residues():
            # Attempt to match one of the existing templates.
            [template, matches] = self._getResidueTemplateMatches(residue, bondedToAtom, ignoreExternalBonds=ignoreExternalBonds)
            # Raise an exception if we have found no templates that match.
            if matches is None:
                raise ValueError('No template found for chainid <%s> resid <%s> resname <%s> (residue index within topology %d).\n%s' % (residue.chain.id, residue.id, residue.name, residue.index, _findMatchErrors(self, residue)))
            else:
                templates.append(template)

        return templates

    def generateTemplatesForUnmatchedResidues(self, topology):
        """Generate forcefield residue templates for residues in specified topology for which no forcefield templates are available.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        topology : Topology
            The Topology whose residues are to be checked against the forcefield residue templates.

        Returns
        -------
        templates : list of _TemplateData
            List of forcefield residue templates corresponding to residues in `topology` for which no forcefield templates are currently available.
            Atom types will be set to `None`, but template name, atom names, elements, and connectivity will be taken from corresponding Residue objects.
        residues : list of Residue
            List of Residue objects that were used to generate the templates.
            `residues[index]` is the Residue that was used to generate the template `templates[index]`

        """
        # Get a non-unique list of unmatched residues.
        unmatched_residues = self.getUnmatchedResidues(topology)
        # Generate a unique list of unmatched residues by comparing fingerprints.
        bondedToAtom = self._buildBondedToAtomList(topology)
        unique_unmatched_residues = list() # list of unique unmatched Residue objects from topology
        templates = list() # corresponding _TemplateData templates
        signatures = set()
        for residue in unmatched_residues:
            signature = _createResidueSignature([ atom.element for atom in residue.atoms() ])
            template = _createResidueTemplate(residue)
            is_unique = True
            if signature in signatures:
                # Signature is the same as an existing residue; check connectivity.
                for check_residue in unique_unmatched_residues:
                    matches = compiled.matchResidueToTemplate(check_residue, template, bondedToAtom, False)
                    if matches is not None:
                        is_unique = False
            if is_unique:
                # Residue is unique.
                unique_unmatched_residues.append(residue)
                signatures.add(signature)
                templates.append(template)

        return [templates, unique_unmatched_residues]

    def createSystem(self, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=None, removeCMMotion=True, hydrogenMass=None, residueTemplates=dict(),
                     ignoreExternalBonds=False, switchDistance=None, flexibleConstraints=False, drudeMass=0.4*unit.amu, **args):
        """Construct an OpenMM System representing a Topology with this force field.

        Parameters
        ----------
        topology : Topology
            The Topology for which to create a System
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, or LJPME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=None
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument.  If None (the default), it uses the
            default behavior for this force field's water model.
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep
            their total mass the same.  If rigidWater is used to make water molecules
            rigid, then water hydrogens are not altered.
        residueTemplates : dict=dict()
            Specifies which template to use for particular residues.  The keys should be Residue
            objects from the Topology, and the values should be the names of the templates to
            use for them.  This is useful when a ForceField contains multiple templates that
            can match the same residue (e.g Fe2+ and Fe3+ templates in the ForceField for a
            monoatomic iron ion in the Topology).
        ignoreExternalBonds : boolean=False
            If true, ignore external bonds when matching residues to templates.  This is
            useful when the Topology represents one piece of a larger molecule, so chains are
            not terminated properly.  This option can create ambiguities where multiple
            templates match the same residue.  If that happens, use the residueTemplates
            argument to specify which one to use.
        switchDistance : float=None
            The distance at which the potential energy switching function is turned on for
            Lennard-Jones interactions. If this is None, no switching function will be used.
        flexibleConstraints : boolean=False
            If True, parameters for constrained degrees of freedom will be added to the System
        drudeMass : mass=0.4*amu
            The mass to use for Drude particles.  Any mass added to a Drude particle is
            subtracted from its parent atom to keep their total mass the same.
        args
            Arbitrary additional keyword arguments may also be specified.
            This allows extra parameters to be specified that are specific to
            particular force fields.

        Returns
        -------
        system
            the newly created System
        """
        args['switchDistance'] = switchDistance
        args['flexibleConstraints'] = flexibleConstraints
        args['drudeMass'] = drudeMass
        args = ArgTracker(args)
        data = ForceField._SystemData(topology)
        rigidResidue = [False]*topology.getNumResidues()

        # Find the template matching each residue and assign atom types.

        templateForResidue = self._matchAllResiduesToTemplates(data, topology, residueTemplates, ignoreExternalBonds)
        for res in topology.residues():
            if res.name == 'HOH':
                # Determine whether this should be a rigid water.

                if rigidWater is None:
                    rigidResidue[res.index] = templateForResidue[res.index].rigidWater
                elif rigidWater:
                    rigidResidue[res.index] = True

        # Create the System and add atoms

        sys = mm.System()
        for atom in topology.atoms():
            # Look up the atom type name, returning a helpful error message if it cannot be found.
            if atom not in data.atomType:
                raise Exception("Could not identify atom type for atom '%s'." % str(atom))
            typename = data.atomType[atom]

            # Look up the type name in the list of registered atom types, returning a helpful error message if it cannot be found.
            if typename not in self._atomTypes:
                msg  = "Could not find typename '%s' for atom '%s' in list of known atom types.\n" % (typename, str(atom))
                msg += "Known atom types are: %s" % str(self._atomTypes.keys())
                raise Exception(msg)

            # Add the particle to the OpenMM system.
            mass = self._atomTypes[typename].mass
            sys.addParticle(mass)

        # Adjust hydrogen masses if requested.

        if hydrogenMass is not None:
            if not unit.is_quantity(hydrogenMass):
                hydrogenMass *= unit.dalton
            for atom1, atom2 in topology.bonds():
                if atom1.element is elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element is elem.hydrogen and atom1.element not in (elem.hydrogen, None) and not rigidResidue[atom2.residue.index]:
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Set periodic boundary conditions.

        boxVectors = topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
        elif nonbondedMethod not in [NoCutoff, CutoffNonPeriodic]:
            raise ValueError('Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')

        # Make a list of all unique angles

        uniqueAngles = set()
        for bond in data.bonds:
            for atom in data.bondedToAtom[bond.atom1]:
                if atom != bond.atom2:
                    if atom < bond.atom2:
                        uniqueAngles.add((atom, bond.atom1, bond.atom2))
                    else:
                        uniqueAngles.add((bond.atom2, bond.atom1, atom))
            for atom in data.bondedToAtom[bond.atom2]:
                if atom != bond.atom1:
                    if atom > bond.atom1:
                        uniqueAngles.add((bond.atom1, bond.atom2, atom))
                    else:
                        uniqueAngles.add((atom, bond.atom2, bond.atom1))
        data.angles = sorted(list(uniqueAngles))

        # Make a list of all unique proper torsions

        uniquePropers = set()
        for angle in data.angles:
            for atom in data.bondedToAtom[angle[0]]:
                if atom not in angle:
                    if atom < angle[2]:
                        uniquePropers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniquePropers.add((angle[2], angle[1], angle[0], atom))
            for atom in data.bondedToAtom[angle[2]]:
                if atom not in angle:
                    if atom > angle[0]:
                        uniquePropers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniquePropers.add((atom, angle[2], angle[1], angle[0]))
        data.propers = sorted(list(uniquePropers))

        # Make a list of all unique improper torsions

        for atom in range(len(data.bondedToAtom)):
            bondedTo = data.bondedToAtom[atom]
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
                bond.isConstrained = atom1.element is elem.hydrogen or atom2.element is elem.hydrogen
        for bond in data.bonds:
            atom1 = data.atoms[bond.atom1]
            atom2 = data.atoms[bond.atom2]
            if rigidResidue[atom1.residue.index] and rigidResidue[atom2.residue.index]:
                bond.isConstrained = True

        # Identify angles that should be implemented with constraints

        if constraints == HAngles:
            for angle in data.angles:
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                numH = 0
                if atom1.element is elem.hydrogen:
                    numH += 1
                if atom3.element is elem.hydrogen:
                    numH += 1
                data.isAngleConstrained.append(numH == 2 or (numH == 1 and atom2.element is elem.oxygen))
        else:
            data.isAngleConstrained = len(data.angles)*[False]
        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if rigidResidue[atom1.residue.index] and rigidResidue[atom2.residue.index] and rigidResidue[atom3.residue.index]:
                data.isAngleConstrained[i] = True

        # Add virtual sites

        for atom in data.virtualSites:
            (site, atoms, excludeWith) = data.virtualSites[atom]
            index = atom.index
            data.excludeAtomWith[excludeWith].append(index)
            if site.type == 'average2':
                sys.setVirtualSite(index, mm.TwoParticleAverageSite(atoms[0], atoms[1], site.weights[0], site.weights[1]))
            elif site.type == 'average3':
                sys.setVirtualSite(index, mm.ThreeParticleAverageSite(atoms[0], atoms[1], atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'outOfPlane':
                sys.setVirtualSite(index, mm.OutOfPlaneSite(atoms[0], atoms[1], atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'localCoords':
                sys.setVirtualSite(index, mm.LocalCoordinatesSite(atoms, site.originWeights, site.xWeights, site.yWeights, site.localPos))

        # Add forces to the System

        for force in self._forces:
            force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        # Let force generators do postprocessing

        for force in self._forces:
            if 'postprocessSystem' in dir(force):
                force.postprocessSystem(sys, data, args)

        # Execute scripts found in the XML files.

        for script in self._scripts:
            exec(script, locals())
        args.checkArgs(self.createSystem)
        return sys


    def _matchAllResiduesToTemplates(self, data, topology, residueTemplates, ignoreExternalBonds, ignoreExtraParticles=False, recordParameters=True):
        """Return a list of which template matches each residue in the topology, and assign atom types."""
        templateForResidue = [None]*topology.getNumResidues()
        unmatchedResidues = []
        for chain in topology.chains():
            for res in chain.residues():
                if res in residueTemplates:
                    tname = residueTemplates[res]
                    template = self._templates[tname]
                    matches = compiled.matchResidueToTemplate(res, template, data.bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                    if matches is None:
                        raise Exception('User-supplied template %s does not match the residue %d (%s)' % (tname, res.index, res.name))
                else:
                    # Attempt to match one of the existing templates.
                    [template, matches] = self._getResidueTemplateMatches(res, data.bondedToAtom, ignoreExternalBonds=ignoreExternalBonds, ignoreExtraParticles=ignoreExtraParticles)
                if matches is None:
                    unmatchedResidues.append(res)
                else:
                    if recordParameters:
                        data.recordMatchedAtomParameters(res, template, matches)
                    templateForResidue[res.index] = template

        # Try to apply patches to find matches for any unmatched residues.

        if len(unmatchedResidues) > 0:
            unmatchedResidues = _applyPatchesToMatchResidues(self, data, unmatchedResidues, templateForResidue, data.bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)

        # If we still haven't found a match for a residue, attempt to use residue template generators to create
        # new templates (and potentially atom types/parameters).

        for res in unmatchedResidues:
            # A template might have been generated on an earlier iteration of this loop.
            [template, matches] = self._getResidueTemplateMatches(res, data.bondedToAtom, ignoreExternalBonds=ignoreExternalBonds, ignoreExtraParticles=ignoreExtraParticles)
            if matches is None:
                # Try all generators.
                for generator in self._templateGenerators:
                    if generator(self, res):
                        # This generator has registered a new residue template that should match.
                        [template, matches] = self._getResidueTemplateMatches(res, data.bondedToAtom, ignoreExternalBonds=ignoreExternalBonds, ignoreExtraParticles=ignoreExtraParticles)
                        if matches is None:
                            # Something went wrong because the generated template does not match the residue signature.
                            raise Exception('The residue handler %s indicated it had correctly parameterized residue %s, but the generated template did not match the residue signature.' % (generator.__class__.__name__, str(res)))
                        else:
                            # We successfully generated a residue template.  Break out of the for loop.
                            break
            if matches is None:
                raise ValueError('No template found for residue %d (%s).  %s  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#template' % (res.index, res.name, _findMatchErrors(self, res)))
            else:
                if recordParameters:
                    data.recordMatchedAtomParameters(res, template, matches)
                templateForResidue[res.index] = template
        return templateForResidue


def _findBondsForExclusions(data, sys):
    """Create a list of bonds to use when identifying exclusions."""
    bondIndices = []
    for bond in data.bonds:
        bondIndices.append((bond.atom1, bond.atom2))

    # If a virtual site does *not* share exclusions with another atom, add a bond between it and its first parent atom.

    for i in range(sys.getNumParticles()):
        if sys.isVirtualSite(i):
            (site, atoms, excludeWith) = data.virtualSites[data.atoms[i]]
            if excludeWith is None:
                bondIndices.append((i, site.getParticle(0)))

    # Certain particles, such as lone pairs and Drude particles, share exclusions with a parent atom.
    # If the parent atom does not interact with an atom, the child particle does not either.

    for atom1, atom2 in bondIndices:
        for child1 in data.excludeAtomWith[atom1]:
            bondIndices.append((child1, atom2))
            for child2 in data.excludeAtomWith[atom2]:
                bondIndices.append((child1, child2))
        for child2 in data.excludeAtomWith[atom2]:
            bondIndices.append((atom1, child2))
    for atom in data.atoms:
        for child in data.excludeAtomWith[atom.index]:
            bondIndices.append((child, atom.index))
    return bondIndices

def _findExclusions(bondIndices, maxSeparation, numAtoms):
    """Identify pairs of atoms in the same molecule separated by no more than maxSeparation bonds."""
    bondedTo = [set() for i in range(numAtoms)]
    for i, j in bondIndices:
        bondedTo[i].add(j)
        bondedTo[j].add(i)

    # Identify all neighbors of each atom with each separation.

    bondedWithSeparation = [bondedTo]
    for i in range(maxSeparation-1):
        lastBonds = bondedWithSeparation[-1]
        newBonds = deepcopy(lastBonds)
        for atom in range(numAtoms):
            for a1 in lastBonds[atom]:
                for a2 in bondedTo[a1]:
                    newBonds[atom].add(a2)
        bondedWithSeparation.append(newBonds)

    # Build the list of pairs.

    pairs = []
    for atom in range(numAtoms):
        for otherAtom in bondedWithSeparation[-1][atom]:
            if otherAtom > atom:
                # Determine the minimum number of bonds between them.
                sep = maxSeparation
                for i in reversed(range(maxSeparation-1)):
                    if otherAtom in bondedWithSeparation[i][atom]:
                        sep -= 1
                    else:
                        break
                pairs.append((atom, otherAtom, sep))
    return pairs


def _findGroups(bondedTo):
    """Given bonds that connect atoms, identify the connected groups."""
    atomGroup = [None]*len(bondedTo)
    numGroups = 0
    for i in range(len(bondedTo)):
        if atomGroup[i] is None:
            # Start a new group.

            atomStack = [i]
            neighborStack = [0]
            group = numGroups
            numGroups += 1

            # Recursively tag all the bonded atoms.

            while len(atomStack) > 0:
                atom = atomStack[-1]
                atomGroup[atom] = group
                while neighborStack[-1] < len(bondedTo[atom]) and atomGroup[bondedTo[atom][neighborStack[-1]]] is not None:
                    neighborStack[-1] += 1
                if neighborStack[-1] < len(bondedTo[atom]):
                    atomStack.append(bondedTo[atom][neighborStack[-1]])
                    neighborStack.append(0)
                else:
                    atomStack.pop()
                    neighborStack.pop()
    return atomGroup

def _countResidueAtoms(elements):
    """Count the number of atoms of each element in a residue."""
    counts = {}
    for element in elements:
        if element in counts:
            counts[element] += 1
        else:
            counts[element] = 1
    return counts


def _createResidueSignature(elements):
    """Create a signature for a residue based on the elements of the atoms it contains."""
    counts = _countResidueAtoms(elements)
    sig = []
    for c in counts:
        if c is not None:
            sig.append((c, counts[c]))
    sig.sort(key=lambda x: -x[0].mass)

    # Convert it to a string.

    s = ''
    for element, count in sig:
        s += element.symbol+str(count)
    return s


def _applyPatchesToMatchResidues(forcefield, data, residues, templateForResidue, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles):
    """Try to apply patches to find matches for residues."""
    # Start by creating all templates than can be created by applying a combination of one-residue patches
    # to a single template.  The number of these is usually not too large, and they often cover a large fraction
    # of residues.

    patchedTemplateSignatures = {}
    patchedTemplates = {}
    for name, template in forcefield._templates.items():
        if name in forcefield._templatePatches:
            patches = sorted([forcefield._patches[patchName] for patchName, patchResidueIndex in forcefield._templatePatches[name] if forcefield._patches[patchName].numResidues == 1])
            if len(patches) > 0:
                newTemplates = []
                patchedTemplates[name] = newTemplates
                _generatePatchedSingleResidueTemplates(template, patches, 0, newTemplates, set())
                for patchedTemplate in newTemplates:
                    signature = _createResidueSignature([atom.element for atom in patchedTemplate.atoms])
                    if signature in patchedTemplateSignatures:
                        patchedTemplateSignatures[signature].append(patchedTemplate)
                    else:
                        patchedTemplateSignatures[signature] = [patchedTemplate]

    # Now see if any of those templates matches any of the residues.

    unmatchedResidues = []
    for res in residues:
        [template, matches] = forcefield._getResidueTemplateMatches(res, bondedToAtom, patchedTemplateSignatures, ignoreExternalBonds, ignoreExtraParticles)
        if matches is None:
            unmatchedResidues.append(res)
        else:
            data.recordMatchedAtomParameters(res, template, matches)
            templateForResidue[res.index] = template
    if len(unmatchedResidues) == 0:
        return []

    # We need to consider multi-residue patches.  This can easily lead to a combinatorial explosion, so we make a simplifying
    # assumption: that no residue is affected by more than one multi-residue patch (in addition to any number of single-residue
    # patches).  Record all multi-residue patches, and the templates they can be applied to.

    patches = {}
    maxPatchSize = 0
    for patch in forcefield._patches.values():
        if patch.numResidues > 1:
            patches[patch.name] = [[] for i in range(patch.numResidues)]
            maxPatchSize = max(maxPatchSize, patch.numResidues)
    if maxPatchSize == 0:
        return unmatchedResidues # There aren't any multi-residue patches
    for templateName in forcefield._templatePatches:
        for patchName, patchResidueIndex in forcefield._templatePatches[templateName]:
            if patchName in patches:
                # The patch should accept this template, *and* all patched versions of it generated above.
                patches[patchName][patchResidueIndex].append(forcefield._templates[templateName])
                if templateName in patchedTemplates:
                    patches[patchName][patchResidueIndex] += patchedTemplates[templateName]

    # Record which unmatched residues are bonded to each other.

    bonds = set()
    topology = residues[0].chain.topology
    for atom1, atom2 in topology.bonds():
        if atom1.residue != atom2.residue:
            res1 = atom1.residue
            res2 = atom2.residue
            if res1 in unmatchedResidues and res2 in unmatchedResidues:
                bond = tuple(sorted((res1, res2), key=lambda x: x.index))
                if bond not in bonds:
                    bonds.add(bond)

    # Identify clusters of unmatched residues that are all bonded to each other.  These are the ones we'll
    # try to apply multi-residue patches to.

    clusterSize = 2
    clusters = bonds
    while clusterSize <= maxPatchSize:
        # Try to apply patches to clusters of this size.

        for patchName in patches:
            patch = forcefield._patches[patchName]
            if patch.numResidues == clusterSize:
                matchedClusters = _matchToMultiResiduePatchedTemplates(data, clusters, patch, patches[patchName], bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                for cluster in matchedClusters:
                    for residue in cluster:
                        unmatchedResidues.remove(residue)
                bonds = set(bond for bond in bonds if bond[0] in unmatchedResidues and bond[1] in unmatchedResidues)

        # Now extend the clusters to find ones of the next size up.

        largerClusters = set()
        for cluster in clusters:
            for bond in bonds:
                if bond[0] in cluster and bond[1] not in cluster:
                    newCluster = tuple(sorted(cluster+(bond[1],), key=lambda x: x.index))
                    largerClusters.add(newCluster)
                elif bond[1] in cluster and bond[0] not in cluster:
                    newCluster = tuple(sorted(cluster+(bond[0],), key=lambda x: x.index))
                    largerClusters.add(newCluster)
        if len(largerClusters) == 0:
            # There are no clusters of this size or larger
            break
        clusters = largerClusters
        clusterSize += 1

    return unmatchedResidues


def _generatePatchedSingleResidueTemplates(template, patches, index, newTemplates, alteredAtoms):
    """Apply all possible combinations of a set of single-residue patches to a template."""
    try:
        if len(alteredAtoms.intersection(patches[index].allAtomNames)) > 0:
            # This patch would alter an atom that another patch has already altered,
            # so don't apply it.
            patchedTemplate = None
        else:
            patchedTemplate = patches[index].createPatchedTemplates([template])[0]
            newTemplates.append(patchedTemplate)
    except:
        # This probably means the patch is inconsistent with another one that has already been applied,
        # so just ignore it.
        patchedTemplate = None

    # Call this function recursively to generate combinations of patches.

    if index+1 < len(patches):
        _generatePatchedSingleResidueTemplates(template, patches, index+1, newTemplates, alteredAtoms)
        if patchedTemplate is not None:
            newAlteredAtoms = alteredAtoms.union(patches[index].allAtomNames)
            _generatePatchedSingleResidueTemplates(patchedTemplate, patches, index+1, newTemplates, newAlteredAtoms)


def _matchToMultiResiduePatchedTemplates(data, clusters, patch, residueTemplates, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles):
    """Apply a multi-residue patch to templates, then try to match them against clusters of residues."""
    matchedClusters = []
    selectedTemplates = [None]*patch.numResidues
    _applyMultiResiduePatch(data, clusters, patch, residueTemplates, selectedTemplates, 0, matchedClusters, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
    return matchedClusters


def _applyMultiResiduePatch(data, clusters, patch, candidateTemplates, selectedTemplates, index, matchedClusters, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles):
    """This is called recursively to apply a multi-residue patch to all possible combinations of templates."""

    if index < patch.numResidues:
        for template in candidateTemplates[index]:
            selectedTemplates[index] = template
            _applyMultiResiduePatch(data, clusters, patch, candidateTemplates, selectedTemplates, index+1, matchedClusters, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
    else:
        # We're at the deepest level of the recursion.  We've selected a template for each residue, so apply the patch,
        # then try to match it against clusters.

        try:
            patchedTemplates = patch.createPatchedTemplates(selectedTemplates)
        except:
            # This probably means the patch is inconsistent with another one that has already been applied,
            # so just ignore it.
            return
        newlyMatchedClusters = []
        for cluster in clusters:
            for residues in itertools.permutations(cluster):
                residueMatches = []
                for residue, template in zip(residues, patchedTemplates):
                    matches = compiled.matchResidueToTemplate(residue, template, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles)
                    if matches is None:
                        residueMatches = None
                        break
                    else:
                        residueMatches.append(matches)
                if residueMatches is not None:
                    # Each residue individually matches.  Now make sure they're bonded in the correct way.

                    bondsMatch = True
                    for a1, a2 in patch.addedBonds:
                        res1 = a1.residue
                        res2 = a2.residue
                        if res1 != res2:
                            # The patch adds a bond between residues.  Make sure that bond exists.

                            atoms1 = patchedTemplates[res1].atoms
                            atoms2 = patchedTemplates[res2].atoms
                            index1 = next(i for i in range(len(atoms1)) if atoms1[residueMatches[res1][i]].name == a1.name)
                            index2 = next(i for i in range(len(atoms2)) if atoms2[residueMatches[res2][i]].name == a2.name)
                            atom1 = list(residues[res1].atoms())[index1]
                            atom2 = list(residues[res2].atoms())[index2]
                            bondsMatch &= atom2.index in bondedToAtom[atom1.index]
                    if bondsMatch:
                        # We successfully matched the template to the residues.  Record the parameters.

                        for i in range(patch.numResidues):
                            data.recordMatchedAtomParameters(residues[i], patchedTemplates[i], residueMatches[i])
                        newlyMatchedClusters.append(cluster)
                        break

        # Record which clusters were successfully matched.

        matchedClusters += newlyMatchedClusters
        for cluster in newlyMatchedClusters:
            clusters.remove(cluster)


def _findMatchErrors(forcefield, res):
    """Try to guess why a residue failed to match any template and return an error message."""

    def counterSubtract(counter1, counter2):
        """Subtracts two Counter objects (equivalent to counter1 - counter2 but preserves negatives)."""
        difference = counter1.copy()
        difference.subtract(counter2)
        return difference

    def elemToNum(elem):
        """Returns the atomic number of an element or 0 for None."""
        return 0 if elem is None else elem.atomic_number

    def makeIntBondSpec(bond):
        """Converts an internal bond to (x, y) where x and y are atomic numbers of the bonded atoms' elements, and x < y."""
        elemNum1 = elemToNum(bond[0].element)
        elemNum2 = elemToNum(bond[1].element)
        return min(elemNum1, elemNum2), max(elemNum1, elemNum2)

    def makeExtBondSpec(bond):
        """Converts an external bond to the atomic number of the element of the non-external atom."""
        return elemToNum(bond.atom1.element) if bond.atom1.residue is res else elemToNum(bond.atom2.element)

    def bestMatchAtom(templatesDiffs, templateNames, useHeavy):
        """Finds the best matching templates based on atoms, optionally ignoring non-heavy atoms."""
        bestMatches = []
        bestScore = None
        for templateName in templateNames:
            # Find the (signed) differences in atom counts for all elements,
            # optionally considering only heavy atoms.
            allDiffs = templatesDiffs[templateName].items()
            diffs = [(elemNum, diff) for elemNum, diff in allDiffs if not useHeavy or elemNum > 1]
            # Build a score as (x, y), where y is the sum of the magnitudes of
            # the differences in element counts, and x is a flag set if the
            # residue has more of any element than the template.  Minimizing the
            # score will, when these tuples are compared, favor templates where
            # the residue is missing atoms first, and only if there are none,
            # include templates for which the residue has extra atoms.
            score = (any(diff > 0 for _, diff in allDiffs), sum(abs(diff) for _, diff in diffs))
            if bestScore is None or score <= bestScore:
                # Keep a list of every template with the lowest score seen.
                if score != bestScore:
                    # If bestScore is None or score < bestScore, clear the list.
                    bestMatches.clear()
                    bestScore = score
                bestMatches.append(templateName)
        return bestMatches, bestScore

    def bestMatchBond(templatesDiffs, templateNames):
        """Finds the best matching templates based on bonds."""
        bestMatches = []
        bestScore = None
        for templateName in templateNames:
            diffs = templatesDiffs[templateName].items()
            # Favor templates where the residue is missing bonds first, and if
            # there are none, include templates where it has extra bonds.
            score = (any(diff > 0 for _, diff in diffs), sum(abs(diff) for _, diff in diffs))
            if bestScore is None or score <= bestScore:
                if score != bestScore:
                    bestMatches.clear()
                    bestScore = score
                bestMatches.append(templateName)
        return bestMatches, bestScore

    def joinMessages(messages):
        """Produce a human-readable message from the individual message strings passed."""
        messages = list(messages)
        if len(messages) < 3:
            return ' and '.join(messages)
        messages[-1] = f'and {messages[-1]}'
        return ', '.join(messages)

    def formatDiffMessage(diffs, formatter):
        """Formats a message describing a difference between a residue and a template."""
        missing, extra = [], []
        for key, diff in diffs.items():
            if diff < 0:
                missing.append((key, -diff))
            if diff > 0:
                extra.append((key, diff))
        messages = []
        if missing:
            message = joinMessages(f'{diff} {formatter(key, diff)}' for key, diff in sorted(missing))
            messages.append(f'is missing {message}')
        if extra:
            message = joinMessages(f'{diff} {formatter(key, diff)}' for key, diff in sorted(extra))
            messages.append(f'has {message} too many')
        return joinMessages(messages)

    def formatAtomDiff(key, diff):
        """Formats a string describing an element associated with a different number of atoms."""
        return (f'{elem.Element.getByAtomicNumber(key).symbol} atom' if key else 'extra site') + ('' if diff == 1 else 's')

    def formatBondDiff(key, diff):
        """Formats a string describing elements associated with a different number of bonds."""
        name1 = elem.Element.getByAtomicNumber(key[0]).symbol if key[0] else 'extra site'
        name2 = elem.Element.getByAtomicNumber(key[1]).symbol if key[1] else 'extra site'
        return f'{name1}-{name2} bond' + ('' if diff == 1 else 's')

    def pickBestMatch(bestMatches):
        """If there are multiple best-scoring matches, pick one with the closest name to the residue.  Call only with a non-empty list."""
        if not res.name:
            return bestMatches[0]
        return max(bestMatches, key=lambda match: SequenceMatcher(a=res.name.strip(), b=match.strip(), autojunk=False).ratio())

    # First check to see if there are no templates in the force field.

    if not forcefield._templates:
        return f'The force field contains no residue templates.'

    # Get all elements in all templates in the force field and see if this
    # residue uses an element not supported.  Otherwise, prepare fingerprints of
    # the residue and templates based on counts of the elements of their atoms.

    supportedElements = set(elemToNum(atom.element) for template in forcefield._templates.values() for atom in template.atoms)
    residueAtomCounts = Counter(elemToNum(atom.element) for atom in res.atoms())
    unsupportedElements = set(residueAtomCounts.keys()) - supportedElements
    if unsupportedElements:
        unsupportedMessage = joinMessages(formatAtomDiff(elemNum, 0) for elemNum in sorted(unsupportedElements))
        return f'The residue contains {unsupportedMessage}, which are not supported by any template in the force field.'

    # It will be useful to provide a special message if this is a terminal
    # residue of the chain, since this is a common cause of topology problems.

    chainResidues = list(res.chain.residues())
    if len(chainResidues) > 1 and (res == chainResidues[0] or res == chainResidues[-1]):
        terminalMessage = '  Is the chain terminated in a way that is unsupported by the force field?'
    else:
        terminalMessage = ''

    def makeTemplateAtomDiff(templateName):
        """Prepares a Counter describing the difference between a residue's and a template's atoms."""
        template = forcefield._templates[templateName]
        return counterSubtract(residueAtomCounts, Counter(elemToNum(atom.element) for atom in template.atoms))

    templatesAtomDiffs = {templateName: makeTemplateAtomDiff(templateName) for templateName in forcefield._templates}

    # Compare the residue with templates, first based on heavy atoms, then if
    # templates are found where all heavy atoms match, on all atoms.

    bestMatches = templatesAtomDiffs.keys()

    bestMatches, bestScore = bestMatchAtom(templatesAtomDiffs, bestMatches, True)
    if bestMatches and bestScore[1]:
        # If there are matches found and the best score's sum is non-zero
        # (an imperfect match), return.  Otherwise, if the best matches are
        # perfect, keep filtering with additional criteria.
        bestMatch = pickBestMatch(bestMatches)
        return f'The set of atoms is similar to {bestMatch}, but {formatDiffMessage(templatesAtomDiffs[bestMatch], formatAtomDiff)}.{terminalMessage}'
    bestMatches, bestScore = bestMatchAtom(templatesAtomDiffs, bestMatches, False)
    if bestMatches and bestScore[1]:
        bestMatch = pickBestMatch(bestMatches)
        bestMatchDiffs = templatesAtomDiffs[bestMatch]
        # Give additional help in the special cases where the residue is missing
        # sites or hydrogen atoms and nothing else.
        if bestMatchDiffs[0] < 0 and all(diff == 0 for key, diff in bestMatchDiffs.items() if key != 0):
            specialLabel = 'it' if bestMatchDiffs[0] == -1 else 'them'
            specialMessage = f'  You may be able to add {specialLabel} with Modeller.addExtraParticles().'
        elif bestMatchDiffs[1] < 0 and all(diff == 0 for key, diff in bestMatchDiffs.items() if key != 1):
            specialLabel = 'it' if bestMatchDiffs[1] == -1 else 'them'
            specialMessage = f'  You may be able to add {specialLabel} with Modeller.addHydrogens().'
        else:
            specialMessage = terminalMessage
        return f'The set of heavy atoms matches {bestMatch}, but the residue {formatDiffMessage(bestMatchDiffs, formatAtomDiff)}.{specialMessage}'

    # We found templates that are atom-for-atom matches to the residue, so now
    # prepare fingerprints of the residue and templates based on their bonds.
    # The compare the residue with templates based on bonds.

    residueIntBondCounts = Counter(makeIntBondSpec(bond) for bond in res.internal_bonds())

    def makeTemplateIntBondDiff(templateName):
        """Prepares a Counter describing the difference between a residue's and a template's bonds."""
        template = forcefield._templates[templateName]
        return counterSubtract(residueIntBondCounts, Counter(makeIntBondSpec((template.atoms[atom1], template.atoms[atom2])) for atom1, atom2 in template.bonds))

    # We only need to prepare data for the templates remaining in bestMatches.

    templatesIntBondDiffs = {templateName: makeTemplateIntBondDiff(templateName) for templateName in bestMatches}

    bestMatches, bestScore = bestMatchBond(templatesIntBondDiffs, bestMatches)
    if bestMatches and bestScore[1]:
        bestMatch = pickBestMatch(bestMatches)
        if not tuple(res.internal_bonds()):
            # Special message when the residue is missing all internal bonds.
            return (f'The set of atoms matches {bestMatch}, but the residue has no bonds between its atoms.  '
                'If the topology was read from a PDB, it may contain non-standard residues/names and/or be missing CONECT records.')
        return f'The set of atoms matches {bestMatch}, but the residue {formatDiffMessage(templatesIntBondDiffs[bestMatch], formatBondDiff)}.'

    # Finally, check external bonds.  Normally these should always be from heavy
    # atoms, so don't do separate checks excluding vs. including non-heavy atoms
    # (but the check will work regardless).

    residueExtBondCounts = Counter(makeExtBondSpec(bond) for bond in res.external_bonds())

    def makeTemplateExtBondDiff(templateName):
        """Prepares a Counter describing the difference between a residue's and a template's external bonds."""
        template = forcefield._templates[templateName]
        return counterSubtract(residueExtBondCounts, Counter(elemToNum(template.atoms[atom].element) for atom in template.externalBonds))

    templatesExtBondDiffs = {templateName: makeTemplateExtBondDiff(templateName) for templateName in bestMatches}

    bestMatches, bestScore = bestMatchAtom(templatesExtBondDiffs, bestMatches, False)
    if bestMatches and bestScore[1]:
        bestMatch = pickBestMatch(bestMatches)
        bestMatchDiffs = templatesExtBondDiffs[bestMatch]
        if all(value <= 0 for value in bestMatchDiffs.values()):
            # Special message if external bonds are missing only, not different.
            specialMessage = '  Is the chain missing a terminal capping group?'
        else:
            specialMessage = terminalMessage
        return f'The atoms and bonds in the residue match {bestMatch}, but the set of externally bonded atoms {formatDiffMessage(bestMatchDiffs, formatAtomDiff)}.{specialMessage}'

    # If we have matches at this point, atoms and bonds match perfectly, so the
    # connectivity must be different.  If bestMatches is empty, something else
    # went wrong, so return a generic error message.

    if bestMatches:
        # Display all possible matching templates at this point to try to help,
        # since we can't give any more detailed information.
        return f'The atoms and bonds in the residue match {joinMessages(bestMatches)}, but the connectivity is different.'

    return 'This might mean your input topology is missing some atoms or bonds, or possibly that you are using the wrong force field.'

def _createResidueTemplate(residue):
    """Create a _TemplateData template from a Residue object.

    Parameters
    ----------
    residue : Residue
        The Residue from which the template is to be constructed.

    Returns
    -------
    template : _TemplateData
        The residue template, with atom types set to None.

    This method may be useful in creating new residue templates for residues without templates defined by the ForceField.

    """
    template = ForceField._TemplateData(residue.name)
    for atom in residue.atoms():
        template.addAtom(ForceField._TemplateAtomData(atom.name, None, atom.element))
    for (atom1,atom2) in residue.internal_bonds():
        template.addBondByName(atom1.name, atom2.name)
    residue_atoms = [ atom for atom in residue.atoms() ]
    for (atom1,atom2) in residue.external_bonds():
        if atom1 in residue_atoms:
            template.addExternalBondByName(atom1.name)
        elif atom2 in residue_atoms:
            template.addExternalBondByName(atom2.name)
    return template

def _matchImproper(data, torsion, generator):
    type1 = data.atomType[data.atoms[torsion[0]]]
    type2 = data.atomType[data.atoms[torsion[1]]]
    type3 = data.atomType[data.atoms[torsion[2]]]
    type4 = data.atomType[data.atoms[torsion[3]]]
    wildcard = generator.ff._atomClasses['']
    match = None
    for tordef in generator.improper:
        types1 = tordef.types1
        types2 = tordef.types2
        types3 = tordef.types3
        types4 = tordef.types4
        hasWildcard = (wildcard in (types1, types2, types3, types4))
        if match is not None and hasWildcard:
            # Prefer specific definitions over ones with wildcards
            continue
        if type1 in types1:
            for (t2, t3, t4) in itertools.permutations(((type2, 1), (type3, 2), (type4, 3))):
                if t2[0] in types2 and t3[0] in types3 and t4[0] in types4:
                    if tordef.ordering == 'default':
                        # Workaround to be more consistent with AMBER.  It uses wildcards to define most of its
                        # impropers, which leaves the ordering ambiguous.  It then follows some bizarre rules
                        # to pick the order.
                        a1 = torsion[t2[1]]
                        a2 = torsion[t3[1]]
                        e1 = data.atoms[a1].element
                        e2 = data.atoms[a2].element
                        if e1 == e2 and a1 > a2:
                            (a1, a2) = (a2, a1)
                        elif e1 != elem.carbon and (e2 == elem.carbon or e1.mass < e2.mass):
                            (a1, a2) = (a2, a1)
                        match = (a1, a2, torsion[0], torsion[t4[1]], tordef)
                        break
                    elif tordef.ordering == 'charmm':
                        if hasWildcard:
                            # Workaround to be more consistent with AMBER.  It uses wildcards to define most of its
                            # impropers, which leaves the ordering ambiguous.  It then follows some bizarre rules
                            # to pick the order.
                            a1 = torsion[t2[1]]
                            a2 = torsion[t3[1]]
                            e1 = data.atoms[a1].element
                            e2 = data.atoms[a2].element
                            if e1 == e2 and a1 > a2:
                                (a1, a2) = (a2, a1)
                            elif e1 != elem.carbon and (e2 == elem.carbon or e1.mass < e2.mass):
                                (a1, a2) = (a2, a1)
                            match = (a1, a2, torsion[0], torsion[t4[1]], tordef)
                        else:
                            # There are no wildcards, so the order is unambiguous.
                            match = (torsion[0], torsion[t2[1]], torsion[t3[1]], torsion[t4[1]], tordef)
                        break
                    elif tordef.ordering == 'amber':
                        # topology atom indexes
                        a2 = torsion[t2[1]]
                        a3 = torsion[t3[1]]
                        a4 = torsion[t4[1]]
                        # residue indexes
                        r2 = data.atoms[a2].residue.index
                        r3 = data.atoms[a3].residue.index
                        r4 = data.atoms[a4].residue.index
                        # template atom indexes
                        ta2 = data.atomTemplateIndexes[data.atoms[a2]]
                        ta3 = data.atomTemplateIndexes[data.atoms[a3]]
                        ta4 = data.atomTemplateIndexes[data.atoms[a4]]
                        # elements
                        e2 = data.atoms[a2].element
                        e3 = data.atoms[a3].element
                        e4 = data.atoms[a4].element
                        if not hasWildcard:
                            if t2[0] == t4[0] and (r2 > r4 or (r2 == r4 and ta2 > ta4)):
                                (a2, a4) = (a4, a2)
                                r2 = data.atoms[a2].residue.index
                                r4 = data.atoms[a4].residue.index
                                ta2 = data.atomTemplateIndexes[data.atoms[a2]]
                                ta4 = data.atomTemplateIndexes[data.atoms[a4]]
                            if t3[0] == t4[0] and (r3 > r4 or (r3 == r4 and ta3 > ta4)):
                                (a3, a4) = (a4, a3)
                                r3 = data.atoms[a3].residue.index
                                r4 = data.atoms[a4].residue.index
                                ta3 = data.atomTemplateIndexes[data.atoms[a3]]
                                ta4 = data.atomTemplateIndexes[data.atoms[a4]]
                            if t2[0] == t3[0] and (r2 > r3 or (r2 == r3 and ta2 > ta3)):
                                (a2, a3) = (a3, a2)
                        else:
                            if e2 == e4 and (r2 > r4 or (r2 == r4 and ta2 > ta4)):
                                (a2, a4) = (a4, a2)
                                r2 = data.atoms[a2].residue.index
                                r4 = data.atoms[a4].residue.index
                                ta2 = data.atomTemplateIndexes[data.atoms[a2]]
                                ta4 = data.atomTemplateIndexes[data.atoms[a4]]
                            if e3 == e4 and (r3 > r4 or (r3 == r4 and ta3 > ta4)):
                                (a3, a4) = (a4, a3)
                                r3 = data.atoms[a3].residue.index
                                r4 = data.atoms[a4].residue.index
                                ta3 = data.atomTemplateIndexes[data.atoms[a3]]
                                ta4 = data.atomTemplateIndexes[data.atoms[a4]]
                            if r2 > r3 or (r2 == r3 and ta2 > ta3):
                                (a2, a3) = (a3, a2)
                        match = (a2, a3, torsion[0], a4, tordef)
                        break
                    elif tordef.ordering == 'smirnoff':
                        # topology atom indexes
                        a1 = torsion[0]
                        a2 = torsion[t2[1]]
                        a3 = torsion[t3[1]]
                        a4 = torsion[t4[1]]
                        # enforce exact match
                        match = (a1, a2, a3, a4, tordef)
                        break

    return match


# The following classes are generators that know how to create Force subclasses and add them to a System that is being
# created.  Each generator class must define two methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding generator object; 2) a createForce() method that constructs the Force object and adds it
# to the System.  The static method should be added to the parsers map.
## @private
class HarmonicBondGenerator(object):
    """A HarmonicBondGenerator constructs a HarmonicBondForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.bondsForAtomType = defaultdict(set)
        self.types1 = []
        self.types2 = []
        self.length = []
        self.k = []

    def registerBond(self, parameters):
        types = self.ff._findAtomTypes(parameters, 2)
        if None not in types:
            index = len(self.types1)
            self.types1.append(types[0])
            self.types2.append(types[1])
            for t in types[0]:
                self.bondsForAtomType[t].add(index)
            for t in types[1]:
                self.bondsForAtomType[t].add(index)
            self.length.append(_convertParameterToNumber(parameters['length']))
            self.k.append(_convertParameterToNumber(parameters['k']))

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, HarmonicBondGenerator)]
        if len(existing) == 0:
            generator = HarmonicBondGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]
        for bond in element.findall('Bond'):
            generator.registerBond(bond.attrib)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [f for f in sys.getForces() if type(f) == mm.HarmonicBondForce]
        if len(existing) == 0:
            force = mm.HarmonicBondForce()
            sys.addForce(force)
        else:
            force = existing[0]
        for bond in data.bonds:
            type1 = data.atomType[data.atoms[bond.atom1]]
            type2 = data.atomType[data.atoms[bond.atom2]]
            for i in self.bondsForAtomType[type1]:
                types1 = self.types1[i]
                types2 = self.types2[i]
                if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
                    bond.length = self.length[i]
                    if bond.isConstrained:
                        data.addConstraint(sys, bond.atom1, bond.atom2, self.length[i])
                    if self.k[i] != 0:
                        # flexibleConstraints allows us to add parameters even if the DOF is
                        # constrained
                        if not bond.isConstrained or args.get('flexibleConstraints', False):
                            force.addBond(bond.atom1, bond.atom2, self.length[i], self.k[i])
                    break

parsers["HarmonicBondForce"] = HarmonicBondGenerator.parseElement

## @private
class HarmonicAngleGenerator(object):
    """A HarmonicAngleGenerator constructs a HarmonicAngleForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.anglesForAtom2Type = defaultdict(list)
        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.angle = []
        self.k = []

    def registerAngle(self, parameters):
        types = self.ff._findAtomTypes(parameters, 3)
        if None not in types:
            index = len(self.types1)
            self.types1.append(types[0])
            self.types2.append(types[1])
            self.types3.append(types[2])
            for t in types[1]:
                self.anglesForAtom2Type[t].append(index)
            self.angle.append(_convertParameterToNumber(parameters['angle']))
            self.k.append(_convertParameterToNumber(parameters['k']))

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, HarmonicAngleGenerator)]
        if len(existing) == 0:
            generator = HarmonicAngleGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]
        for angle in element.findall('Angle'):
            generator.registerAngle(angle.attrib)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        pass

    def postprocessSystem(self, sys, data, args):
        # We need to wait until after all bonds have been added so their lengths will be set correctly.

        existing = [f for f in sys.getForces() if type(f) == mm.HarmonicAngleForce]
        if len(existing) == 0:
            force = mm.HarmonicAngleForce()
            sys.addForce(force)
        else:
            force = existing[0]
        for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):
            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]
            for i in self.anglesForAtom2Type[type2]:
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
                                data.addConstraint(sys, angle[0], angle[2], length)
                    if self.k[i] != 0:
                        if not isConstrained or args.get('flexibleConstraints', False):
                            force.addAngle(angle[0], angle[1], angle[2], self.angle[i], self.k[i])
                    break

parsers["HarmonicAngleForce"] = HarmonicAngleGenerator.parseElement

## @private
class PeriodicTorsion(object):
    """A PeriodicTorsion records the information for a periodic torsion definition."""

    def __init__(self, types):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.periodicity = []
        self.phase = []
        self.k = []
        self.ordering = 'default'

## @private
class PeriodicTorsionGenerator(object):
    """A PeriodicTorsionGenerator constructs a PeriodicTorsionForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.proper = []
        self.improper = []
        self.propersForAtomType = defaultdict(set)

    def registerProperTorsion(self, parameters):
        torsion = self.ff._parseTorsion(parameters)
        if torsion is not None:
            index = len(self.proper)
            self.proper.append(torsion)
            for t in torsion.types2:
                self.propersForAtomType[t].add(index)
            for t in torsion.types3:
                self.propersForAtomType[t].add(index)

    def registerImproperTorsion(self, parameters, ordering='default'):
        torsion = self.ff._parseTorsion(parameters)
        if torsion is not None:
            if ordering in ['default', 'charmm', 'amber', 'smirnoff']:
                torsion.ordering = ordering
            else:
                raise ValueError('Illegal ordering type %s for improper torsion %s' % (ordering, torsion))
            self.improper.append(torsion)

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, PeriodicTorsionGenerator)]
        if len(existing) == 0:
            generator = PeriodicTorsionGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]
        for torsion in element.findall('Proper'):
            generator.registerProperTorsion(torsion.attrib)
        for torsion in element.findall('Improper'):
            if 'ordering' in element.attrib:
                generator.registerImproperTorsion(torsion.attrib, element.attrib['ordering'])
            else:
                generator.registerImproperTorsion(torsion.attrib)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [f for f in sys.getForces() if type(f) == mm.PeriodicTorsionForce]
        if len(existing) == 0:
            force = mm.PeriodicTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]
        wildcard = self.ff._atomClasses['']
        proper_cache = {}
        for torsion in data.propers:
            type1, type2, type3, type4 = [data.atomType[data.atoms[torsion[i]]] for i in range(4)]
            sig = (type1, type2, type3, type4)
            sig = frozenset((sig, sig[::-1]))
            match = proper_cache.get(sig, None)
            if match == -1:
                continue
            if match is None:
                for index in self.propersForAtomType[type2]:
                    tordef = self.proper[index]
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
                if match is None:
                    proper_cache[sig] = -1
                else:
                    proper_cache[sig] = match
            if match is not None:
                for i in range(len(match.phase)):
                    if match.k[i] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], match.periodicity[i], match.phase[i], match.k[i])
        impr_cache = {}
        for torsion in data.impropers:
            t1, t2, t3, t4 = [data.atomType[data.atoms[torsion[i]]] for i in range(4)]
            sig = (t1, t2, t3, t4)
            match = impr_cache.get(sig, None)
            if match == -1:
                # Previously checked, and doesn't appear in the database
                continue
            elif match:
                i1, i2, i3, i4, tordef = match
                a1, a2, a3, a4 = (torsion[i] for i in (i1, i2, i3, i4))
                match = (a1, a2, a3, a4, tordef)
            if match is None:
                match = _matchImproper(data, torsion, self)
                if match is not None:
                    order = match[:4]
                    i1, i2, i3, i4 = tuple(torsion.index(a) for a in order)
                    impr_cache[sig] = (i1, i2, i3, i4, match[-1])
                else:
                    impr_cache[sig] = -1
            if match is not None:
                (a1, a2, a3, a4, tordef) = match
                for i in range(len(tordef.phase)):
                    if tordef.k[i] != 0:
                        if tordef.ordering == 'smirnoff':
                            # Add all torsions in trefoil
                            force.addTorsion(a1, a2, a3, a4, tordef.periodicity[i], tordef.phase[i], tordef.k[i])
                            force.addTorsion(a1, a3, a4, a2, tordef.periodicity[i], tordef.phase[i], tordef.k[i])
                            force.addTorsion(a1, a4, a2, a3, tordef.periodicity[i], tordef.phase[i], tordef.k[i])
                        else:
                            force.addTorsion(a1, a2, a3, a4, tordef.periodicity[i], tordef.phase[i], tordef.k[i])
parsers["PeriodicTorsionForce"] = PeriodicTorsionGenerator.parseElement

## @private
class RBTorsion(object):
    """An RBTorsion records the information for a Ryckaert-Bellemans torsion definition."""

    def __init__(self, types, c, ordering='charmm'):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.c = c
        if ordering in ['default', 'charmm', 'amber']:
            self.ordering = ordering
        else:
            raise ValueError('Illegal ordering type %s for RBTorsion (%s,%s,%s,%s)' % (ordering, types[0], types[1], types[2], types[3]))

## @private
class RBTorsionGenerator(object):
    """An RBTorsionGenerator constructs an RBTorsionForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.proper = []
        self.improper = []

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, RBTorsionGenerator)]
        if len(existing) == 0:
            generator = RBTorsionGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]
        for torsion in element.findall('Proper'):
            types = ff._findAtomTypes(torsion.attrib, 4)
            if None not in types:
                generator.proper.append(RBTorsion(types, [float(torsion.attrib['c'+str(i)]) for i in range(6)]))
        for torsion in element.findall('Improper'):
            types = ff._findAtomTypes(torsion.attrib, 4)
            if None not in types:
                if 'ordering' in element.attrib:
                    generator.improper.append(RBTorsion(types, [float(torsion.attrib['c'+str(i)]) for i in range(6)], element.attrib['ordering']))
                else:
                    generator.improper.append(RBTorsion(types, [float(torsion.attrib['c'+str(i)]) for i in range(6)]))

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [f for f in sys.getForces() if type(f) == mm.RBTorsionForce]
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
            match = _matchImproper(data, torsion, self)
            if match is not None:
                (a1, a2, a3, a4, tordef) = match
                force.addTorsion(a1, a2, a3, a4, tordef.c[0], tordef.c[1], tordef.c[2], tordef.c[3], tordef.c[4], tordef.c[5])

parsers["RBTorsionForce"] = RBTorsionGenerator.parseElement

## @private
class CMAPTorsion(object):
    """A CMAPTorsion records the information for a CMAP torsion definition."""

    def __init__(self, types, map):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.types5 = types[4]
        self.map = map

## @private
class CMAPTorsionGenerator(object):
    """A CMAPTorsionGenerator constructs a CMAPTorsionForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.torsions = []
        self.maps = []

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, CMAPTorsionGenerator)]
        if len(existing) == 0:
            generator = CMAPTorsionGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]
        mapOffset = len(generator.maps)
        for map in element.findall('Map'):
            values = [float(x) for x in map.text.split()]
            size = sqrt(len(values))
            if size*size != len(values):
                raise ValueError('CMAP must have the same number of elements along each dimension')
            generator.maps.append(values)
        for torsion in element.findall('Torsion'):
            types = ff._findAtomTypes(torsion.attrib, 5)
            if None not in types:
                generator.torsions.append(CMAPTorsion(types, int(torsion.attrib['map']) + mapOffset))

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [f for f in sys.getForces() if type(f) == mm.CMAPTorsionForce]
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
class NonbondedGenerator(object):
    """A NonbondedGenerator constructs a NonbondedForce."""

    SCALETOL = 1e-5

    def __init__(self, forcefield, coulomb14scale, lj14scale, useDispersionCorrection):
        self.ff = forcefield
        self.coulomb14scale = coulomb14scale
        self.lj14scale = lj14scale
        self.useDispersionCorrection = useDispersionCorrection
        self.params = ForceField._AtomTypeParameters(forcefield, 'NonbondedForce', 'Atom', ('charge', 'sigma', 'epsilon'))

    def registerAtom(self, parameters):
        self.params.registerAtom(parameters)

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, NonbondedGenerator)]
        if 'useDispersionCorrection' in element.attrib:
            useDispersionCorrection = bool(eval(element.attrib.get('useDispersionCorrection')))
        else:
            useDispersionCorrection = None
        if len(existing) == 0:
            generator = NonbondedGenerator(
                ff,
                float(element.attrib['coulomb14scale']),
                float(element.attrib['lj14scale']),
                useDispersionCorrection
            )
            ff.registerGenerator(generator)
        else:
            # Multiple <NonbondedForce> tags were found, probably in different files.  Simply add more types to the existing one.
            generator = existing[0]
            if (abs(generator.coulomb14scale - float(element.attrib['coulomb14scale'])) > NonbondedGenerator.SCALETOL
                or abs(generator.lj14scale - float(element.attrib['lj14scale'])) > NonbondedGenerator.SCALETOL
            ):
                raise ValueError('Found multiple NonbondedForce tags with different 1-4 scales')
            if (
                    generator.useDispersionCorrection is not None
                    and useDispersionCorrection is not None
                    and generator.useDispersionCorrection != useDispersionCorrection
            ):
                raise ValueError('Found multiple NonbondedForce tags with different useDispersionCorrection settings.')
        generator.params.parseDefinitions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic,
                     Ewald:mm.NonbondedForce.Ewald,
                     PME:mm.NonbondedForce.PME,
                     LJPME:mm.NonbondedForce.LJPME}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for NonbondedForce')
        force = mm.NonbondedForce()
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            force.addParticle(values[0], values[1], values[2])
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        if args['switchDistance'] is not None:
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(args['switchDistance'])
        if 'ewaldErrorTolerance' in args:
            force.setEwaldErrorTolerance(args['ewaldErrorTolerance'])
        if 'useDispersionCorrection' in args:
            lrc_keyword = bool(args['useDispersionCorrection'])
            if self.useDispersionCorrection is not None and lrc_keyword != self.useDispersionCorrection:
                warnings.warn(
                    "Conflicting settings for useDispersionCorrection in createSystem() and forcefield file. "
                    "Using the one specified in createSystem()."
                )
            force.setUseDispersionCorrection(lrc_keyword)
        elif self.useDispersionCorrection is not None:
            force.setUseDispersionCorrection(self.useDispersionCorrection)
        else:
            # by default
            force.setUseDispersionCorrection(True)
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # Create the exceptions.

        bondIndices = _findBondsForExclusions(data, sys)
        nonbonded = [f for f in sys.getForces() if isinstance(f, mm.NonbondedForce)][0]
        nonbonded.createExceptionsFromBonds(bondIndices, self.coulomb14scale, self.lj14scale)

parsers["NonbondedForce"] = NonbondedGenerator.parseElement

## @private
class LennardJonesGenerator(object):
    """A NBFix generator to construct the L-J force with NBFIX implemented as a lookup table"""

    def __init__(self, forcefield, lj14scale, useDispersionCorrection):
        self.ff = forcefield
        self.nbfixParameters = []
        self.nbfixTypes1 = defaultdict(set)
        self.nbfixTypes2 = defaultdict(set)
        self.lj14scale = lj14scale
        self.useDispersionCorrection = useDispersionCorrection
        self.ljTypes = ForceField._AtomTypeParameters(forcefield, 'LennardJonesForce', 'Atom', ('sigma', 'epsilon'))

    def registerNBFIX(self, parameters):
        types = self.ff._findAtomTypes(parameters, 2)

        if None not in types:
            sigma = _convertParameterToNumber(parameters['sigma'])
            epsilon = _convertParameterToNumber(parameters['epsilon'])

            # Retrieve the index of nbfixParameters into which this sigma and
            # epsilon will be stored, then register this index with the atom
            # types that should have this sigma and epsilon applied.
            nbfixIndex = len(self.nbfixParameters)
            self.nbfixParameters.append([sigma, epsilon])
            for type1 in types[0]:
                self.nbfixTypes1[type1].add(nbfixIndex)
            for type2 in types[1]:
                self.nbfixTypes2[type2].add(nbfixIndex)

    def getNBFIX(self, type1, type2):
        nbfixIndices = (self.nbfixTypes1[type1] & self.nbfixTypes2[type2]) | (self.nbfixTypes2[type1] & self.nbfixTypes1[type2])
        if nbfixIndices:
            if len(nbfixIndices) > 1:
                raise ValueError('Multiple NBFixPair entries match atom types %s-%s.' % (type1, type2))
            return self.nbfixParameters[nbfixIndices.pop()]
        else:
            return None

    def registerLennardJones(self, parameters):
        self.ljTypes.registerAtom(parameters)

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, LennardJonesGenerator)]
        if 'useDispersionCorrection' in element.attrib:
            useDispersionCorrection = bool(eval(element.attrib.get('useDispersionCorrection')))
        else:
            useDispersionCorrection = None
        if len(existing) == 0:
            generator = LennardJonesGenerator(
                ff,
                float(element.attrib['lj14scale']),
                useDispersionCorrection=useDispersionCorrection
            )
            ff.registerGenerator(generator)
        else:
            # Multiple <LennardJonesForce> tags were found, probably in different files
            generator = existing[0]
            if abs(generator.lj14scale - float(element.attrib['lj14scale'])) > NonbondedGenerator.SCALETOL:
                raise ValueError('Found multiple LennardJonesForce tags with different 1-4 scales')
            if (
                    generator.useDispersionCorrection is not None
                    and useDispersionCorrection is not None
                    and generator.useDispersionCorrection != useDispersionCorrection
            ):
                raise ValueError('Found multiple LennardJonesForce tags with different useDispersionCorrection settings.')
        for LJ in element.findall('Atom'):
            generator.registerLennardJones(LJ.attrib)
        for Nbfix in element.findall('NBFixPair'):
            generator.registerNBFIX(Nbfix.attrib)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        # First derive the lookup tables.  We need to include entries for every type
        # that a) appears in the system and b) has unique parameters.

        nbfixTypeSet = {t for nbfixTypes in (self.nbfixTypes1, self.nbfixTypes2) for t in nbfixTypes if nbfixTypes[t]}
        allTypes = set(data.atomType[atom] for atom in data.atoms)
        mergedTypes = []
        mergedTypeParams = []
        paramsToMergedType = {}
        typeToMergedType = {}
        for t in allTypes:
            typeParams = self.ljTypes.paramsForType[t]
            params = (typeParams['sigma'], typeParams['epsilon'])
            if t in nbfixTypeSet:
                # NBFIX types cannot be merged.
                typeToMergedType[t] = len(mergedTypes)
                mergedTypes.append(t)
                mergedTypeParams.append(params)
            elif params in paramsToMergedType:
                # We can merge this with another type.
                typeToMergedType[t] = paramsToMergedType[params]
            else:
                # This is a new type.
                typeToMergedType[t] = len(mergedTypes)
                paramsToMergedType[params] = len(mergedTypes)
                mergedTypes.append(t)
                mergedTypeParams.append(params)

        # Now everything is assigned. Create the A- and B-coefficient arrays

        numLjTypes = len(mergedTypes)
        acoef = [0]*(numLjTypes*numLjTypes)
        bcoef = acoef[:]
        for m in range(numLjTypes):
            for n in range(numLjTypes):
                nbfix = self.getNBFIX(mergedTypes[m], mergedTypes[n])
                if nbfix is not None:
                    sigma, epsilon = nbfix
                    sigma6 = sigma**6
                    acoef[m+numLjTypes*n] = 4*epsilon*sigma6*sigma6
                    bcoef[m+numLjTypes*n] = 4*epsilon*sigma6
                    continue
                else:
                    sigma = 0.5*(mergedTypeParams[m][0]+mergedTypeParams[n][0])
                    sigma6 = sigma**6
                    epsilon = math.sqrt(mergedTypeParams[m][1]*mergedTypeParams[n][1])
                    acoef[m+numLjTypes*n] = 4*epsilon*sigma6*sigma6
                    bcoef[m+numLjTypes*n] = 4*epsilon*sigma6

        self.force = mm.CustomNonbondedForce('acoef(type1, type2)/r^12 - bcoef(type1, type2)/r^6;')
        self.force.addTabulatedFunction('acoef', mm.Discrete2DFunction(numLjTypes, numLjTypes, acoef))
        self.force.addTabulatedFunction('bcoef', mm.Discrete2DFunction(numLjTypes, numLjTypes, bcoef))
        self.force.addPerParticleParameter('type')
        self.force.setName('LennardJones')
        if nonbondedMethod in [CutoffPeriodic, Ewald, PME]:
            self.force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        elif nonbondedMethod is NoCutoff:
            self.force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is CutoffNonPeriodic:
            self.force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        elif nonbondedMethod is LJPME:
            raise ValueError('LJPME is not supported by LennardJonesForce')
        else:
            raise AssertionError('Unrecognized nonbonded method [%s]' % nonbondedMethod)
        if args['switchDistance'] is not None:
            self.force.setUseSwitchingFunction(True)
            self.force.setSwitchingDistance(args['switchDistance'])
        if 'useDispersionCorrection' in args:
            lrc_keyword = bool(args['useDispersionCorrection'])
            if self.useDispersionCorrection is not None and lrc_keyword != self.useDispersionCorrection:
                warnings.warn(
                    "Conflicting settings for useDispersionCorrection in createSystem() and forcefield file. "
                    "Using the one specified in createSystem()."
                )
            self.force.setUseLongRangeCorrection(lrc_keyword)
        elif self.useDispersionCorrection is not None:
            self.force.setUseLongRangeCorrection(self.useDispersionCorrection)
        else:
            # by default
            self.force.setUseLongRangeCorrection(True)

        # Add the particles

        for atom in data.atoms:
            self.force.addParticle((typeToMergedType[data.atomType[atom]],))
        self.force.setCutoffDistance(nonbondedCutoff)
        sys.addForce(self.force)

    def postprocessSystem(self, sys, data, args):
        # Create the exceptions.

        bondIndices = _findBondsForExclusions(data, sys)
        forceCopy = deepcopy(self.force)
        forceCopy.createExclusionsFromBonds(bondIndices, 2)
        self.force.createExclusionsFromBonds(bondIndices, 3)
        if self.force.getNumExclusions() > forceCopy.getNumExclusions() and self.lj14scale != 0:
            # We need to create a CustomBondForce and use it to implement the scaled 1-4 interactions.

            bonded = mm.CustomBondForce('%g*epsilon*((sigma/r)^12-(sigma/r)^6)' % (4*self.lj14scale))
            bonded.addPerBondParameter('sigma')
            bonded.addPerBondParameter('epsilon')
            bonded.setName('LennardJones14')
            sys.addForce(bonded)
            skip = set(tuple(forceCopy.getExclusionParticles(i)) for i in range(forceCopy.getNumExclusions()))
            for i in range(self.force.getNumExclusions()):
                p1,p2 = self.force.getExclusionParticles(i)
                a1 = data.atoms[p1]
                a2 = data.atoms[p2]
                if (p1,p2) not in skip and (p2,p1) not in skip:
                    nbfix = self.getNBFIX(data.atomType[a1], data.atomType[a2])
                    if nbfix is not None:
                        sigma, epsilon = nbfix
                    else:
                        values1 = self.ljTypes.getAtomParameters(a1, data)
                        values2 = self.ljTypes.getAtomParameters(a2, data)
                        extra1 = self.ljTypes.getExtraParameters(a1, data)
                        extra2 = self.ljTypes.getExtraParameters(a2, data)
                        sigma1 = float(extra1['sigma14']) if 'sigma14' in extra1 else values1[0]
                        sigma2 = float(extra2['sigma14']) if 'sigma14' in extra2 else values2[0]
                        epsilon1 = float(extra1['epsilon14']) if 'epsilon14' in extra1 else values1[1]
                        epsilon2 = float(extra2['epsilon14']) if 'epsilon14' in extra2 else values2[1]
                        sigma = 0.5*(sigma1+sigma2)
                        epsilon = sqrt(epsilon1*epsilon2)
                    bonded.addBond(p1, p2, (sigma, epsilon))

parsers["LennardJonesForce"] = LennardJonesGenerator.parseElement

## @private
class GBSAOBCGenerator(object):
    """A GBSAOBCGenerator constructs a GBSAOBCForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.params = ForceField._AtomTypeParameters(forcefield, 'GBSAOBCForce', 'Atom', ('charge', 'radius', 'scale'))

    def registerAtom(self, parameters):
        self.params.registerAtom(parameters)

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, GBSAOBCGenerator)]
        if len(existing) == 0:
            generator = GBSAOBCGenerator(ff)
            ff.registerGenerator(generator)
        else:
            # Multiple <GBSAOBCForce> tags were found, probably in different files.  Simply add more types to the existing one.
            generator = existing[0]
        generator.params.parseDefinitions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.GBSAOBCForce.NoCutoff,
                     CutoffNonPeriodic:mm.GBSAOBCForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.GBSAOBCForce.CutoffPeriodic,
                     Ewald:mm.GBSAOBCForce.CutoffPeriodic,
                     PME:mm.GBSAOBCForce.CutoffPeriodic,
                     LJPME:mm.GBSAOBCForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for GBSAOBCForce')
        force = mm.GBSAOBCForce()
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            force.addParticle(values[0], values[1], values[2])
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        if 'soluteDielectric' in args:
            force.setSoluteDielectric(float(args['soluteDielectric']))
        if 'solventDielectric' in args:
            force.setSolventDielectric(float(args['solventDielectric']))
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # Disable the reaction field approximation, since it produces bad results when combined with GB.

        for force in sys.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setReactionFieldDielectric(1.0)

parsers["GBSAOBCForce"] = GBSAOBCGenerator.parseElement

## @private
class CustomBondGenerator(object):
    """A CustomBondGenerator constructs a CustomBondForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.types1 = []
        self.types2 = []
        self.globalParams = {}
        self.perBondParams = []
        self.paramValues = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomBondGenerator(ff)
        ff.registerGenerator(generator)
        generator.energy = element.attrib['energy']
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerBondParameter'):
            generator.perBondParams.append(param.attrib['name'])
        for bond in element.findall('Bond'):
            types = ff._findAtomTypes(bond.attrib, 2)
            if None not in types:
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
class CustomAngleGenerator(object):
    """A CustomAngleGenerator constructs a CustomAngleForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.globalParams = {}
        self.perAngleParams = []
        self.paramValues = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomAngleGenerator(ff)
        ff.registerGenerator(generator)
        generator.energy = element.attrib['energy']
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerAngleParameter'):
            generator.perAngleParams.append(param.attrib['name'])
        for angle in element.findall('Angle'):
            types = ff._findAtomTypes(angle.attrib, 3)
            if None not in types:
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
class CustomTorsion(object):
    """A CustomTorsion records the information for a custom torsion definition."""

    def __init__(self, types, paramValues, ordering='charmm'):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.paramValues = paramValues
        if ordering in ['default', 'charmm', 'amber']:
            self.ordering = ordering
        else:
            raise ValueError('Illegal ordering type %s for CustomTorsion (%s,%s,%s,%s)' % (ordering, types[0], types[1], types[2], types[3]))

## @private
class CustomTorsionGenerator(object):
    """A CustomTorsionGenerator constructs a CustomTorsionForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.proper = []
        self.improper = []
        self.globalParams = {}
        self.perTorsionParams = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomTorsionGenerator(ff)
        ff.registerGenerator(generator)
        generator.energy = element.attrib['energy']
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerTorsionParameter'):
            generator.perTorsionParams.append(param.attrib['name'])
        for torsion in element.findall('Proper'):
            types = ff._findAtomTypes(torsion.attrib, 4)
            if None not in types:
                generator.proper.append(CustomTorsion(types, [float(torsion.attrib[param]) for param in generator.perTorsionParams]))
        for torsion in element.findall('Improper'):
            types = ff._findAtomTypes(torsion.attrib, 4)
            if None not in types:
                if 'ordering' in element.attrib:
                    generator.improper.append(CustomTorsion(types, [float(torsion.attrib[param]) for param in generator.perTorsionParams], element.attrib['ordering']))
                else:
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
            match = _matchImproper(data, torsion, self)
            if match is not None:
                (a1, a2, a3, a4, tordef) = match
                force.addTorsion(a1, a2, a3, a4, tordef.paramValues)

parsers["CustomTorsionForce"] = CustomTorsionGenerator.parseElement

## @private
class CustomNonbondedGenerator(object):
    """A CustomNonbondedGenerator constructs a CustomNonbondedForce."""

    def __init__(self, forcefield, energy, bondCutoff):
        self.ff = forcefield
        self.energy = energy
        self.bondCutoff = bondCutoff
        self.globalParams = {}
        self.perParticleParams = []
        self.computedValues = {}
        self.functions = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomNonbondedGenerator(ff, element.attrib['energy'], int(element.attrib['bondCutoff']))
        ff.registerGenerator(generator)
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerParticleParameter'):
            generator.perParticleParams.append(param.attrib['name'])
        for value in element.findall('ComputedValue'):
            generator.computedValues[value.attrib['name']] = value.attrib['expression']
        generator.params = ForceField._AtomTypeParameters(ff, 'CustomNonbondedForce', 'Atom', generator.perParticleParams)
        generator.params.parseDefinitions(element)
        generator.functions += _parseFunctions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.CustomNonbondedForce.NoCutoff,
                     CutoffNonPeriodic:mm.CustomNonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.CustomNonbondedForce.CutoffPeriodic,
                     Ewald:mm.CustomNonbondedForce.CutoffPeriodic,
                     PME:mm.CustomNonbondedForce.CutoffPeriodic,
                     LJPME:mm.CustomNonbondedForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for CustomNonbondedForce')
        force = mm.CustomNonbondedForce(self.energy)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perParticleParams:
            force.addPerParticleParameter(param)
        for name in self.computedValues:
            force.addComputedValue(name, self.computedValues[name])
        _createFunctions(force, self.functions)
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            force.addParticle(values)
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # Create the exclusions.

        bondIndices = _findBondsForExclusions(data, sys)
        for f in sys.getForces():
            if isinstance(f, mm.CustomNonbondedForce) and f.getEnergyFunction() == self.energy:
                f.createExclusionsFromBonds(bondIndices, self.bondCutoff)

parsers["CustomNonbondedForce"] = CustomNonbondedGenerator.parseElement

## @private
class CustomGBGenerator(object):
    """A CustomGBGenerator constructs a CustomGBForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.globalParams = {}
        self.perParticleParams = []
        self.computedValues = []
        self.energyTerms = []
        self.functions = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomGBGenerator(ff)
        ff.registerGenerator(generator)
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerParticleParameter'):
            generator.perParticleParams.append(param.attrib['name'])
        generator.params = ForceField._AtomTypeParameters(ff, 'CustomGBForce', 'Atom', generator.perParticleParams)
        generator.params.parseDefinitions(element)
        computationMap = {"SingleParticle" : mm.CustomGBForce.SingleParticle,
                          "ParticlePair" : mm.CustomGBForce.ParticlePair,
                          "ParticlePairNoExclusions" : mm.CustomGBForce.ParticlePairNoExclusions}
        for value in element.findall('ComputedValue'):
            generator.computedValues.append((value.attrib['name'], value.text, computationMap[value.attrib['type']]))
        for term in element.findall('EnergyTerm'):
            generator.energyTerms.append((term.text, computationMap[term.attrib['type']]))
        generator.functions += _parseFunctions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.CustomGBForce.NoCutoff,
                     CutoffNonPeriodic:mm.CustomGBForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.CustomGBForce.CutoffPeriodic,
                     Ewald:mm.CustomGBForce.CutoffPeriodic,
                     PME:mm.CustomGBForce.CutoffPeriodic,
                     LJPME:mm.CustomGBForce.CutoffPeriodic}
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
        _createFunctions(force, self.functions)
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            force.addParticle(values)
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        sys.addForce(force)

parsers["CustomGBForce"] = CustomGBGenerator.parseElement

## @private
class CustomHbondGenerator(object):
    """A CustomHbondGenerator constructs a CustomHbondForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.donorTypes1 = []
        self.donorTypes2 = []
        self.donorTypes3 = []
        self.acceptorTypes1 = []
        self.acceptorTypes2 = []
        self.acceptorTypes3 = []
        self.globalParams = {}
        self.perDonorParams = []
        self.perAcceptorParams = []
        self.donorParamValues = []
        self.acceptorParamValues = []
        self.functions = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomHbondGenerator(ff)
        ff.registerGenerator(generator)
        generator.energy = element.attrib['energy']
        generator.bondCutoff = int(element.attrib['bondCutoff'])
        generator.particlesPerDonor = int(element.attrib['particlesPerDonor'])
        generator.particlesPerAcceptor = int(element.attrib['particlesPerAcceptor'])
        if generator.particlesPerDonor < 1 or generator.particlesPerDonor > 3:
            raise ValueError('Illegal value for particlesPerDonor for CustomHbondForce')
        if generator.particlesPerAcceptor < 1 or generator.particlesPerAcceptor > 3:
            raise ValueError('Illegal value for particlesPerAcceptor for CustomHbondForce')
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerDonorParameter'):
            generator.perDonorParams.append(param.attrib['name'])
        for param in element.findall('PerAcceptorParameter'):
            generator.perAcceptorParams.append(param.attrib['name'])
        for donor in element.findall('Donor'):
            types = ff._findAtomTypes(donor.attrib, 3)[:generator.particlesPerDonor]
            if None not in types:
                generator.donorTypes1.append(types[0])
                if len(types) > 1:
                    generator.donorTypes2.append(types[1])
                if len(types) > 2:
                    generator.donorTypes3.append(types[2])
                generator.donorParamValues.append([float(donor.attrib[param]) for param in generator.perDonorParams])
        for acceptor in element.findall('Acceptor'):
            types = ff._findAtomTypes(acceptor.attrib, 3)[:generator.particlesPerAcceptor]
            if None not in types:
                generator.acceptorTypes1.append(types[0])
                if len(types) > 1:
                    generator.acceptorTypes2.append(types[1])
                if len(types) > 2:
                    generator.acceptorTypes3.append(types[2])
                generator.acceptorParamValues.append([float(acceptor.attrib[param]) for param in generator.perAcceptorParams])
        generator.functions += _parseFunctions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.CustomHbondForce.NoCutoff,
                     CutoffNonPeriodic:mm.CustomHbondForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.CustomHbondForce.CutoffPeriodic,
                     Ewald:mm.CustomHbondForce.CutoffPeriodic,
                     PME:mm.CustomHbondForce.CutoffPeriodic,
                     LJPME:mm.CustomHbondForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for CustomNonbondedForce')
        force = mm.CustomHbondForce(self.energy)
        sys.addForce(force)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perDonorParams:
            force.addPerDonorParameter(param)
        for param in self.perAcceptorParams:
            force.addPerAcceptorParameter(param)
        _createFunctions(force, self.functions)
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)

        # Add donors.

        if self.particlesPerDonor == 1:
            for atom in data.atoms:
                type1 = data.atomType[atom]
                for i in range(len(self.donorTypes1)):
                    types1 = self.donorTypes1[i]
                    if type1 in types1:
                        force.addDonor(atom.index, -1, -1, self.donorParamValues[i])
        elif self.particlesPerDonor == 2:
            for bond in data.bonds:
                type1 = data.atomType[data.atoms[bond.atom1]]
                type2 = data.atomType[data.atoms[bond.atom2]]
                for i in range(len(self.donorTypes1)):
                    types1 = self.donorTypes1[i]
                    types2 = self.donorTypes2[i]
                    if type1 in types1 and type2 in types2:
                        force.addDonor(bond.atom1, bond.atom2, -1, self.donorParamValues[i])
                    elif type1 in types2 and type2 in types1:
                        force.addDonor(bond.atom2, bond.atom1, -1, self.donorParamValues[i])
        else:
            for angle in data.angles:
                type1 = data.atomType[data.atoms[angle[0]]]
                type2 = data.atomType[data.atoms[angle[1]]]
                type3 = data.atomType[data.atoms[angle[2]]]
                for i in range(len(self.donorTypes1)):
                    types1 = self.donorTypes1[i]
                    types2 = self.donorTypes2[i]
                    types3 = self.donorTypes3[i]
                    if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                        force.addDonor(angle[0], angle[1], angle[2], self.donorParamValues[i])

        # Add acceptors.

        if self.particlesPerAcceptor == 1:
            for atom in data.atoms:
                type1 = data.atomType[atom]
                for i in range(len(self.acceptorTypes1)):
                    types1 = self.acceptorTypes1[i]
                    if type1 in types1:
                        force.addAcceptor(atom.index, -1, -1, self.acceptorParamValues[i])
        elif self.particlesPerAcceptor == 2:
            for bond in data.bonds:
                type1 = data.atomType[data.atoms[bond.atom1]]
                type2 = data.atomType[data.atoms[bond.atom2]]
                for i in range(len(self.acceptorTypes1)):
                    types1 = self.acceptorTypes1[i]
                    types2 = self.acceptorTypes2[i]
                    if type1 in types1 and type2 in types2:
                        force.addAcceptor(bond.atom1, bond.atom2, -1, self.acceptorParamValues[i])
                    elif type1 in types2 and type2 in types1:
                        force.addAcceptor(bond.atom2, bond.atom1, -1, self.acceptorParamValues[i])
        else:
            for angle in data.angles:
                type1 = data.atomType[data.atoms[angle[0]]]
                type2 = data.atomType[data.atoms[angle[1]]]
                type3 = data.atomType[data.atoms[angle[2]]]
                for i in range(len(self.acceptorTypes1)):
                    types1 = self.acceptorTypes1[i]
                    types2 = self.acceptorTypes2[i]
                    types3 = self.acceptorTypes3[i]
                    if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                        force.addAcceptor(angle[0], angle[1], angle[2], self.acceptorParamValues[i])

        # Add exclusions.

        for donor in range(force.getNumDonors()):
            (d1, d2, d3, params) = force.getDonorParameters(donor)
            outerAtoms = set((d1, d2, d3))
            if -1 in outerAtoms:
                outerAtoms.remove(-1)
            excludedAtoms = set(outerAtoms)
            for i in range(self.bondCutoff):
                newOuterAtoms = set()
                for atom in outerAtoms:
                    for bond in data.atomBonds[atom]:
                        b = data.bonds[bond]
                        bondedAtom = (b.atom2 if b.atom1 == atom else b.atom1)
                        if bondedAtom not in excludedAtoms:
                            newOuterAtoms.add(bondedAtom)
                            excludedAtoms.add(bondedAtom)
                outerAtoms = newOuterAtoms
            for acceptor in range(force.getNumAcceptors()):
                (a1, a2, a3, params) = force.getAcceptorParameters(acceptor)
                if a1 in excludedAtoms or a2 in excludedAtoms or a3 in excludedAtoms:
                    force.addExclusion(donor, acceptor)

parsers["CustomHbondForce"] = CustomHbondGenerator.parseElement

## @private
class CustomManyParticleGenerator(object):
    """A CustomManyParticleGenerator constructs a CustomManyParticleForce."""

    def __init__(self, forcefield, particlesPerSet, energy, permutationMode, bondCutoff):
        self.ff = forcefield
        self.particlesPerSet = particlesPerSet
        self.energy = energy
        self.permutationMode = permutationMode
        self.bondCutoff = bondCutoff
        self.globalParams = {}
        self.perParticleParams = []
        self.functions = []
        self.typeFilters = []

    @staticmethod
    def parseElement(element, ff):
        permutationMap = {"SinglePermutation" : mm.CustomManyParticleForce.SinglePermutation,
                          "UniqueCentralParticle" : mm.CustomManyParticleForce.UniqueCentralParticle}
        generator = CustomManyParticleGenerator(ff, int(element.attrib['particlesPerSet']), element.attrib['energy'], permutationMap[element.attrib['permutationMode']], int(element.attrib['bondCutoff']))
        ff.registerGenerator(generator)
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerParticleParameter'):
            generator.perParticleParams.append(param.attrib['name'])
        for param in element.findall('TypeFilter'):
            generator.typeFilters.append((int(param.attrib['index']), [int(x) for x in param.attrib['types'].split(',')]))
        generator.params = ForceField._AtomTypeParameters(ff, 'CustomManyParticleForce', 'Atom', generator.perParticleParams)
        generator.params.parseDefinitions(element)
        generator.functions += _parseFunctions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.CustomManyParticleForce.NoCutoff,
                     CutoffNonPeriodic:mm.CustomManyParticleForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.CustomManyParticleForce.CutoffPeriodic,
                     Ewald:mm.CustomManyParticleForce.CutoffPeriodic,
                     PME:mm.CustomManyParticleForce.CutoffPeriodic,
                     LJPME:mm.CustomManyParticleForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for CustomManyParticleForce')
        force = mm.CustomManyParticleForce(self.particlesPerSet, self.energy)
        force.setPermutationMode(self.permutationMode)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perParticleParams:
            force.addPerParticleParameter(param)
        for index, types in self.typeFilters:
            force.setTypeFilter(index, types)
        for (name, type, values, params) in self.functions:
            if type == 'Continuous1D':
                force.addTabulatedFunction(name, mm.Continuous1DFunction(values, params['min'], params['max'], params['periodic']))
            elif type == 'Continuous2D':
                force.addTabulatedFunction(name, mm.Continuous2DFunction(params['xsize'], params['ysize'], values, params['xmin'], params['xmax'], params['ymin'], params['ymax']))
            elif type == 'Continuous3D':
                force.addTabulatedFunction(name, mm.Continuous2DFunction(params['xsize'], params['ysize'], params['zsize'], values, params['xmin'], params['xmax'], params['ymin'], params['ymax'], params['zmin'], params['zmax']))
            elif type == 'Discrete1D':
                force.addTabulatedFunction(name, mm.Discrete1DFunction(values))
            elif type == 'Discrete2D':
                force.addTabulatedFunction(name, mm.Discrete2DFunction(params['xsize'], params['ysize'], values))
            elif type == 'Discrete3D':
                force.addTabulatedFunction(name, mm.Discrete2DFunction(params['xsize'], params['ysize'], params['zsize'], values))
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            type = int(self.params.getExtraParameters(atom, data)['filterType'])
            force.addParticle(values, type)
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # Create exclusions based on bonds.

        bondIndices = []
        for bond in data.bonds:
            bondIndices.append((bond.atom1, bond.atom2))

        # If a virtual site does *not* share exclusions with another atom, add a bond between it and its first parent atom.

        for i in range(sys.getNumParticles()):
            if sys.isVirtualSite(i):
                (site, atoms, excludeWith) = data.virtualSites[data.atoms[i]]
                if excludeWith is None:
                    bondIndices.append((i, site.getParticle(0)))

        # Certain particles, such as lone pairs and Drude particles, share exclusions with a parent atom.
        # If the parent atom does not interact with an atom, the child particle does not either.

        for atom1, atom2 in bondIndices:
            for child1 in data.excludeAtomWith[atom1]:
                bondIndices.append((child1, atom2))
                for child2 in data.excludeAtomWith[atom2]:
                    bondIndices.append((child1, child2))
            for child2 in data.excludeAtomWith[atom2]:
                bondIndices.append((atom1, child2))

        # Create the exclusions.

        nonbonded = [f for f in sys.getForces() if isinstance(f, mm.CustomManyParticleForce)][0]
        nonbonded.createExclusionsFromBonds(bondIndices, self.bondCutoff)

parsers["CustomManyParticleForce"] = CustomManyParticleGenerator.parseElement

def getAtomPrint(data, atomIndex):

    if (atomIndex < len(data.atoms)):
        atom = data.atoms[atomIndex]
        returnString = "%4s %4s %5d" % (atom.name, atom.residue.name, atom.residue.index)
    else:
        returnString = "NA"

    return returnString


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

    print("Constraints bond=%d angle=%d  total=%d" % (bondCount, angleCount, (bondCount+angleCount)))

## @private
class AmoebaBondGenerator(object):

    #=============================================================================================

    """An AmoebaBondGenerator constructs a AmoebaBondForce."""

    #=============================================================================================

    def __init__(self, cubic, quartic):
        self.builder = amoebaforces.AmoebaBondForceBuilder(cubic, quartic)

    @staticmethod
    def parseElement(element, forceField):
        # <AmoebaBondForce bond-cubic="-25.5" bond-quartic="379.3125">
        # <Bond class1="1" class2="2" length="0.1437" k="156900.0"/>
        generator = AmoebaBondGenerator(element.attrib['bond-cubic'], element.attrib['bond-quartic'])
        forceField._forces.append(generator)
        for bond in element.findall('Bond'):
            try:
                generator.builder.registerParams((bond.attrib['class1'], bond.attrib['class2']), float(bond.attrib['length']), float(bond.attrib['k']))
            except:
                outputString = "AmoebaBondGenerator: error getting types: %s %s" % (bond.attrib['class1'], bond.attrib['class2'])
                raise ValueError(outputString)
        generator.classNameForType = dict((t.name, int(t.atomClass)) for t in forceField._atomTypes.values())

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = self.builder.getForce(sys)
        bondsConstraints = []
        atomClasses = [str(self.classNameForType[data.atomType[atom]]) for atom in data.atoms]
        for bond in data.bonds:
            bondsConstraints.append(bond.isConstrained)
            bondType = (atomClasses[bond.atom1], atomClasses[bond.atom2])
            params = self.builder._findMatchingParams(self.builder.bondParams, bondType)
            bond.length = params[0]
            if bond.isConstrained:
                data.addConstraint(sys, bond.atom1, bond.atom2, params[0])
        self.builder.addBonds(force, atomClasses, _findBondsForExclusions(data, sys), bondsConstraints, args.get('flexibleConstraints', False))

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
                data.addConstraint(sys, angle[0], angle[2], length)
                return

## @private
class AmoebaAngleGenerator(object):
    """An AmoebaAngleGenerator constructs a AmoebaAngleForce."""

    def __init__(self, forceField, cubic, quartic, pentic, sextic):
        self.angleBuilder = amoebaforces.AmoebaAngleForceBuilder(float(cubic), float(quartic), float(pentic), float(sextic))
        self.inPlaneAngleBuilder = amoebaforces.AmoebaInPlaneAngleForceBuilder(float(cubic), float(quartic), float(pentic), float(sextic))
        self.forceField = forceField


    @staticmethod
    def parseElement(element, forceField):
        # <AmoebaAngleForce angle-cubic="-0.014" angle-quartic="5.6e-05" angle-pentic="-7e-07" angle-sextic="2.2e-08">
        #   <Angle class1="2" class2="1" class3="3" k="0.0637259642196" angle1="122.00"  />
        existing = [f for f in forceField._forces if isinstance(f, AmoebaAngleGenerator)]
        if len(existing) == 0:
            generator = AmoebaAngleGenerator(forceField, element.attrib['angle-cubic'], element.attrib['angle-quartic'],  element.attrib['angle-pentic'], element.attrib['angle-sextic'])
            forceField.registerGenerator(generator)
        else:
            generator = existing[0]
            if tuple(element.attrib[x] for x in ('angle-cubic', 'angle-quartic', 'angle-pentic', 'angle-sextic')) != (generator.cubic, generator.quartic, generator.pentic, generator.sextic):
                raise ValueError('All <AmoebaAngleForce> tags must use identical scale factors')
        generator.angleParams = {}
        for angle in element.findall('Angle'):
            try:
                theta0 = [float(angle.attrib['angle1'])]
                if 'angle2' in angle.attrib:
                    theta0.append(float(angle.attrib['angle2']))
                if 'angle3' in angle.attrib:
                    theta0.append(float(angle.attrib['angle3']))
                generator.angleParams[(angle.attrib['class1'], angle.attrib['class2'], angle.attrib['class3'])] = {"k": float(angle.attrib['k']), 
                                                                                                                   "theta0": theta0, 
                                                                                                                   "inPlane": angle.attrib.get('inPlane', None)}
            except Exception as e:
                outputString = "AmoebaAngleGenerator: error getting types: %s %s %s" % (
                                    angle.attrib['class1'],
                                    angle.attrib['class2'],
                                    angle.attrib['class3'])
                raise ValueError(outputString)
            
        generator.classNameForType = dict((t.name, int(t.atomClass)) for t in forceField._atomTypes.values())

    #=============================================================================================
    # createForce is bypassed here since the AmoebaOutOfPlaneBendForce generator must first execute
    # and partition angles into in-plane and non-in-plane angles
    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        if not any(isinstance(f, AmoebaOutOfPlaneBendGenerator) for f in self.forceField.getGenerators()):
            raise ValueError('A ForceField containing an <AmoebaAngleForce> must also contain an <AmoebaOutOfPlaneBendForce>')

    #=============================================================================================
    # createForcePostOpBendAngle is called by AmoebaOutOfPlaneBendForce with the list of
    # non-in-plane angles
    #=============================================================================================

    def createForcePostOpBendAngle(self, sys, data, nonbondedMethod, nonbondedCutoff, angleList, args):
        anglesConstraints = []
        angles = []
        atomClasses = [str(self.classNameForType[data.atomType[atom]]) for atom in data.atoms]
        genericAngleParams = {k: v for k, v in self.angleParams.items() if v.get('inPlane') != 'True'}

        for angleDict in angleList:
            angle = angleDict['angle']
            if angleDict['inPlane']:
                print("in-plane angle found in non-in-plane list", angle)
                continue

            isConstrained = angleDict['isConstrained']
            angleClasses = (atomClasses[angle[0]], atomClasses[angle[1]], atomClasses[angle[2]])

            params = genericAngleParams.get(angleClasses) or genericAngleParams.get(angleClasses[::-1])
            if params is None:
                outputString = "AmoebaAngleGenerator: no parameters found for angle: %s-%s-%s (atom %s, %s, %s)" % (
                    angleClasses[0],
                    angleClasses[2],
                    getAtomPrint(data, angle[0]),
                    getAtomPrint(data, angle[1]),
                    getAtomPrint(data, angle[2])
                )
                raise ValueError(outputString)

            angleTuple = (angle[0], angle[1], angle[2])
            idealAngle = None

            if isConstrained and params["k"] != 0.0:
                idealAngle = params["theta0"][0] * math.pi / 180.0
                addAngleConstraint(angle, idealAngle, data, sys)

            if params["k"] != 0.0 and (not isConstrained or args.get('flexibleConstraints', False)):
                idealAngle = self.angleBuilder.getIdealAngle(angleTuple, params["theta0"], data) * math.pi / 180
                self.angleBuilder.registerParams(angleTuple, (idealAngle, params["k"]))

            anglesConstraints.append(isConstrained)
            angles.append(angle)
            angleDict['idealAngle'] = idealAngle

        force = self.angleBuilder.getForce(sys)
        self.angleBuilder.addAngles(force, angles, anglesConstraints, args.get('flexibleConstraints', False))



        """
        type1 = data.atomType[data.atoms[angle[0]]]
        type2 = data.atomType[data.atoms[angle[1]]]
        type3 = data.atomType[data.atoms[angle[2]]]
        for i in range(len(self.types1)):
            # self.inPlane is used for modern force fields.  inPlane is used for legacy ones that don't specify it.
            if self.inPlane[i] or (self.inPlane[i] is None and inPlane):
                continue
            types1 = self.types1[i]
            types2 = self.types2[i]
            types3 = self.types3[i]
            if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                if isConstrained and self.k[i] != 0.0:
                    angleDict['idealAngle'] = self.angle[i][0]
                    addAngleConstraint(angle, self.angle[i][0]*DEG_TO_RAD, data, sys)


                if self.k[i] != 0 and (not isConstrained or args.get('flexibleConstraints', False)):
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
            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]

                            angleValue =  self.angle[i][numberOfHydrogens]
                        else:
                            outputString = "AmoebaAngleGenerator angle index=%d is out of range: [0, %5d] " % (numberOfHydrogens, lenAngle)
                            raise ValueError(outputString)
                    else:
                        angleValue =  self.angle[i][0]

                    angleDict['idealAngle'] = angleValue
                    force.addAngle(angle[0], angle[1], angle[2], [angleValue, self.k[i]])
                break
        """

    #=============================================================================================
    # createForcePostOpBendInPlaneAngle is called by AmoebaOutOfPlaneBendForce with the list of
    # in-plane angles
    #=============================================================================================

    def createForcePostOpBendInPlaneAngle(self, sys, data, nonbondedMethod, nonbondedCutoff, angleList, args):
        inPlaneAngleParams = {k: v for k, v in self.angleParams.items() if v.get('inPlane') == 'True'}
        atomClasses = [str(self.classNameForType[data.atomType[atom]]) for atom in data.atoms]
        inPlaneAnglesConstraints = []
        inPlaneAngles = []
        for angleDict in angleList:
            angle = angleDict['angle']
            isConstrained = angleDict['isConstrained']
            inPlane = angleDict['inPlane']

            if not inPlane:
                continue

            angleClasses = (atomClasses[angle[0]], atomClasses[angle[1]], atomClasses[angle[2]], atomClasses[angle[3]])
            params = inPlaneAngleParams.get(angleClasses[:3], None) or inPlaneAngleParams.get(angleClasses[:3][::-1], None)

            if params is None:
                outputString = "AmoebaAngleGenerator: no parameters found for in-plane angle: %s-%s-%s (atom %s, %s, %s)" % (
                                    angleClasses[0],
                                    angleClasses[2],
                                    getAtomPrint(data, angle[0]),
                                    getAtomPrint(data, angle[1]),
                                    getAtomPrint(data, angle[2]))
                raise ValueError(outputString)

            if isConstrained and params["k"] != 0.0:
                addAngleConstraint(angle, idealAngle, data, sys)

            idealAngle = params["theta0"][0]*math.pi/180.0
            inPlaneAnglesConstraints.append(isConstrained)
            inPlaneAngles.append(angle)
            angleDict['idealAngle'] = idealAngle
            self.inPlaneAngleBuilder.registerParams(angleClasses, (idealAngle, params["k"]))

        force = self.inPlaneAngleBuilder.getForce(sys)
        self.inPlaneAngleBuilder.addInPlaneAngles(force, atomClasses, inPlaneAngles, inPlaneAnglesConstraints, args.get('flexibleConstraints', False))

        """
        type1 = data.atomType[data.atoms[angle[0]]]
        type2 = data.atomType[data.atoms[angle[1]]]
        type3 = data.atomType[data.atoms[angle[2]]]

        for i in range(len(self.types1)):
            # self.inPlane is used for modern force fields.  inPlane is used for legacy ones that don't specify it.
            if self.inPlane[i] == False or (self.inPlane[i] is None and not inPlane):
                continue
            types1 = self.types1[i]
            types2 = self.types2[i]
            types3 = self.types3[i]

            if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                angleDict['idealAngle'] = self.angle[i][0]
                if (isConstrained and self.k[i] != 0.0):
                    addAngleConstraint(angle, self.angle[i][0]*math.pi/180.0, data, sys)
                if self.k[i] != 0.0 and (not isConstrained or args.get('flexibleConstraints', False)):
                    force.addBond((angle[0], angle[1], angle[2], angle[3]), (self.angle[i][0], self.k[i]))
                break
        """

parsers["AmoebaAngleForce"] = AmoebaAngleGenerator.parseElement

#=============================================================================================
# Generator for the AmoebaOutOfPlaneBend covalent force; also calls methods in the
# AmoebaAngleGenerator to generate the AmoebaAngleForce and
# AmoebaInPlaneAngleForce
#=============================================================================================
## @private
class AmoebaOutOfPlaneBendGenerator(object):
    """An AmoebaOutOfPlaneBendGenerator constructs a AmoebaOutOfPlaneBendForce."""

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

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaOutOfPlaneBendForce type="ALLINGER" opbend-cubic="-0.014" opbend-quartic="5.6e-05" opbend-pentic="-7e-07" opbend-sextic="2.2e-08">
        #   <Angle class1="2" class2="1" class3="0" class4="0" k="0.0531474541591"/>
        #   <Angle class1="3" class2="1" class3="0" class4="0" k="0.0898536095496"/>

        # get global scalar parameters

        existing = [f for f in forceField._forces if isinstance(f, AmoebaOutOfPlaneBendGenerator)]
        if len(existing) == 0:
            generator = AmoebaOutOfPlaneBendGenerator(forceField, element.attrib['type'], element.attrib['opbend-cubic'], element.attrib['opbend-quartic'],  element.attrib['opbend-pentic'], element.attrib['opbend-sextic'])
            forceField.registerGenerator(generator)
        else:
            generator = existing[0]
            if tuple(element.attrib[x] for x in ('type', 'opbend-cubic', 'opbend-quartic', 'opbend-pentic', 'opbend-sextic')) != (generator.type, generator.cubic, generator.quartic, generator.pentic, generator.sextic):
                raise ValueError('All <AmoebaOutOfPlaneBendForce> tags must use identical scale factors')

        for angle in element.findall('Angle'):
            if 'class3' in angle.attrib and 'class4' in angle.attrib and angle.attrib['class3'] == '0' and angle.attrib['class4'] == '0':
                # This is needed for backward compatibility with old AMOEBA force fields that specified wildcards in a nonstandard way.
                angle.attrib['class3'] = ''
                angle.attrib['class4'] = ''
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

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        self._nonbondedMethod = nonbondedMethod
        self._nonbondedCutoff = nonbondedCutoff

    def postprocessSystem(self, sys, data, args):
        # We need to wait until after all bonds have been added so their lengths will be set correctly.

        builder = amoebaforces.AmoebaOutOfPlaneBendForceBuilder(self.cubic, self.quartic, self.pentic, self.sextic)
        force = builder.getForce(sys)

        # this hash is used to ensure the out-of-plane-bend bonds
        # are only added once

        skipAtoms = dict()
        angles = []

        def addBond(particles):
            types = [data.atomType[data.atoms[p]] for p in particles]
            for i in range(len(self.types1)):
                if types[1] in self.types2[i] and types[3] in self.types1[i]:
                    if (types[0] in self.types3[i] and types[2] in self.types4[i]) or (types[2] in self.types3[i] and types[0] in self.types4[i]):
                        force.addBond(particles, [self.ks[i]])
                        return

        for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):

            middleAtom = angle[1]
            middleType = data.atomType[data.atoms[middleAtom]]
            middleCovalency = len(data.atomBonds[middleAtom])

            # if middle atom has covalency of 3 and
            # the types of the middle atom and the partner atom (atom bonded to
            # middle atom, but not in angle) match types1 and types2, then
            # three out-of-plane bend angles are generated. Three in-plane angle
            # are also generated. If the conditions are not satisfied, the angle is marked as 'generic' angle (not a in-plane angle)

            if middleCovalency == 3 and middleAtom not in skipAtoms:

                partners = []

                for bond in data.atomBonds[middleAtom]:
                    atom1 = data.bonds[bond].atom1
                    atom2 = data.bonds[bond].atom2
                    if atom1 != middleAtom:
                        partner = atom1
                    else:
                        partner = atom2

                    partnerType = data.atomType[data.atoms[partner]]
                    for i in range(len(self.types1)):
                        types1 = self.types1[i]
                        types2 = self.types2[i]
                        if middleType in types2 and partnerType in types1:
                            partners.append(partner)
                            break

                if len(partners) == 3:

                    addBond([partners[0], middleAtom, partners[1], partners[2]])
                    addBond([partners[2], middleAtom, partners[0], partners[1]])
                    addBond([partners[1], middleAtom, partners[2], partners[0]])

                    # skipAtoms is used to ensure angles are only included once

                    skipAtoms[middleAtom] = set(partners[:3])

                    # in-plane angle

                    angleDict = {}
                    angleList = list(angle[:3])
                    for atomIndex in partners:
                        if atomIndex not in angleList:
                            angleList.append(atomIndex)
                    angleDict['angle'] = angleList
                    angleDict['isConstrained'] = 0
                    angleDict['inPlane'] = True
                    angles.append(angleDict)

                else:
                    angleDict = {}
                    angleList = list(angle[:3])
                    for atomIndex in partners:
                        if atomIndex not in angleList:
                            angleList.append(atomIndex)
                    angleDict['angle'] = angleList
                    angleDict['isConstrained'] = isConstrained
                    angleDict['inPlane'] = False
                    angles.append(angleDict)
            elif middleCovalency == 3 and middleAtom in skipAtoms:

                angleDict = {}
                angleList = list(angle[:3])
                for atomIndex in skipAtoms[middleAtom]:
                    if atomIndex not in angleList:
                        angleList.append(atomIndex)
                angleDict['angle'] = angleList
                angleDict['isConstrained'] = isConstrained
                angleDict['inPlane'] = True
                angles.append(angleDict)

            else:
                angleDict = {}
                angleDict['angle'] = angle
                angleDict['isConstrained'] = isConstrained
                angleDict['inPlane'] = False
                angles.append(angleDict)

        # get AmoebaAngleGenerator and add AmoebaAngle and AmoebaInPlaneAngle forces

        for force in self.forceField._forces:
            if (force.__class__.__name__ == 'AmoebaAngleGenerator'):
                force.createForcePostOpBendAngle(sys, data, self._nonbondedMethod, self._nonbondedCutoff, angles, args)
                force.createForcePostOpBendInPlaneAngle(sys, data, self._nonbondedMethod, self._nonbondedCutoff, angles, args)

        for force in self.forceField._forces:
            if (force.__class__.__name__ == 'AmoebaStretchBendGenerator'):
                force.createForcePostAmoebaBondForce(sys, data, self._nonbondedMethod, self._nonbondedCutoff, angles, args)

parsers["AmoebaOutOfPlaneBendForce"] = AmoebaOutOfPlaneBendGenerator.parseElement


## @private
class AmoebaTorsionGenerator(object):
    """An AmoebaTorsionGenerator constructs a AmoebaTorsionForce."""

    def __init__(self, torsionUnit):
        self.torsionUnit = torsionUnit

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.types4 = []

        self.t1 = []
        self.t2 = []
        self.t3 = []

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
            types = forceField._findAtomTypes(torsion.attrib, 4)
            if None not in types:

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
                                    torsion.attrib['class1'],
                                    torsion.attrib['class2'],
                                    torsion.attrib['class3'],
                                    torsion.attrib['class4'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nontorsionedMethod, nontorsionedCutoff, args):
        builder = amoebaforces.AmoebaTorsionForceBuilder()
        force = builder.getForce(sys)

        for torsion in data.propers:

            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]

            for i in range(len(self.types1)):

                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                types4 = self.types4[i]

                # match types in forward or reverse direction

                if (type1 in types1 and type2 in types2 and type3 in types3 and type4 in types4) or (type4 in types1 and type3 in types2 and type2 in types3 and type1 in types4):
                    if self.t1[i][0] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], 1, self.t1[i][1], self.t1[i][0])
                    if self.t2[i][0] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], 2, self.t2[i][1], self.t2[i][0])
                    if self.t3[i][0] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], 3, self.t3[i][1], self.t3[i][0])
                    break

parsers["AmoebaTorsionForce"] = AmoebaTorsionGenerator.parseElement

## @private
class AmoebaPiTorsionGenerator(object):
    """An AmoebaPiTorsionGenerator constructs a AmoebaPiTorsionForce."""

    def __init__(self):
        self.builder = amoebaforces.AmoebaPiTorsionForceBuilder()

    @staticmethod
    def parseElement(element, forceField):
        generator = AmoebaPiTorsionGenerator()
        forceField._forces.append(generator)

        for piTorsion in element.findall('PiTorsion'):
            # TODO: make it read pitorsionunit
            try:
                generator.builder.registerParams((piTorsion.attrib['class1'], piTorsion.attrib['class2']), (float(piTorsion.attrib['k']),))
            except:
                outputString = "AmoebaPiTorsionGenerator: error getting types: %s %s " % (
                                    piTorsion.attrib['class1'],
                                    piTorsion.attrib['class2'])
                raise ValueError(outputString)
        generator.classNameForType = dict((t.name, int(t.atomClass)) for t in forceField._atomTypes.values())

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        atomClasses = [str(self.classNameForType[data.atomType[atom]]) for atom in data.atoms]
        bondedToAtom = [bondedToAtom for bondedToAtom in data.bondedToAtom]
        processedPiTorsions = self.builder.getAllPiTorsions(atomClasses, bondedToAtom, _findBondsForExclusions(data, sys))
        force = self.builder.getForce(sys)
        self.builder.addPiTorsions(force, atomClasses, processedPiTorsions)

parsers["AmoebaPiTorsionForce"] = AmoebaPiTorsionGenerator.parseElement

## @private
class AmoebaStretchTorsionGenerator(object):
    """An AmoebaStretchTorsionGenerator constructs a AmoebaStretchTorsionForce."""

    def __init__(self):
        self.builder = amoebaforces.AmoebaStretchTorsionForceBuilder()

    @staticmethod
    def parseElement(element, forceField):
        # <Torsion class1="44" class2="46" class3="68" class4="65" v11="0.0" v12="0.0" v13="62.760000000000005" v21="0.0" v22="0.0" v23="-167.36" v31="0.0" v32="0.0" v33="217.568"/>
        generator = AmoebaStretchTorsionGenerator()
        forceField._forces.append(generator)
        for torsion in element.findall('Torsion'):
            try:
                params = tuple(float(torsion.attrib[p]) for p in ('v11', 'v12', 'v13', 'v21', 'v22', 'v23', 'v31', 'v32', 'v33'))
                generator.builder.registerParams((torsion.attrib['class1'], torsion.attrib['class2'], torsion.attrib['class3'], torsion.attrib['class4']), params)
            except Exception as e:
                outputString = "AmoebaStretchTorsionGenerator: error getting types: %s %s %s %s" % (
                                    torsion.attrib['class1'],
                                    torsion.attrib['class2'],
                                    torsion.attrib['class3'],
                                    torsion.attrib['class4'])
                raise ValueError(outputString)
     
        generator.classNameForType = dict((t.name, int(t.atomClass)) for t in forceField._atomTypes.values())

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        pass

    def postprocessSystem(self, sys, data, args):
        # We need to wait until after all bonds and torsions have been added before adding the stretch-torsions,
        # since it needs parameters from them.
        atomClasses = [str(self.classNameForType[data.atomType[atom]]) for atom in data.atoms]
        force = self.builder.getForce(sys)
        self.builder.addStretchTorsions(sys, force, atomClasses, data.propers)

parsers["AmoebaStretchTorsionForce"] = AmoebaStretchTorsionGenerator.parseElement

## @private
class AmoebaAngleTorsionGenerator(object):
    """An AmoebaAngleTorsionGenerator constructs a AmoebaAngleTorsionForce."""

    def __init__(self):
        self.builder = amoebaforces.AmoebaAngleTorsionForceBuilder()

    @staticmethod
    def parseElement(element, forceField):
        # <AmoebaAngleTorsionForce>
        #  <Torsion class1="44" class2="46" class3="68" class4="65" v11="3.3555680000000003" v12="0.0" v13="-13.903432" v21="0.0" v22="0.0" v23="-2.63592"/>
        generator = AmoebaAngleTorsionGenerator()
        forceField._forces.append(generator)
        for torsion in element.findall('Torsion'):
            try:
                params = tuple(float(torsion.attrib[p]) for p in ('v11', 'v12', 'v13', 'v21', 'v22', 'v23'))
                generator.builder.registerParams((torsion.attrib['class1'], torsion.attrib['class2'], torsion.attrib['class3'], torsion.attrib['class4']), params)                            
            except:
                outputString = "AmoebaAngleTorsionGenerator: error getting types: %s %s %s %s" % (
                                    torsion.attrib['class1'],
                                    torsion.attrib['class2'],
                                    torsion.attrib['class3'],
                                    torsion.attrib['class4'])
                raise ValueError(outputString)
        generator.classNameForType = dict((t.name, int(t.atomClass)) for t in forceField._atomTypes.values())

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        pass

    def postprocessSystem(self, sys, data, args):
        # We need to wait until after all angles and torsions have been added before adding the angle-torsions,
        # since it needs parameters from them.
        atomClasses = [str(self.classNameForType[data.atomType[atom]]) for atom in data.atoms]
        force = self.builder.getForce(sys)
        self.builder.addAngleTorsions(sys, force, atomClasses, data.propers)

parsers["AmoebaAngleTorsionForce"] = AmoebaAngleTorsionGenerator.parseElement

## @private
class AmoebaTorsionTorsionGenerator(object):
    """An AmoebaTorsionTorsionGenerator constructs a AmoebaTorsionTorsionForce."""

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.types4 = []
        self.types5 = []

        self.gridIndex = []

        self.grids = []

    @staticmethod
    def parseElement(element, forceField):

        generator = AmoebaTorsionTorsionGenerator()
        forceField._forces.append(generator)
        maxGridIndex = -1

        # <AmoebaTorsionTorsionForce >
        # <TorsionTorsion class1="3" class2="1" class3="2" class4="3" class5="1" grid="0" nx="25" ny="25" />

        for torsionTorsion in element.findall('TorsionTorsion'):
            types = forceField._findAtomTypes(torsionTorsion.attrib, 5)
            if None not in types:

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
                if 'fx' in gridEntry.attrib:
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

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        builder = amoebaforces.AmoebaTorsionTorsionForceBuilder()
        force = builder.getForce(sys)
        
        # Add torsion-torsion interactions
        builder.addTorsionTorsionInteractions(force, data, self.types1, self.types2, self.types3, 
                                              self.types4, self.types5, self.gridIndex, sys)
        
        # Set grids
        for (index, grid) in enumerate(self.grids):
            builder.setTorsionTorsionGrid(force, index, grid)

parsers["AmoebaTorsionTorsionForce"] = AmoebaTorsionTorsionGenerator.parseElement

## @private
class AmoebaStretchBendGenerator(object):
    """An AmoebaStretchBendGenerator constructs a AmoebaStretchBendForce."""

    def __init__(self, forcefield):
        self.forcefield = forcefield
        self.builder = amoebaforces.AmoebaStretchBendForceBuilder()

    @staticmethod
    def parseElement(element, forceField):
        generator = AmoebaStretchBendGenerator(forceField)
        forceField._forces.append(generator)

        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />
        generator.stretchBendParams = {}
        for stretchBend in element.findall('StretchBend'):
            try:
                class1 = stretchBend.attrib['class1']
                class2 = stretchBend.attrib['class2']
                class3 = stretchBend.attrib['class3']
                k1 = float(stretchBend.attrib['k1'])
                k2 = float(stretchBend.attrib['k2'])
                generator.stretchBendParams[(class1, class2, class3)] = {'k1': k1, 'k2': k2}
            except:
                outputString = "AmoebaStretchBendGenerator : error getting types: %s %s %s" % (
                                    stretchBend.attrib['class1'],
                                    stretchBend.attrib['class2'],
                                    stretchBend.attrib['class3'])
                raise ValueError(outputString)
        generator.classNameForType = dict((t.name, int(t.atomClass)) for t in forceField._atomTypes.values())

    #=============================================================================================

    # The setup of this force is dependent on AmoebaBondForce and AmoebaAngleForce
    # having been called since the ideal bond lengths and angle are needed here.
    # As a conseqeunce, createForce() is not implemented since it is not guaranteed that the generator for
    # AmoebaBondForce and AmoebaAngleForce have been called prior to AmoebaStretchBendGenerator().
    # Instead, createForcePostAmoebaBondForce() is called
    # after the generators for AmoebaBondForce and AmoebaAngleForce have been called

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        if not any(isinstance(f, AmoebaOutOfPlaneBendGenerator) for f in self.forcefield.getGenerators()):
            raise ValueError('A ForceField containing an <AmoebaStretchBendForce> must also contain an <AmoebaOutOfPlaneBendForce>')

    #=============================================================================================

    # Note: request for constrained bonds is ignored.

    #=============================================================================================

    def createForcePostAmoebaBondForce(self, sys, data, nonbondedMethod, nonbondedCutoff, angleList, args):
        atomClasses = [str(self.classNameForType[data.atomType[atom]]) for atom in data.atoms]
        angles = []
        idealAngles = {}
        for angleDict in angleList:
            angle = tuple(angleDict['angle'][:3])
            angles.append(angle)
            idealAngle = angleDict.get('idealAngle')
            if idealAngle is None:
                outputString = "AmoebaStretchBendGenerator: ideal angle is not set for following entry:\n"
                outputString += "   types: %5s %5s %5s atoms: " % (data.atomType[data.atoms[angle[0]]],
                                                                data.atomType[data.atoms[angle[1]]],
                                                                data.atomType[data.atoms[angle[2]]])
                outputString += getAtomPrint( data, angle[0] ) + ' '
                outputString += getAtomPrint( data, angle[1] ) + ' '
                outputString += getAtomPrint( data, angle[2] )
                raise ValueError(outputString)
            idealAngles[angle] = idealAngle*math.pi/180.0

        bondParams = {
            (atomClasses[bond.atom1], atomClasses[bond.atom2]): {"r0": bond.length}
            for bond in data.bonds
        }

        processedAngles = self.builder.registerAllStretchBendParams(atomClasses, angles, self.stretchBendParams, bondParams, idealAngles)
        force = self.builder.getForce(sys)
        self.builder.addStretchBends(force, processedAngles)

parsers["AmoebaStretchBendForce"] = AmoebaStretchBendGenerator.parseElement


## @private
class AmoebaVdwGenerator(object):
    """A AmoebaVdwGenerator constructs a AmoebaVdwForce."""

    def __init__(self, type, radiusrule, radiustype, radiussize, epsilonrule, vdw13Scale, vdw14Scale, vdw15Scale):
        self.builder = amoebaforces.AmoebaVdwForceBuilder(str(type), 
                                                          str(radiusrule), 
                                                          str(radiustype), 
                                                          str(radiussize), 
                                                          str(epsilonrule), 
                                                          float(vdw13Scale), 
                                                          float(vdw14Scale), 
                                                          float(vdw15Scale))

    @staticmethod
    def parseElement(element, forceField):
        # <AmoebaVdwForce type="BUFFERED-14-7" radiusrule="CUBIC-MEAN" radiustype="R-MIN" radiussize="DIAMETER" epsilonrule="HHG" vdw-13-scale="0.0" vdw-14-scale="1.0" vdw-15-scale="1.0" >
        #  <Vdw class="1" sigma="0.371" epsilon="0.46024000000000004" reduction="1.0"/>
        #  <Vdw class="2" sigma="0.382" epsilon="0.42258400000000007" reduction="1.0"/>

        existing = [f for f in forceField._forces if isinstance(f, AmoebaVdwGenerator)]
        if len(existing) == 0:
            generator = AmoebaVdwGenerator(element.attrib['type'], element.attrib['radiusrule'], element.attrib['radiustype'], element.attrib['radiussize'], element.attrib['epsilonrule'],
                                           element.attrib['vdw-13-scale'], element.attrib['vdw-14-scale'], element.attrib['vdw-15-scale'])
            forceField.registerGenerator(generator)
        else:
            # Multiple <AmoebaVdwForce> tags were found, probably in different files. Simply add more types to the existing one.
            generator = existing[0]
            if abs(generator.vdw13Scale - float(element.attrib['vdw-13-scale'])) > NonbondedGenerator.SCALETOL or \
                    abs(generator.vdw14Scale - float(element.attrib['vdw-14-scale'])) > NonbondedGenerator.SCALETOL or \
                    abs(generator.vdw15Scale - float(element.attrib['vdw-15-scale'])) > NonbondedGenerator.SCALETOL:
                raise ValueError('Found multiple AmoebaVdwForce tags with different scale factors')
            if generator.radiusrule != element.attrib['radiusrule'] or generator.epsilonrule != element.attrib['epsilonrule'] or \
                    generator.radiustype != element.attrib['radiustype'] or generator.radiussize != element.attrib['radiussize']:
                raise ValueError('Found multiple AmoebaVdwForce tags with different combining rules')
        for vdw in element.findall('Vdw'):
            generator.builder.registerClassParams(vdw.attrib['class'], float(vdw.attrib['sigma']), float(vdw.attrib['epsilon']), float(vdw.attrib['reduction']))
        for pair in element.findall('Pair'):
            generator.builder.registerPairParams(pair.attrib['class1'], pair.attrib['class2'], float(pair.attrib['sigma']), float(pair.attrib['epsilon']))
        generator.classNameForType = dict((t.name, int(t.atomClass)) for t in forceField._atomTypes.values())

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = self.builder.getForce(sys, nonbondedMethod, args.get('vdwCutoff', nonbondedCutoff), args.get('useDispersionCorrection', True))
        atomClasses = [str(self.classNameForType[data.atomType[atom]]) for atom in data.atoms]
        self.builder.addParticles(force, atomClasses, data.atoms, _findBondsForExclusions(data, sys))

parsers["AmoebaVdwForce"] = AmoebaVdwGenerator.parseElement

## @private
class AmoebaMultipoleGenerator(object):
    """A AmoebaMultipoleGenerator constructs an AmoebaMultipoleForce."""

    def __init__(self, forceField):
        self.builder = amoebaforces.AmoebaMultipoleForceBuilder()
        self.multipoleType = defaultdict(list)
        self.polarizationType = {}

    @staticmethod
    def parseElement(element, forceField):
        existing = [f for f in forceField._forces if isinstance(f, AmoebaMultipoleGenerator)]
        if len(existing) == 0:
            generator = AmoebaMultipoleGenerator(forceField)
            forceField.registerGenerator(generator)
        else:
            # Multiple <AmoebaMultipoleForce> tags were found, probably in different files. Simply add more types to the existing one.
            generator = existing[0]

        # set type map: [ kIndices, multipoles, AMOEBA/OpenMM axis type]

        for atom in element.findall('Multipole'):
            types = forceField._findAtomTypes(atom.attrib, 1)
            if None not in types:
                # k-indices not provided default to 0

                kIndices = [int(atom.attrib['type'])]
                kStrings = [ 'kz', 'kx', 'ky' ]
                for kString in kStrings:
                    kIndices.append(int(atom.attrib.get(kString,0)))

                # set multipole

                charge = float(atom.attrib['c0'])
                conversion = 1.0
                dipole = [ conversion*float(atom.attrib['d1']), conversion*float(atom.attrib['d2']), conversion*float(atom.attrib['d3'])]
                quadrupole = [conversion*float(atom.attrib[a]) for a in ['q11', 'q21', 'q31', 'q21', 'q22', 'q32', 'q31', 'q32', 'q33']]
                for t in types[0]:
                    generator.multipoleType[t].append({'kIndices':kIndices, 'charge':charge, 'dipole':dipole, 'quadrupole':quadrupole})
            else:
                outputString = "AmoebaMultipoleGenerator: error getting type for multipole: %s" % (atom.attrib['class'])
                raise ValueError(outputString)

        # polarization parameters

        for atom in element.findall('Polarize'):
            types = forceField._findAtomTypes(atom.attrib, 1)
            if None not in types:
                polarizability = float(atom.attrib['polarizability'])
                thole = float(atom.attrib['thole'])
                groupTypes = []
                for index in range(1, 7):
                    pgrp = f'pgrp{index}'
                    if pgrp in atom.attrib:
                        groupTypes.append(int(atom.attrib[pgrp]))
                for t in types[0]:
                    if t in generator.polarizationType:
                        raise ValueError(f'Multipole Polarize tags found for type {t}')
                    generator.polarizationType[t] = {'polarizability': polarizability, 'thole':thole, 'groupTypes':groupTypes}
            else:
                outputString = "AmoebaMultipoleGenerator: error getting type for polarize: %s" % (atom.attrib['class'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        builder = amoebaforces.AmoebaMultipoleForceBuilder()
        for type, params in self.multipoleType.items():
            for p in params:
                builder.registerMultipoleParams(int(type), amoebaforces.MultipoleParams(p['kIndices'], p['charge'], p['dipole'], p['quadrupole']))
        for type, params in self.polarizationType.items():
            builder.registerPolarizationParams(int(type), amoebaforces.PolarizationParams(params['polarizability'], params['thole'], params['groupTypes']))
        force = builder.getForce(sys, nonbondedMethod, nonbondedCutoff, args.get('ewaldErrorTolerance', 1e-4), args.get('polarization', 'mutual'),
                                 args.get('mutualInducedTargetEpsilon', 1e-5), args.get('mutualInducedMaxIterations', 60))
        atomTypes = [int(data.atomType[atom]) for atom in data.atoms]
        builder.addMultipoles(force, atomTypes, data.atoms, _findBondsForExclusions(data, sys))

parsers["AmoebaMultipoleForce"] = AmoebaMultipoleGenerator.parseElement

## @private
class AmoebaWcaDispersionGenerator(object):
    """A AmoebaWcaDispersionGenerator constructs a AmoebaWcaDispersionForce."""

    def __init__(self, epso, epsh, rmino, rminh, awater, slevy, dispoff, shctd):
        self.builder = amoebaforces.AmoebaWcaDispersionForceBuilder(float(epso), 
                                                                    float(epsh), 
                                                                    float(rmino), 
                                                                    float(rminh), 
                                                                    float(awater), 
                                                                    float(slevy), 
                                                                    float(dispoff), 
                                                                    float(shctd))
    @staticmethod
    def parseElement(element, forceField):
        #  <AmoebaWcaDispersionForce epso="0.46024" epsh="0.056484" rmino="0.17025" rminh="0.13275" awater="33.428" slevy="1.0"  dispoff="0.026" shctd="0.81" >
        #   <WcaDispersion class="1" radius="0.1855" epsilon="0.46024" />
        #   <WcaDispersion class="2" radius="0.191" epsilon="0.422584" />
        existing = [f for f in forceField._forces if isinstance(f, AmoebaWcaDispersionGenerator)]
        if len(existing) == 0:
            generator = AmoebaWcaDispersionGenerator(element.attrib['epso'],
                                                     element.attrib['epsh'],
                                                     element.attrib['rmino'],
                                                     element.attrib['rminh'],
                                                     element.attrib['awater'],
                                                     element.attrib['slevy'],
                                                     element.attrib['dispoff'],
                                                     element.attrib['shctd'])
            forceField.registerGenerator(generator)
        else:
            # Multiple <AmoebaWcaDispersionForce> tags were found, probably in different files. Simply add more types to the existing one.
            generator = existing[0]

        for wca in element.findall('WcaDispersion'):
            generator.builder.registerParams(int(wca.attrib['class']), float(wca.attrib['radius']), float(wca.attrib['epsilon']))
    
        generator.classNameForType = dict((t.name, int(t.atomClass)) for t in forceField._atomTypes.values())

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        atomClasses = [self.classNameForType[data.atomType[atom]] for atom in data.atoms]
        force = self.builder.getForce(sys)
        self.builder.addParticles(force, atomClasses)

parsers["AmoebaWcaDispersionForce"] = AmoebaWcaDispersionGenerator.parseElement

## @private
class AmoebaGeneralizedKirkwoodGenerator(object):
    """A AmoebaGeneralizedKirkwoodGenerator constructs a AmoebaGeneralizedKirkwoodForce."""

    def __init__(self, forceField, solventDielectric, soluteDielectric, includeCavityTerm, probeRadius, surfaceAreaFactor):
        self.builder = amoebaforces.AmoebaGeneralizedKirkwoodForceBuilder(float(solventDielectric), 
                                                                          float(soluteDielectric),
                                                                          bool(includeCavityTerm),
                                                                          float(probeRadius),
                                                                          float(surfaceAreaFactor))
        self.forceField = forceField
        
       
    @staticmethod
    def parseElement(element, forceField):
        #  <AmoebaGeneralizedKirkwoodForce solventDielectric="78.3" soluteDielectric="1.0" includeCavityTerm="1" probeRadius="0.14" surfaceAreaFactor="-170.351730663">
        #   <GeneralizedKirkwood type="1" charge="-0.22620" shct="0.79"  />
        #   <GeneralizedKirkwood type="2" charge="-0.15245" shct="0.72"  />
        existing = [f for f in forceField._forces if isinstance(f, AmoebaGeneralizedKirkwoodGenerator)]
        if len(existing) == 0:
            generator = AmoebaGeneralizedKirkwoodGenerator(forceField, element.attrib['solventDielectric'],
                                                           element.attrib['soluteDielectric'],
                                                           element.attrib['includeCavityTerm'],
                                                           element.attrib['probeRadius'],
                                                           element.attrib['surfaceAreaFactor'])
            forceField.registerGenerator(generator)
        else:
            # Multiple <AmoebaGeneralizedKirkwoodFprce> tags were found, probably in different files.  Simply add more types to the existing one.
            generator = existing[0]

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        if( nonbondedMethod != NoCutoff ):
            raise ValueError( "Only the nonbondedMethod=NoCutoff option is available for implicit solvent simulations." )

        # check if AmoebaMultipoleForce exists since charges needed
        # if it has not been created, raise an error

        amoebaMultipoleForceList = [f for f in sys.getForces() if type(f) == mm.AmoebaMultipoleForce]
        if (len(amoebaMultipoleForceList) > 0):
            amoebaMultipoleForce = amoebaMultipoleForceList[0]
        else:
            # call AmoebaMultipoleForceGenerator.createForce() to ensure charges have been set
            for force in self.forceField._forces:
                if (force.__class__.__name__ == 'AmoebaMultipoleGenerator'):
                    force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
                    break
      
        # Get the charges from the AmoebaMultipoleForce and register them with the builder
        for atomIndex in range(amoebaMultipoleForce.getNumMultipoles()):
            multipoleParameters = amoebaMultipoleForce.getMultipoleParameters(atomIndex)
            self.builder.registerParams(multipoleParameters[0])

        force = self.builder.getForce(sys, implicitSolvent=True)
        if force is not None:
            self.builder.addParticles(force, data.atoms, _findBondsForExclusions(data, sys))

parsers["AmoebaGeneralizedKirkwoodForce"] = AmoebaGeneralizedKirkwoodGenerator.parseElement

## @private
class AmoebaUreyBradleyGenerator(object):
    """An AmoebaUreyBradleyGenerator constructs a AmoebaUreyBradleyForce."""

    def __init__(self):
        self.builder = amoebaforces.AmoebaUreyBradleyForceBuilder()

    @staticmethod
    def parseElement(element, forceField):
        #  <AmoebaUreyBradleyForce>
        #   <UreyBradley class1="74" class2="73" class3="74" k="16003.8" d="0.15537" />

        generator = AmoebaUreyBradleyGenerator()
        forceField._forces.append(generator)
        for bond in element.findall('UreyBradley'):
            try:
                generator.builder.registerParams((bond.attrib['class1'], bond.attrib['class2'], bond.attrib['class3']), 
                                                 (float(bond.attrib['d']), float(bond.attrib['k'])))
            except:
                outputString = "AmoebaUreyBradleyGenerator : error getting types: %s %s %s" % (
                                    bond.attrib['class1'],
                                    bond.attrib['class2'],
                                    bond.attrib['class3'])
                raise ValueError(outputString)
        generator.classNameForType = dict((t.name, int(t.atomClass)) for t in forceField._atomTypes.values())

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = self.builder.getForce(sys)
        atomClasses = [str(self.classNameForType[data.atomType[atom]]) for atom in data.atoms]
        self.builder.addUreyBradleys(force, atomClasses, data.angles, data.isAngleConstrained, args.get('flexibleConstraints', False))

parsers["AmoebaUreyBradleyForce"] = AmoebaUreyBradleyGenerator.parseElement


## @private
class HippoNonbondedGenerator(object):
    """A HippoNonbondedGenerator constructs a HippoNonbondedForce."""

    def __init__(self, forcefield, extrapCoeff):
        self.ff = forcefield
        self.extrapCoeff = extrapCoeff
        self.exceptions = {}

    @staticmethod
    def parseElement(element, ff):
        extrapCoeff = [float(c) for c in element.attrib['extrapolationCoefficients'].split(',')]
        generator = HippoNonbondedGenerator(ff, extrapCoeff)
        ff.registerGenerator(generator)
        scaleNames = ('mmScale', 'dmScale', 'ddScale', 'dispScale', 'repScale', 'ctScale')
        paramNames = ('charge', 'coreCharge', 'alpha', 'epsilon', 'damping', 'c6', 'pauliK', 'pauliQ', 'pauliAlpha', 'polarizability', 'axisType', 'd0', 'd1', 'd2', 'q11', 'q12', 'q13', 'q21', 'q22', 'q23', 'q31', 'q32', 'q33')
        for ex in element.findall('Exception'):
            separation = int(ex.attrib['separation'])
            ingroup = ex.attrib['ingroup'].lower() == 'true'
            key = (separation, ingroup)
            if key in generator.exceptions:
                raise ValueError('HippoNonbondedForce: multiple exceptions with separation=%d ingroup=%s' % (separation, ingroup))
            generator.exceptions[key] = [float(ex.attrib[s]) for s in scaleNames]
        generator.params = ForceField._AtomTypeParameters(ff, 'HippoNonbondedForce', 'Atom', paramNames)
        generator.params.parseDefinitions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.HippoNonbondedForce.NoCutoff,
                     PME:mm.HippoNonbondedForce.PME}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for HippoNonbondedForce')

        # Build data structures we'll need for building local coordinate frames.

        bondIndices = _findBondsForExclusions(data, sys)
        pairs = _findExclusions(bondIndices, 2, len(data.atoms))
        bonded12 = [set() for i in range(len(data.atoms))]
        bonded13 = [set() for i in range(len(data.atoms))]
        for atom1, atom2, sep in pairs:
            if sep == 1:
                bonded12[atom1].add(data.atoms[atom2])
                bonded12[atom2].add(data.atoms[atom1])
            else:
                bonded13[atom1].add(data.atoms[atom2])
                bonded13[atom2].add(data.atoms[atom1])

        # Create the force.

        force = mm.HippoNonbondedForce()
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            params = [float(v) for v in values[:10]]
            axisType = int(values[10])
            dipole = [float(v) for v in values[11:14]]
            quadrupole = [float(v) for v in values[14:23]]
            extra = self.params.getExtraParameters(atom, data)
            zAtom = self._findAxisAtom('zAtomType', extra, bonded12[atom.index], None, data, [])
            xAtom = self._findAxisAtom('xAtomType', extra, bonded12[atom.index], bonded13[atom.index], data, [zAtom])
            yAtom = self._findAxisAtom('yAtomType', extra, bonded12[atom.index], bonded13[atom.index], data, [zAtom, xAtom])
            force.addParticle(params[0], dipole, quadrupole, *params[1:], axisType=axisType, multipoleAtomZ=zAtom, multipoleAtomX=xAtom, multipoleAtomY=yAtom)
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setExtrapolationCoefficients(self.extrapCoeff)
        force.setCutoffDistance(nonbondedCutoff)
        if args['switchDistance'] is not None:
            force.setSwitchingDistance(args['switchDistance'])
        if 'ewaldErrorTolerance' in args:
            force.setEwaldErrorTolerance(args['ewaldErrorTolerance'])
        sys.addForce(force)

    def _findAxisAtom(self, paramName, params, bonded12, bonded13, data, exclude):
        if paramName not in params:
            return -1
        atomType = params[paramName]
        for atom in bonded12:
            if data.atomType[atom] == atomType and atom.index not in exclude:
                return atom.index
        if bonded13 is not None:
            for atom in bonded13:
                if data.atomType[atom] == atomType and atom.index not in exclude:
                    return atom.index
        raise ValueError('No bonded atom of type %s' % atomType)

    def postprocessSystem(self, sys, data, args):
        # Identify polarization groups.

        bondIndices = _findBondsForExclusions(data, sys)
        groupBondTypes = [self.params.getExtraParameters(atom, data)['groupTypes'].split(',') for atom in data.atoms]
        groupBonds = [[] for i in range(len(data.atoms))]
        for i,j in bondIndices:
            if data.atomType[data.atoms[i]] in groupBondTypes[j]:
                groupBonds[i].append(j)
                groupBonds[j].append(i)
        polarizationGroup = _findGroups(groupBonds)

        # Create the exclusions.

        maxSeparation = max(e[0] for e in self.exceptions)
        hippo = [f for f in sys.getForces() if isinstance(f, mm.HippoNonbondedForce)][0]
        pairs = _findExclusions(bondIndices, maxSeparation, hippo.getNumParticles())
        for atom1, atom2, sep in pairs:
            params = self.exceptions[(sep, polarizationGroup[atom1] == polarizationGroup[atom2])]
            hippo.addException(atom1, atom2, *params)

parsers["HippoNonbondedForce"] = HippoNonbondedGenerator.parseElement

## @private
class DrudeGenerator(object):
    """A DrudeGenerator constructs a DrudeForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.typeMap = {}

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, DrudeGenerator)]
        if len(existing) == 0:
            generator = DrudeGenerator(ff)
            ff.registerGenerator(generator)
        else:
            # Multiple <DrudeForce> tags were found, probably in different files.  Simply add more types to the existing one.
            generator = existing[0]
        for particle in element.findall('Particle'):
            types = ff._findAtomTypes(particle.attrib, 5)
            if None not in types[:2]:
                aniso12 = 0.0
                aniso34 = 0.0
                if 'aniso12' in particle.attrib:
                    aniso12 = float(particle.attrib['aniso12'])
                if 'aniso34' in particle.attrib:
                    aniso34 = float(particle.attrib['aniso34'])
                values = (types[1], types[2], types[3], types[4], float(particle.attrib['charge']), float(particle.attrib['polarizability']), aniso12, aniso34, float(particle.attrib['thole']))
                for t in types[0]:
                    generator.typeMap[t] = values

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = mm.DrudeForce()
        if not any(isinstance(f, mm.NonbondedForce) for f in sys.getForces()):
            raise ValueError('<DrudeForce> must come after <NonbondedForce> in XML file')

        # Add Drude particles.

        for atom in data.atoms:
            t = data.atomType[atom]
            if t in self.typeMap:
                # Find other atoms in the residue that affect the Drude particle.
                p = [-1, -1, -1, -1]
                values = self.typeMap[t]
                for atom2 in atom.residue.atoms():
                    type2 = data.atomType[atom2]
                    if type2 in values[0]:
                        p[0] = atom2.index
                    elif values[1] is not None and type2 in values[1]:
                        p[1] = atom2.index
                    elif values[2] is not None and type2 in values[2]:
                        p[2] = atom2.index
                    elif values[3] is not None and type2 in values[3]:
                        p[3] = atom2.index
                force.addParticle(atom.index, p[0], p[1], p[2], p[3], values[4], values[5], values[6], values[7])
                data.excludeAtomWith[p[0]].append(atom.index)
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # For every nonbonded exclusion between Drude particles, add a screened pair.

        drude = [f for f in sys.getForces() if isinstance(f, mm.DrudeForce)][0]
        nonbonded = [f for f in sys.getForces() if isinstance(f, mm.NonbondedForce)][0]
        particleMap = {}
        for i in range(drude.getNumParticles()):
            particleMap[drude.getParticleParameters(i)[0]] = i
        for i in range(nonbonded.getNumExceptions()):
            (particle1, particle2, charge, sigma, epsilon) = nonbonded.getExceptionParameters(i)
            if charge._value == 0 and epsilon._value == 0:
                # This is an exclusion.
                if particle1 in particleMap and particle2 in particleMap:
                    # It connects two Drude particles, so add a screened pair.
                    drude1 = particleMap[particle1]
                    drude2 = particleMap[particle2]
                    type1 = data.atomType[data.atoms[particle1]]
                    type2 = data.atomType[data.atoms[particle2]]
                    thole1 = self.typeMap[type1][8]
                    thole2 = self.typeMap[type2][8]
                    drude.addScreenedPair(drude1, drude2, thole1+thole2)

        # Set the masses of Drude particles.

        drudeMass = args['drudeMass']
        if not unit.is_quantity(drudeMass):
            drudeMass *= unit.dalton
        for i in range(drude.getNumParticles()):
            params = drude.getParticleParameters(i)
            particle = params[0]
            parent = params[1]
            transferMass = drudeMass-sys.getParticleMass(particle)
            sys.setParticleMass(particle, drudeMass)
            sys.setParticleMass(parent, sys.getParticleMass(parent)-transferMass)

parsers["DrudeForce"] = DrudeGenerator.parseElement