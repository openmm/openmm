#!/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tools for constructing systems from AMBER prmtop/crd files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Randall J. Radmer, John D. Chodera, Peter Eastman
Contributors: Christoph Klein, Michael R. Shirts

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

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import os.path
import copy
import re
import math

try:
    import numpy
except:
    pass

import simtk.unit as units
import simtk.openmm
import customgbforces as customgb

#=============================================================================================
# AMBER parmtop loader (from 'zander', by Randall J. Radmer)
#=============================================================================================

# A regex for extracting print format info from the FORMAT lines.
FORMAT_RE_PATTERN=re.compile("([0-9]+)([a-zA-Z]+)([0-9]+)\.?([0-9]*)")

# Pointer labels which map to pointer numbers at top of prmtop files
POINTER_LABELS  = """
              NATOM,  NTYPES, NBONH,  MBONA,  NTHETH, MTHETA,
              NPHIH,  MPHIA,  NHPARM, NPARM,  NEXT,   NRES,
              NBONA,  NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA,
              NATYP,  NPHB,   IFPERT, NBPER,  NGPER,  NDPER,
              MBPER,  MGPER,  MDPER,  IFBOX,  NMXRS,  IFCAP
"""

# Pointer labels (above) as a list, not string.
POINTER_LABEL_LIST = POINTER_LABELS.replace(',', '').split()

class PrmtopLoader(object):
    """Parsed AMBER prmtop file.

    ParmtopLoader reads, parses and manages content from a AMBER prmtop file.

    EXAMPLES

    Parse a prmtop file of alanine dipeptide in implicit solvent.

    >>> import os, os.path
    >>> directory = os.path.join(os.getenv('YANK_INSTALL_DIR'), 'test', 'systems', 'alanine-dipeptide-gbsa')
    >>> prmtop_filename = os.path.join(directory, 'alanine-dipeptide.prmtop')
    >>> prmtop = PrmtopLoader(prmtop_filename)

    Parse a prmtop file of alanine dipeptide in explicit solvent.

    >>> import os, os.path
    >>> directory = os.path.join(os.getenv('YANK_INSTALL_DIR'), 'test', 'systems', 'alanine-dipeptide-explicit')
    >>> prmtop_filename = os.path.join(directory, 'alanine-dipeptide.prmtop')
    >>> prmtop = PrmtopLoader(prmtop_filename)    

    """
    def __init__(self, inFilename):
        """
        Create a PrmtopLoader object from an AMBER prmtop file.

        ARGUMENTS

        inFilename (string) - AMBER 'new-style' prmtop file, probably generated with one of the AMBER tleap/xleap/sleap        

        """

        self._prmtopVersion=None
        self._flags=[]
        self._raw_format={}
        self._raw_data={}

        fIn=open(inFilename)
        for line in fIn:
            if line.startswith('%VERSION'):
                tag, self._prmtopVersion = line.rstrip().split(None, 1)
            elif line.startswith('%FLAG'):
                tag, flag = line.rstrip().split(None, 1)
                self._flags.append(flag)
                self._raw_data[flag] = []
            elif line.startswith('%FORMAT'):
                format = line.rstrip()
                index0=format.index('(')
                index1=format.index(')')
                format = format[index0+1:index1]
                m = FORMAT_RE_PATTERN.search(format)
                self._raw_format[self._flags[-1]] = (format, m.group(1), m.group(2), m.group(3), m.group(4))
            elif self._flags \
                 and 'TITLE'==self._flags[-1] \
                 and not self._raw_data['TITLE']:
                self._raw_data['TITLE'] = line.rstrip()
            else:
                flag=self._flags[-1]
                (format, numItems, itemType,
                 itemLength, itemPrecision) = self._getFormat(flag)
                iLength=int(itemLength)
                line = line.rstrip()
                for index in range(0, len(line), iLength):
                    item = line[index:index+iLength]
                    if item:
                        self._raw_data[flag].append(item.strip())
        fIn.close()

    def _getFormat(self, flag=None):
        if not flag:
            flag=self._flags[-1]
        return self._raw_format[flag]

    def _getPointerValue(self, pointerLabel):
        """Return pointer value given pointer label

           Parameter:
            - pointerLabel: a string matching one of the following:

            NATOM  : total number of atoms 
            NTYPES : total number of distinct atom types
            NBONH  : number of bonds containing hydrogen
            MBONA  : number of bonds not containing hydrogen
            NTHETH : number of angles containing hydrogen
            MTHETA : number of angles not containing hydrogen
            NPHIH  : number of dihedrals containing hydrogen
            MPHIA  : number of dihedrals not containing hydrogen
            NHPARM : currently not used
            NPARM  : currently not used
            NEXT   : number of excluded atoms
            NRES   : number of residues
            NBONA  : MBONA + number of constraint bonds
            NTHETA : MTHETA + number of constraint angles
            NPHIA  : MPHIA + number of constraint dihedrals
            NUMBND : number of unique bond types
            NUMANG : number of unique angle types
            NPTRA  : number of unique dihedral types
            NATYP  : number of atom types in parameter file, see SOLTY below
            NPHB   : number of distinct 10-12 hydrogen bond pair types
            IFPERT : set to 1 if perturbation info is to be read in
            NBPER  : number of bonds to be perturbed
            NGPER  : number of angles to be perturbed
            NDPER  : number of dihedrals to be perturbed
            MBPER  : number of bonds with atoms completely in perturbed group
            MGPER  : number of angles with atoms completely in perturbed group
            MDPER  : number of dihedrals with atoms completely in perturbed groups
            IFBOX  : set to 1 if standard periodic box, 2 when truncated octahedral
            NMXRS  : number of atoms in the largest residue
            IFCAP  : set to 1 if the CAP option from edit was specified
        """
        index = POINTER_LABEL_LIST.index(pointerLabel) 
        return float(self._raw_data['POINTERS'][index])

    def getNumAtoms(self):
        """Return the number of atoms in the system"""
        return int(self._getPointerValue('NATOM'))

    def getNumTypes(self):
        """Return the number of AMBER atoms types in the system"""
        return int(self._getPointerValue('NTYPES'))

    def getIfBox(self):
        """Return True if the system was build with periodic boundary conditions (PBC)"""
        return int(self._getPointerValue('IFBOX'))

    def getIfCap(self):
        """Return True if the system was build with the cap option)"""
        return int(self._getPointerValue('IFCAP'))

    def getIfPert(self):
        """Return True if the system was build with the perturbation parameters)"""
        return int(self._getPointerValue('IFPERT'))

    def getMasses(self):
        """Return a list of atomic masses in the system"""
        try:
            return self._massList
        except AttributeError:
            pass

        self._massList=[]
        raw_masses=self._raw_data['MASS']
        for ii in range(self.getNumAtoms()):
            self._massList.append(float(raw_masses[ii]))
        self._massList = units.Quantity(self._massList, units.amu)
        return self._massList

    def getCharges(self):
        """Return a list of atomic charges in the system"""
        try:
            return self._chargeList
        except AttributeError:
            pass

        self._chargeList=[]
        raw_charges=self._raw_data['CHARGE']
        for ii in range(self.getNumAtoms()):
            self._chargeList.append(float(raw_charges[ii])/18.2223)
        self._chargeList = units.Quantity(self._chargeList, units.elementary_charge)
        return self._chargeList

    def getAtomName(self, iAtom):
        """Return the atom name for iAtom"""
        atomNames = self.getAtomNames()
        return atomNames[iAtom]

    def getAtomNames(self):
        """Return the list of the system atom names"""
        return self._raw_data['ATOM_NAME']

    def _getAtomTypeIndexes(self):
        try:
            return self._atomTypeIndexes
        except AttributeError:
            pass
        self._atomTypeIndexes=[]
        for atomTypeIndex in  self._raw_data['ATOM_TYPE_INDEX']:
            self._atomTypeIndexes.append(int(atomTypeIndex))
        return self._atomTypeIndexes

    def getAtomType(self, iAtom):
        """Return the AMBER atom type for iAtom"""
        atomTypes=self.getAtomTypes()
        return atomTypes[iAtom]

    def getAtomTypes(self):
        """Return the list of the AMBER atom types"""
        return self._raw_data['AMBER_ATOM_TYPE']

    def getResidueNumber(self, iAtom):
        """Return iAtom's residue number"""
        return self._getResiduePointer(iAtom)+1

    def getResidueLabel(self, iAtom=None, iRes=None):
        """Return residue label for iAtom OR iRes"""
        if iRes==None and iAtom==None:
            raise Exception("only specify iRes or iAtom, not both")
        if iRes!=None and iAtom!=None:
            raise Exception("iRes or iAtom must be set")
        if iRes!=None:
            return self._raw_data['RESIDUE_LABEL'][iRes]
        else:
            return self.getResidueLabel(iRes=self._getResiduePointer(iAtom))

    def _getResiduePointer(self, iAtom):
        try:
            return self.residuePointerDict[iAtom]
        except:
            pass
        self.residuePointerDict = {}
        resPointers=self._raw_data['RESIDUE_POINTER']
        firstAtom = [int(p)-1 for p in resPointers]
        firstAtom.append(self.getNumAtoms())
        res = 0
        for i in range(self.getNumAtoms()):
            while firstAtom[res+1] <= i:
                res += 1
            self.residuePointerDict[i] = res
        return self.residuePointerDict[iAtom]

    def getNonbondTerms(self):
        """Return list of all rVdw, epsilon pairs for each atom"""
        try:
            return self._nonbondTerms
        except AttributeError:
            pass
        self._nonbondTerms=[]
        for iAtom in range(self.getNumAtoms()):
            numTypes=self.getNumTypes()
            atomTypeIndexes=self._getAtomTypeIndexes()
            index=(numTypes+1)*(atomTypeIndexes[iAtom]-1)
            nbIndex=int(self._raw_data['NONBONDED_PARM_INDEX'][index])-1
            if nbIndex<0:
                raise Exception("10-12 interactions are not supported")
            acoef = float(self._raw_data['LENNARD_JONES_ACOEF'][nbIndex])
            bcoef = float(self._raw_data['LENNARD_JONES_BCOEF'][nbIndex])
            try:
                rMin = (2*acoef/bcoef)**(1/6.0)
                epsilon = 0.25*bcoef*bcoef/acoef
            except ZeroDivisionError:
                rMin = 1.0
                epsilon = 0.0
            rVdw =    units.Quantity(rMin/2.0, units.angstrom)
            epsilon = units.Quantity(epsilon, units.kilocalorie_per_mole)
            self._nonbondTerms.append( (rVdw, epsilon) )
        return self._nonbondTerms

    def _getBonds(self, bondPointers):
        forceConstant=self._raw_data["BOND_FORCE_CONSTANT"]
        bondEquil=self._raw_data["BOND_EQUIL_VALUE"]
        returnList=[]
        for ii in range(0,len(bondPointers),3):
             if int(bondPointers[ii])<0 or \
                int(bondPointers[ii+1])<0:
                 raise Exception("Found negative bonded atom pointers %s"
                                 % ((bondPointers[ii],
                                     bondPointers[ii+1]),))
             iType=int(bondPointers[ii+2])-1
             forceConstUnit=units.kilocalorie_per_mole/(units.angstrom*units.angstrom)
             returnList.append((int(bondPointers[ii])/3,
                                int(bondPointers[ii+1])/3,
                                units.Quantity(float(forceConstant[iType]), forceConstUnit),
                                units.Quantity(float(bondEquil[iType]), units.angstrom)))
        return returnList

    def getBondsWithH(self):
        """Return list of bonded atom pairs, K, and Rmin for each bond with a hydrogen"""
        try:
            return self._bondListWithH
        except AttributeError:
            pass
        bondPointers=self._raw_data["BONDS_INC_HYDROGEN"]
        self._bondListWithH = self._getBonds(bondPointers)
        return self._bondListWithH
        

    def getBondsNoH(self):
        """Return list of bonded atom pairs, K, and Rmin for each bond with no hydrogen"""
        try:
            return self._bondListNoH
        except AttributeError:
            pass
        bondPointers=self._raw_data["BONDS_WITHOUT_HYDROGEN"]
        self._bondListNoH = self._getBonds(bondPointers)
        return self._bondListNoH

    def getAngles(self):
        """Return list of atom triplets, K, and ThetaMin for each bond angle"""
        try:
            return self._angleList
        except AttributeError:
            pass
        forceConstant=self._raw_data["ANGLE_FORCE_CONSTANT"]
        angleEquil=self._raw_data["ANGLE_EQUIL_VALUE"]
        anglePointers = self._raw_data["ANGLES_INC_HYDROGEN"] \
                       +self._raw_data["ANGLES_WITHOUT_HYDROGEN"]
        self._angleList=[]
        for ii in range(0,len(anglePointers),4):
             if int(anglePointers[ii])<0 or \
                int(anglePointers[ii+1])<0 or \
                int(anglePointers[ii+2])<0:
                 raise Exception("Found negative angle atom pointers %s"
                                 % ((anglePointers[ii],
                                     anglePointers[ii+1],
                                     anglePointers[ii+2]),))
             iType=int(anglePointers[ii+3])-1
             forceConstUnit=units.kilocalorie_per_mole/(units.radian*units.radian)
             self._angleList.append((int(anglePointers[ii])/3,
                                int(anglePointers[ii+1])/3,
                                int(anglePointers[ii+2])/3,
                                units.Quantity(float(forceConstant[iType]), forceConstUnit),
                                units.Quantity(float(angleEquil[iType])*180/math.pi, units.degree)))
        return self._angleList

    def getDihedrals(self):
        """Return list of atom quads, K, phase and periodicity for each dihedral angle"""
        try:
            return self._dihedralList
        except AttributeError:
            pass
        forceConstant=self._raw_data["DIHEDRAL_FORCE_CONSTANT"]
        phase=self._raw_data["DIHEDRAL_PHASE"]
        periodicity=self._raw_data["DIHEDRAL_PERIODICITY"]
        dihedralPointers = self._raw_data["DIHEDRALS_INC_HYDROGEN"] \
                          +self._raw_data["DIHEDRALS_WITHOUT_HYDROGEN"]
        self._dihedralList=[]
        for ii in range(0,len(dihedralPointers),5):
             if int(dihedralPointers[ii])<0 or int(dihedralPointers[ii+1])<0:
                 raise Exception("Found negative dihedral atom pointers %s"
                                 % ((dihedralPointers[ii],
                                    dihedralPointers[ii+1],
                                    dihedralPointers[ii+2],
                                    dihedralPointers[ii+3]),))
             iType=int(dihedralPointers[ii+4])-1
             self._dihedralList.append((int(dihedralPointers[ii])/3,
                                int(dihedralPointers[ii+1])/3,
                                abs(int(dihedralPointers[ii+2]))/3,
                                abs(int(dihedralPointers[ii+3]))/3,
                                units.Quantity(float(forceConstant[iType]), units.kilocalorie_per_mole),
                                units.Quantity(float(phase[iType])*180/math.pi, units.degree),
                                int(0.5+float(periodicity[iType]))))
        return self._dihedralList

    def get14Interactions(self):
        """Return list of atom pairs, chargeProduct, rMin and epsilon for each 1-4 interaction"""
        dihedralPointers = self._raw_data["DIHEDRALS_INC_HYDROGEN"] \
                          +self._raw_data["DIHEDRALS_WITHOUT_HYDROGEN"]
        returnList=[]
        charges=self.getCharges()
        nonbondTerms = self.getNonbondTerms()
        for ii in range(0,len(dihedralPointers),5):
             if int(dihedralPointers[ii+2])>0 and int(dihedralPointers[ii+3])>0:
                 iAtom = int(dihedralPointers[ii])/3
                 lAtom = int(dihedralPointers[ii+3])/3
                 chargeProd = charges[iAtom]*charges[lAtom]
                 (rVdwI, epsilonI) = nonbondTerms[iAtom]
                 (rVdwL, epsilonL) = nonbondTerms[lAtom]
                 rMin = (rVdwI+rVdwL)
                 epsilon = units.sqrt(epsilonI*epsilonL)
                 returnList.append((iAtom, lAtom, chargeProd, rMin, epsilon))
        return returnList

    def getExcludedAtoms(self):
        """Return list of lists, giving all pairs of atoms that should have no non-bond interactions"""
        try:
            return self._excludedAtoms
        except AttributeError:
            pass
        self._excludedAtoms=[]
        numExcludedAtomsList=self._raw_data["NUMBER_EXCLUDED_ATOMS"]
        excludedAtomsList=self._raw_data["EXCLUDED_ATOMS_LIST"]
        total=0
        for iAtom in range(self.getNumAtoms()):
            index0=total
            n=int(numExcludedAtomsList[iAtom])
            total+=n
            index1=total
            atomList=[]
            for jAtom in excludedAtomsList[index0:index1]:
                j=int(jAtom)
                if j>0:
                    atomList.append(j-1)
            self._excludedAtoms.append(atomList)
        return self._excludedAtoms

    def getGBParms(self, symbls=None):
        """Return list giving GB params, Radius and screening factor"""
        try:
            return self._gb_List
        except AttributeError:
            pass
        self._gb_List=[]
        radii=self._raw_data["RADII"]
        screen=self._raw_data["SCREEN"]
        # Update screening parameters for GBn if specified
        if symbls:
            for (i, symbl) in enumerate(symbls):
                if symbl[0] == ('c' or 'C'):
                    screen[i] = 0.48435382330
                elif symbl[0] == ('h' or 'H'):
                    screen[i] = 1.09085413633
                elif symbl[0] == ('n' or 'N'):
                    screen[i] = 0.700147318409
                elif symbl[0] == ('o' or 'O'):
                    screen[i] = 1.06557401132
                elif symbl[0] == ('s' or 'S'):
                    screen[i] = 0.602256336067
                else:
                    screen[i] = 0.5
        for iAtom in range(len(radii)):
            self._gb_List.append((float(radii[iAtom])*units.angstrom, float(screen[iAtom])))
        return self._gb_List

    def getBoxBetaAndDimensions(self):
        """Return periodic boundary box beta angle and dimensions"""
        beta=float(self._raw_data["BOX_DIMENSIONS"][0])
        x=float(self._raw_data["BOX_DIMENSIONS"][1])
        y=float(self._raw_data["BOX_DIMENSIONS"][2])
        z=float(self._raw_data["BOX_DIMENSIONS"][3])
        return (units.Quantity(beta, units.degree),
                units.Quantity(x, units.angstrom),
                units.Quantity(y, units.angstrom),
                units.Quantity(z, units.angstrom))

#=============================================================================================
# AMBER System builder (based on, but not identical to, systemManager from 'zander')
#=============================================================================================

def readAmberSystem(prmtop_filename=None, prmtop_loader=None, shake=None, gbmodel=None, soluteDielectric=1.0, solventDielectric=78.5, nonbondedCutoff=None, nonbondedMethod='NoCutoff', scee=1.2, scnb=2.0, mm=None, verbose=False, EwaldErrorTolerance=None, flexibleConstraints=True, rigidWater=True):
    """
    Create an OpenMM System from an Amber prmtop file.
    
    ARGUMENTS (specify  one or the other, but not both)
      prmtop_filename (String) - name of Amber prmtop file (new-style only)
      prmtop_loader (PrmtopLoader) - the loaded prmtop file
      
    OPTIONAL ARGUMENTS
      shake (String) - if 'h-bonds', will SHAKE all bonds to hydrogen and water; if 'all-bonds', will SHAKE all bonds and water (default: None)
      gbmodel (String) - if 'OBC', OBC GBSA will be used; if 'GBVI', GB/VI will be used (default: None)
      soluteDielectric (float) - The solute dielectric constant to use in the implicit solvent model (default: 1.0)
      solventDielectric (float) - The solvent dielectric constant to use in the implicit solvent model (default: 78.5)
      nonbondedCutoff (float) - if specified, will set nonbondedCutoff (default: None)
      scnb (float) - 1-4 Lennard-Jones scaling factor (default: 1.2)
      scee (float) - 1-4 electrostatics scaling factor (default: 2.0)
      mm - if specified, this module will be used in place of pyopenmm (default: None)
      verbose (boolean) - if True, print out information on progress (default: False)
      flexibleConstraints (boolean) - if True, flexible bonds will be added in addition ot constrained bonds
      rigidWater (boolean=True) If true, water molecules will be fully rigid regardless of the value passed for the shake argument

    NOTES

    Even if bonds are SHAKEn, their harmonic stretch terms are still included in the potential.

    TODO

    Should these option names be changed to reflect their 'sander' counterparts?

    EXAMPLES

    Create a system of alanine dipeptide in implicit solvent.

    >>> directory = os.path.join(os.getenv('YANK_INSTALL_DIR'), 'test', 'systems', 'alanine-dipeptide-gbsa')
    >>> prmtop_filename = os.path.join(directory, 'alanine-dipeptide.prmtop')
    >>> system = readAmberSystem(prmtop_filename)

    Parse a prmtop file of alanine dipeptide in explicit solvent.

    >>> directory = os.path.join(os.getenv('YANK_INSTALL_DIR'), 'test', 'systems', 'alanine-dipeptide-explicit')
    >>> prmtop_filename = os.path.join(directory, 'alanine-dipeptide.prmtop')
    >>> system = readAmberSystem(prmtop_filename)    

    """
    
    if prmtop_filename is None and prmtop_loader is None:
        raise Exception("Must specify a filename or loader")
    if prmtop_filename is not None and prmtop_loader is not None:
        raise Exception("Cannot specify both a filename and a loader")
    if prmtop_filename is not None:
        # Load prmtop file.
        if verbose: print "Reading prmtop file '%s'..." % prmtop_filename
        prmtop = PrmtopLoader(prmtop_filename)
    else:
        prmtop = prmtop_loader

    if prmtop.getIfCap()>0:
        raise Exception("CAP option not currently supported")

    if prmtop.getIfPert()>0:
        raise Exception("perturbation not currently supported")
        
    if prmtop.getIfBox()>1:
        raise Exception("only standard periodic boxes are currently supported")

    # Use pyopenmm implementation of OpenMM by default.
    if mm is None:
        mm = simtk.openmm

    # Create OpenMM System.
    if verbose: print "Creating OpenMM system..."
    system = mm.System()

    # Populate system with atomic masses.
    if verbose: print "Adding particles..."
    for mass in prmtop.getMasses():
        system.addParticle(mass)

    # Add constraints.
    isWater = [prmtop.getResidueLabel(i) == 'WAT' for i in range(prmtop.getNumAtoms())]
    if shake in ('h-bonds', 'all-bonds', 'h-angles'):
        for (iAtom, jAtom, k, rMin) in prmtop.getBondsWithH():
            system.addConstraint(iAtom, jAtom, rMin)
    if shake in ('all-bonds', 'h-angles'):
        for (iAtom, jAtom, k, rMin) in prmtop.getBondsNoH():
            system.addConstraint(iAtom, jAtom, rMin)
    if rigidWater and shake == None:
        for (iAtom, jAtom, k, rMin) in prmtop.getBondsWithH():
            if isWater[iAtom] and isWater[jAtom]:
                system.addConstraint(iAtom, jAtom, rMin)
            
    # Add harmonic bonds.
    if verbose: print "Adding bonds..."    
    force = mm.HarmonicBondForce()
    if flexibleConstraints or (shake not in ('h-bonds', 'all-bonds', 'h-angles')):
        for (iAtom, jAtom, k, rMin) in prmtop.getBondsWithH():
            if flexibleConstraints or not (rigidWater and isWater[iAtom] and isWater[jAtom]):
                force.addBond(iAtom, jAtom, rMin, 2*k)
    if flexibleConstraints or (shake not in ('all-bonds', 'h-angles')):
        for (iAtom, jAtom, k, rMin) in prmtop.getBondsNoH():
            force.addBond(iAtom, jAtom, rMin, 2*k)
    system.addForce(force)

    # Add harmonic angles.
    if verbose: print "Adding angles..."    
    force = mm.HarmonicAngleForce()
    if shake == 'h-angles':
        numConstrainedBonds = system.getNumConstraints()
        atomConstraints = [[]]*system.getNumParticles()
        for i in range(system.getNumConstraints()):
            c = system.getConstraintParameters(i)
            atomConstraints[c[0]].append((c[1], c[2]))
            atomConstraints[c[1]].append((c[0], c[2]))
    for (iAtom, jAtom, kAtom, k, aMin) in prmtop.getAngles():
        if shake == 'h-angles':
            type1 = prmtop.getAtomType(iAtom)
            type2 = prmtop.getAtomType(jAtom)
            type3 = prmtop.getAtomType(kAtom)
            numH = len([type for type in (type1, type3) if type.startswith('H')])
            constrained = (numH == 2 or (numH == 1 and type2.startswith('O')))
        else:
            constrained = False
        if constrained:
            # Find the two bonds that make this angle.
            l1 = None
            l2 = None
            for bond in atomConstraints[jAtom]:
                if bond[0] == iAtom:
                    l1 = bond[1]
                elif bond[0] == kAtom:
                    l2 = bond[1]
            
            # Compute the distance between atoms and add a constraint
            length = units.sqrt(l1*l1 + l2*l2 - 2*l1*l2*units.cos(aMin))
            system.addConstraint(iAtom, kAtom, length)
        if flexibleConstraints or not constrained:
            force.addAngle(iAtom, jAtom, kAtom, aMin, 2*k)
    system.addForce(force)

    # Add torsions.
    if verbose: print "Adding torsions..."    
    force = mm.PeriodicTorsionForce()
    for (iAtom, jAtom, kAtom, lAtom, forceConstant, phase, periodicity) in prmtop.getDihedrals():
        force.addTorsion(iAtom, jAtom, kAtom, lAtom, periodicity, phase, forceConstant)
    system.addForce(force)

    # Add nonbonded interactions.
    if verbose: print "Adding nonbonded interactions..."    
    force = mm.NonbondedForce()
    if (prmtop.getIfBox() == 0):
        # System is non-periodic.
        if nonbondedMethod == 'NoCutoff':
            force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        elif nonbondedMethod == 'CutoffNonPeriodic':
            if nonbondedCutoff is None:
                raise Exception("No cutoff value specified")
            force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)            
            force.setCutoffDistance(nonbondedCutoff)
        else:
            raise Exception("Illegal nonbonded method for a non-periodic system")
    else:
        # System is periodic. 
        # Set periodic box vectors for periodic system
        (boxBeta, boxX, boxY, boxZ) = prmtop.getBoxBetaAndDimensions()
        d0 = units.Quantity(0.0, units.angstroms)
        xVec = units.Quantity((boxX, d0,   d0))
        yVec = units.Quantity((d0,   boxY, d0))
        zVec = units.Quantity((d0,   d0,   boxZ))
        system.setDefaultPeriodicBoxVectors(xVec, yVec, zVec)
        
        # Set cutoff.
        if nonbondedCutoff is None:
            # Compute cutoff automatically.
            min_box_width = min([boxX / units.nanometers, boxY / units.nanometers, boxZ / units.nanometers])
            CLEARANCE_FACTOR = 0.97 # reduce the cutoff to be a bit smaller than 1/2 smallest box length            
            nonbondedCutoff = units.Quantity((min_box_width * CLEARANCE_FACTOR) / 2.0, units.nanometers)
        if nonbondedMethod != 'NoCutoff':
            force.setCutoffDistance(nonbondedCutoff)
        
        # Set nonbonded method.
        if nonbondedMethod == 'NoCutoff':
            force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        elif nonbondedMethod == 'CutoffNonPeriodic':
            force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
        elif nonbondedMethod == 'CutoffPeriodic':
            force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        elif nonbondedMethod == 'Ewald':
            force.setNonbondedMethod(mm.NonbondedForce.Ewald)
        elif nonbondedMethod == 'PME':
            force.setNonbondedMethod(mm.NonbondedForce.PME)
        else:
            raise Exception("Cutoff method not understood.")

        if EwaldErrorTolerance is not None:
            force.setEwaldErrorTolerance(EwaldErrorTolerance)

    # Add per-particle nonbonded parameters.
    sigmaScale = 2**(-1./6.) * 2.0
    for (charge, (rVdw, epsilon)) in zip(prmtop.getCharges(), prmtop.getNonbondTerms()):
        sigma = rVdw * sigmaScale
        force.addParticle(charge, sigma, epsilon)

    # Add 1-4 Interactions
    excludedAtomPairs = set()
    sigmaScale = 2**(-1./6.)
    for (iAtom, lAtom, chargeProd, rMin, epsilon) in prmtop.get14Interactions():
        chargeProd /= scee
        epsilon /= scnb
        sigma = rMin * sigmaScale
        force.addException(iAtom, lAtom, chargeProd, sigma, epsilon)
        excludedAtomPairs.add(min((iAtom, lAtom), (lAtom, iAtom)))

    # Add Excluded Atoms
    excludedAtoms=prmtop.getExcludedAtoms()
    excludeParams = (0.0*units.elementary_charge**2, 1.0*units.angstroms, 0.0*units.kilocalories_per_mole)
    for iAtom in range(prmtop.getNumAtoms()):
        for jAtom in excludedAtoms[iAtom]:
            if min((iAtom, jAtom), (jAtom, iAtom)) in excludedAtomPairs: continue            
            force.addException(iAtom, jAtom, excludeParams[0], excludeParams[1], excludeParams[2])

    system.addForce(force)

    # Add virtual sites for water.
    epNames = ['EP', 'LP']
    ep = [i for i in range(prmtop.getNumAtoms()) if isWater[i] and prmtop.getAtomName(i)[:2] in epNames]
    if len(ep) > 0:
        epRes = set((prmtop.getResidueNumber(i) for i in ep))
        numRes = max(epRes)+1
        # For each residue that contains an "extra point", find the oxygen, hydrogens, and points.
        waterO = []
        waterH = []
        waterEP = []
        for i in range(numRes):
            waterO.append([])
            waterH.append([])
            waterEP.append([])
        for i in range(prmtop.getNumAtoms()):
            res = prmtop.getResidueNumber(i)
            if res in epRes:
                name = prmtop.getAtomName(i)
                if name[0] == 'O':
                    waterO[res].append(i)
                if name[0] == 'H':
                    waterH[res].append(i)
                if name[:2] in epNames:
                    waterEP[res].append(i)
        # Record bond lengths for faster access.
        distOH = [None]*numRes
        distHH = [None]*numRes
        distOE = [None]*numRes
        for (atom1, atom2, k, dist) in prmtop.getBondsWithH()+prmtop.getBondsNoH():
            res = prmtop.getResidueNumber(atom1)
            if res in epRes:
                name1 = prmtop.getAtomName(atom1)
                name2 = prmtop.getAtomName(atom2)
                if name1[0] == 'H' or name2[0] == 'H':
                    if name1[0] == 'H' and name2[0] == 'H':
                        distHH[res] = dist
                    if name1[0] == 'O' or name2[0] == 'O':
                        distOH[res] = dist
                elif (name1[0] == 'O' or name2[0] == 'O') and ((name1[:2] in epNames or name2[:2] in epNames)):
                    distOE[res] = dist
        # Loop over residues and add the virtual sites.
        outOfPlaneAngle = 54.735*units.degree
        cosOOP = units.cos(outOfPlaneAngle)
        sinOOP = units.sin(outOfPlaneAngle)
        for res in range(numRes):
            if len(waterO[res]) == 1 and len(waterH[res]) == 2:
                if len(waterEP[res]) == 1:
                    # Four point water
                    weightH = distOE[res]/units.sqrt(distOH[res]**2-(0.5*distHH[res])**2)
                    system.setVirtualSite(waterEP[res][0], mm.ThreeParticleAverageSite(waterO[res][0], waterH[res][0], waterH[res][1], 1-weightH, weightH/2, weightH/2))
                elif len(waterEP[res]) == 2:
                    # Five point water
                    weightH = cosOOP*distOE[res]/units.sqrt(distOH[res]**2-(0.5*distHH[res])**2)
                    angleHOH = 2*math.asin(0.5*distHH[res]/distOH[res])
                    lenCross = (distOH[res].value_in_unit(units.nanometer)**2)*math.sin(angleHOH)*units.nanometer
                    weightCross = sinOOP*distOE[res]/lenCross
                    system.setVirtualSite(waterEP[res][0], mm.OutOfPlaneSite(waterO[res][0], waterH[res][0], waterH[res][1], weightH/2, weightH/2, weightCross))
                    system.setVirtualSite(waterEP[res][1], mm.OutOfPlaneSite(waterO[res][0], waterH[res][0], waterH[res][1], weightH/2, weightH/2, -weightCross))

    # Add GBSA model.
    if gbmodel is not None:
        if verbose: print "Adding GB parameters..."
        charges = prmtop.getCharges()
        symbls = None
        if gbmodel == 'GBn':
            symbls = prmtop.getAtomTypes()
        gb_parms = prmtop.getGBParms(symbls)
        if gbmodel == 'HCT':
            gb = customgb.GBSAHCTForce(solventDielectric, soluteDielectric, 'ACE')
        elif gbmodel == 'OBC1':
            gb = customgb.GBSAOBC1Force(solventDielectric, soluteDielectric, 'ACE')
        elif gbmodel == 'OBC2':
            gb = mm.GBSAOBCForce()
            gb.setSoluteDielectric(soluteDielectric)
            gb.setSolventDielectric(solventDielectric)
        elif gbmodel == 'GBn':
            gb = customgb.GBSAGBnForce(solventDielectric, soluteDielectric, 'ACE')
        else:
            raise Exception("Illegal value specified for implicit solvent model")
        for iAtom in range(prmtop.getNumAtoms()):
            if gbmodel == 'OBC2':
                gb.addParticle(charges[iAtom], gb_parms[iAtom][0], gb_parms[iAtom][1])
            else:
                gb.addParticle([charges[iAtom], gb_parms[iAtom][0], gb_parms[iAtom][1]])        
        system.addForce(gb)
        if nonbondedMethod == 'NoCutoff':
            gb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        elif nonbondedMethod == 'CutoffNonPeriodic':
            gb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
            gb.setCutoffDistance(nonbondedCutoff)
        elif nonbondedMethod == 'CutoffPeriodic':
            gb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            gb.setCutoffDistance(nonbondedCutoff)
        else:
            raise Exception("Illegal nonbonded method for use with GBSA")

    # TODO: Add GBVI terms?

    return system

#=============================================================================================
# AMBER INPCRD loader
#=============================================================================================

def readAmberCoordinates(filename, read_box=False, read_velocities=False, verbose=False, asNumpy=False):
    """
    Read atomic coordinates (and optionally, box vectors) from Amber formatted coordinate file.    

    ARGUMENTS

    filename (string) - name of Amber coordinates file to be read in
    system (simtk.openmm.System) - System object for which coordinates are to be read

    OPTIONAL ARGUMENTS

    verbose (boolean) - if True, will print out verbose information about the file being read
    asNumpy (boolean) - if True, results will be returned as Numpy arrays instead of lists of Vec3s

    EXAMPLES

    Read coordinates in vacuum.

    >>> directory = os.path.join(os.getenv('YANK_INSTALL_DIR'), 'test', 'systems', 'alanine-dipeptide-gbsa')
    >>> crd_filename = os.path.join(directory, 'alanine-dipeptide.inpcrd')    
    >>> coordinates = readAmberCoordinates(crd_filename)

    Read coordinates in solvent.

    >>> directory = os.path.join(os.getenv('YANK_INSTALL_DIR'), 'test', 'systems', 'alanine-dipeptide-explicit')
    >>> crd_filename = os.path.join(directory, 'alanine-dipeptide.inpcrd')    
    >>> [coordinates, box_vectors] = readAmberCoordinates(crd_filename, read_box=True)

    """

    # Open coordinate file for reading.
    infile = open(filename, 'r')

    # Read title
    title = infile.readline().strip()
    if verbose: print "title: '%s'" % title

    # Read number of atoms
    natoms = int(infile.readline().split()[0])
    if verbose: print "%d atoms" % natoms

    # Allocate storage for coordinates
    coordinates = []

    # Read coordinates
    mm = simtk.openmm
    natoms_read = 0
    while (natoms_read < natoms):
        line = infile.readline()
        if len(line) == 0:
            raise ValueError("Unexpected end of file while reading coordinates")
        line = line.strip()
        elements = line.split()
        while (len(elements) > 0):
            coordinates.append(mm.Vec3(float(elements.pop(0)), float(elements.pop(0)), float(elements.pop(0))))
            natoms_read += 1
    if asNumpy:
        newcoords = numpy.zeros([natoms,3], numpy.float32)
        for i in range(len(coordinates)):
            for j in range(3):
                newcoords[i,j] = coordinates[i][j]
        coordinates = newcoords
    # Assign units.
    coordinates = units.Quantity(coordinates, units.angstroms)

    # Read velocities if requested.
    velocities = None
    if (read_velocities):
        # Read velocities
        velocities = []
        natoms_read = 0
        while (natoms_read < natoms):
            line = infile.readline()
            if len(line) == 0:
                raise ValueError("Unexpected end of file while reading velocities")
            line = line.strip()
            elements = line.split()
            while (len(elements) > 0):
                velocities.append(20.455*mm.Vec3(float(elements.pop(0)), float(elements.pop(0)), float(elements.pop(0))))
                natoms_read += 1
        if asNumpy:
            newvel = numpy.zeros([natoms,3], numpy.float32)
            for i in range(len(velocities)):
                for j in range(3):
                    newvel[i,j] = velocities[i][j]
            velocities = newvel
        # Assign units.
        velocities = units.Quantity(velocities, units.angstroms/units.picoseconds)
            
    # Read box size if present
    box_vectors = None
    if (read_box):
        line = infile.readline()
        if len(line) == 0:
            raise ValueError("Unexpected end of file while reading box vectors")
        line = line.strip()
        elements = line.split()
        nelements = len(elements)
        box_dimensions = [0.0]*nelements
        for i in range(nelements):
            box_dimensions[i] = float(elements[i])
        # TODO: Deal with non-standard box sizes.
        if nelements == 6:
            if asNumpy:
                a = units.Quantity(numpy.array([box_dimensions[0], 0.0, 0.0]), units.angstroms)
                b = units.Quantity(numpy.array([0.0, box_dimensions[1], 0.0]), units.angstroms)
                c = units.Quantity(numpy.array([0.0, 0.0, box_dimensions[2]]), units.angstroms)
            else:
                a = units.Quantity(mm.Vec3(box_dimensions[0], 0.0, 0.0), units.angstroms)
                b = units.Quantity(mm.Vec3(0.0, box_dimensions[1], 0.0), units.angstroms)
                c = units.Quantity(mm.Vec3(0.0, 0.0, box_dimensions[2]), units.angstroms)            
            box_vectors = [a,b,c]
        else:
            raise Exception("Don't know what to do with box vectors: %s" % line)

    # Close file
    infile.close()
    
    if box_vectors and velocities:
        return (coordinates, box_vectors, velocities)
    if box_vectors:
        return (coordinates, box_vectors)
    if velocities:
        return (coordinates, velocities)    
    return coordinates

#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

if __name__ == "__main__":
    import doctest
    doctest.testmod()

    
