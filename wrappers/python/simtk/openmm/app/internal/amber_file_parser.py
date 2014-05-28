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

Portions copyright (c) 2012-2013 Stanford University and the Authors.
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
import re
from math import ceil, cos, sin, asin, sqrt, pi
import warnings

try:
    import numpy as np
except:
    np = None

import simtk.unit as units
import simtk.openmm
from simtk.openmm.app import element as elem
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

VELSCALE = 20.455 # velocity conversion factor to angstroms/picosecond
TINY = 1.0e-8

class PrmtopLoader(object):
    """Parsed AMBER prmtop file.

    ParmtopLoader reads, parses and manages content from a AMBER prmtop file.

    EXAMPLES

    Parse a prmtop file of alanine dipeptide in implicit solvent.

    >>> import os
    >>> directory = os.path.join(os.getenv('YANK_INSTALL_DIR'), 'test', 'systems', 'alanine-dipeptide-gbsa')
    >>> prmtop_filename = os.path.join(directory, 'alanine-dipeptide.prmtop')
    >>> prmtop = PrmtopLoader(prmtop_filename)

    Parse a prmtop file of alanine dipeptide in explicit solvent.

    >>> import os
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
        self._massList = self._massList
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
        self._chargeList = self._chargeList
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
        lengthConversionFactor = units.angstrom.conversion_factor_to(units.nanometer)
        energyConversionFactor = units.kilocalorie_per_mole.conversion_factor_to(units.kilojoule_per_mole)
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
            rVdw = rMin/2.0*lengthConversionFactor
            epsilon = epsilon*energyConversionFactor
            self._nonbondTerms.append( (rVdw, epsilon) )
        return self._nonbondTerms

    def _getBonds(self, bondPointers):
        forceConstant=self._raw_data["BOND_FORCE_CONSTANT"]
        bondEquil=self._raw_data["BOND_EQUIL_VALUE"]
        returnList=[]
        forceConstConversionFactor = (units.kilocalorie_per_mole/(units.angstrom*units.angstrom)).conversion_factor_to(units.kilojoule_per_mole/(units.nanometer*units.nanometer))
        lengthConversionFactor = units.angstrom.conversion_factor_to(units.nanometer)
        for ii in range(0,len(bondPointers),3):
             if int(bondPointers[ii])<0 or \
                int(bondPointers[ii+1])<0:
                 raise Exception("Found negative bonded atom pointers %s"
                                 % ((bondPointers[ii],
                                     bondPointers[ii+1]),))
             iType=int(bondPointers[ii+2])-1
             returnList.append((int(bondPointers[ii])//3,
                                int(bondPointers[ii+1])//3,
                                float(forceConstant[iType])*forceConstConversionFactor,
                                float(bondEquil[iType])*lengthConversionFactor))
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
        forceConstConversionFactor = (units.kilocalorie_per_mole/(units.radian*units.radian)).conversion_factor_to(units.kilojoule_per_mole/(units.radian*units.radian))
        for ii in range(0,len(anglePointers),4):
             if int(anglePointers[ii])<0 or \
                int(anglePointers[ii+1])<0 or \
                int(anglePointers[ii+2])<0:
                 raise Exception("Found negative angle atom pointers %s"
                                 % ((anglePointers[ii],
                                     anglePointers[ii+1],
                                     anglePointers[ii+2]),))
             iType=int(anglePointers[ii+3])-1
             self._angleList.append((int(anglePointers[ii])//3,
                                int(anglePointers[ii+1])//3,
                                int(anglePointers[ii+2])//3,
                                float(forceConstant[iType])*forceConstConversionFactor,
                                float(angleEquil[iType])))
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
        forceConstConversionFactor = (units.kilocalorie_per_mole).conversion_factor_to(units.kilojoule_per_mole)
        for ii in range(0,len(dihedralPointers),5):
             if int(dihedralPointers[ii])<0 or int(dihedralPointers[ii+1])<0:
                 raise Exception("Found negative dihedral atom pointers %s"
                                 % ((dihedralPointers[ii],
                                    dihedralPointers[ii+1],
                                    dihedralPointers[ii+2],
                                    dihedralPointers[ii+3]),))
             iType=int(dihedralPointers[ii+4])-1
             self._dihedralList.append((int(dihedralPointers[ii])//3,
                                int(dihedralPointers[ii+1])//3,
                                abs(int(dihedralPointers[ii+2]))//3,
                                abs(int(dihedralPointers[ii+3]))//3,
                                float(forceConstant[iType])*forceConstConversionFactor,
                                float(phase[iType]),
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
                iAtom = int(dihedralPointers[ii])//3
                lAtom = int(dihedralPointers[ii+3])//3
                iidx = int(dihedralPointers[ii+4]) - 1
                chargeProd = charges[iAtom]*charges[lAtom]
                (rVdwI, epsilonI) = nonbondTerms[iAtom]
                (rVdwL, epsilonL) = nonbondTerms[lAtom]
                rMin = (rVdwI+rVdwL)
                epsilon = sqrt(epsilonI*epsilonL)
                try:
                    iScee = float(self._raw_data["SCEE_SCALE_FACTOR"][iidx])
                except KeyError:
                    iScee = 1.2
                try:
                    iScnb = float(self._raw_data["SCNB_SCALE_FACTOR"][iidx])
                except KeyError:
                    iScnb = 2.0

                returnList.append((iAtom, lAtom, chargeProd, rMin, epsilon, iScee, iScnb))
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

    def getGBParms(self, gbmodel, elements):
        """Return list giving GB params, Radius and screening factor"""
        gb_List=[]
        radii=[float(r) for r in self._raw_data["RADII"]]
        screen=[float(s) for s in self._raw_data["SCREEN"]]
        if gbmodel == 'GBn2':
            alpha = [0 for i in self._raw_data['RADII']]
            beta = [0 for i in self._raw_data['RADII']]
            gamma = [0 for i in self._raw_data['RADII']]
        # Update screening parameters for GBn if specified
        if gbmodel == 'GBn':
            if elements is None:
                raise Exception('GBn model requires element information')
            for i, element in enumerate(elements):
                if element is elem.carbon:
                    screen[i] = 0.48435382330
                elif element is elem.hydrogen:
                    screen[i] = 1.09085413633
                elif element is elem.nitrogen:
                    screen[i] = 0.700147318409
                elif element is elem.oxygen:
                    screen[i] = 1.06557401132
                elif element is elem.sulfur:
                    screen[i] = 0.602256336067
                else:
                    screen[i] = 0.5
        if gbmodel == 'GBn2':
            if elements is None:
                raise Exception('GBn2 model requires element information')
            for i, element in enumerate(elements):
                if element is elem.carbon:
                    screen[i] = 1.058554
                    alpha[i] = 0.733756
                    beta[i] = 0.506378
                    gamma[i] = 0.205844
                elif element is elem.hydrogen:
                    screen[i] = 1.425952
                    alpha[i] = 0.788440
                    beta[i] = 0.798699
                    gamma[i] = 0.437334
                elif element is elem.nitrogen:
                    screen[i] = 0.733599
                    alpha[i] = 0.503364
                    beta[i] = 0.316828
                    gamma[i] = 0.192915
                elif element is elem.oxygen:
                    screen[i] = 1.061039
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
                elif element is elem.sulfur:
                    screen[i] = -0.703469
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
                else: # not optimized
                    screen[i] = 0.5
                    alpha[i] = 1.0
                    beta[i] = 0.8
                    gamma[i] = 4.85
        lengthConversionFactor = units.angstrom.conversion_factor_to(units.nanometer)
        if gbmodel == 'GBn2':
            for rad, scr, alp, bet, gam in zip(radii, screen, alpha, beta, gamma):
                gb_List.append((rad*lengthConversionFactor, scr, alp, bet, gam))
        else:
            for rad, scr in zip(radii, screen):
                gb_List.append((rad*lengthConversionFactor, scr))
        return gb_List

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

    @property
    def has_scee_scnb(self):
        return ("SCEE_SCALE_FACTOR" in self._raw_data and "SCNB_SCALE_FACTOR" in self._raw_data)

    @property
    def has_atomic_number(self):
        return 'ATOMIC_NUMBER' in self._raw_data

#=============================================================================================
# AMBER System builder (based on, but not identical to, systemManager from 'zander')
#=============================================================================================

def readAmberSystem(prmtop_filename=None, prmtop_loader=None, shake=None, gbmodel=None,
          soluteDielectric=1.0, solventDielectric=78.5,
          implicitSolventKappa=0.0*(1/units.nanometer), nonbondedCutoff=None,
          nonbondedMethod='NoCutoff', scee=None, scnb=None, mm=None, verbose=False,
          EwaldErrorTolerance=None, flexibleConstraints=True, rigidWater=True, elements=None):
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
      implicitSolventKappa (float) - The Debye screening parameter corresponding to implicit solvent ionic strength
      nonbondedCutoff (float) - if specified, will set nonbondedCutoff (default: None)
      scnb (float) - 1-4 Lennard-Jones scaling factor (default: taken from prmtop or 1.2 if not present there)
      scee (float) - 1-4 electrostatics scaling factor (default: taken from prmtop or 2.0 if not present there)
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

    if prmtop.has_scee_scnb and (scee is not None or scnb is not None):
        warnings.warn("1-4 scaling parameters in topology file are being ignored. "
            "This is not recommended unless you know what you are doing.")

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
        for i in range(numConstrainedBonds):
            c = system.getConstraintParameters(i)
            distance = c[2].value_in_unit(units.nanometer)
            atomConstraints[c[0]].append((c[1], distance))
            atomConstraints[c[1]].append((c[0], distance))
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
            length = sqrt(l1*l1 + l2*l2 - 2*l1*l2*cos(aMin))
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
        boxBeta = boxBeta.value_in_unit(units.degrees)
        boxX = boxX.value_in_unit(units.angstroms)
        boxY = boxY.value_in_unit(units.angstroms)
        boxZ = boxZ.value_in_unit(units.angstroms)
        tmp = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]]
        _box_vectors_from_lengths_angles([boxX, boxY, boxZ],
                                         [boxBeta, boxBeta, boxBeta],
                                         tmp)
        xVec = units.Quantity(tmp[0], units.angstroms)
        yVec = units.Quantity(tmp[1], units.angstroms)
        zVec = units.Quantity(tmp[2], units.angstroms)
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
    _scee, _scnb = scee, scnb
    for (iAtom, lAtom, chargeProd, rMin, epsilon, iScee, iScnb) in prmtop.get14Interactions():
        if scee is None: _scee = iScee
        if scnb is None: _scnb = iScnb
        chargeProd /= _scee
        epsilon /= _scnb
        sigma = rMin * sigmaScale
        force.addException(iAtom, lAtom, chargeProd, sigma, epsilon)
        excludedAtomPairs.add(min((iAtom, lAtom), (lAtom, iAtom)))

    # Add Excluded Atoms
    excludedAtoms=prmtop.getExcludedAtoms()
    excludeParams = (0.0, 0.1, 0.0)
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
                    weightH = distOE[res]/sqrt(distOH[res]**2-(0.5*distHH[res])**2)
                    system.setVirtualSite(waterEP[res][0], mm.ThreeParticleAverageSite(waterO[res][0], waterH[res][0], waterH[res][1], 1-weightH, weightH/2, weightH/2))
                elif len(waterEP[res]) == 2:
                    # Five point water
                    weightH = cosOOP*distOE[res]/sqrt(distOH[res]**2-(0.5*distHH[res])**2)
                    angleHOH = 2*asin(0.5*distHH[res]/distOH[res])
                    lenCross = (distOH[res]**2)*sin(angleHOH)
                    weightCross = sinOOP*distOE[res]/lenCross
                    system.setVirtualSite(waterEP[res][0], mm.OutOfPlaneSite(waterO[res][0], waterH[res][0], waterH[res][1], weightH/2, weightH/2, weightCross))
                    system.setVirtualSite(waterEP[res][1], mm.OutOfPlaneSite(waterO[res][0], waterH[res][0], waterH[res][1], weightH/2, weightH/2, -weightCross))

    # Add GBSA model.
    if gbmodel is not None:
        # Convert implicitSolventKappa to nanometers if it is a unit.
        if units.is_quantity(implicitSolventKappa):
            implicitSolventKappa = implicitSolventKappa.value_in_unit((1/units.nanometers).unit)
        if verbose: print "Adding GB parameters..."
        charges = prmtop.getCharges()
        cutoff = None
        if nonbondedMethod != 'NoCutoff':
            cutoff = nonbondedCutoff
            if units.is_quantity(cutoff):
                cutoff = cutoff.value_in_unit(units.nanometers)
        gb_parms = prmtop.getGBParms(gbmodel, elements)
        if gbmodel == 'HCT':
            gb = customgb.GBSAHCTForce(solventDielectric, soluteDielectric, 'ACE', cutoff, implicitSolventKappa)
        elif gbmodel == 'OBC1':
            gb = customgb.GBSAOBC1Force(solventDielectric, soluteDielectric, 'ACE', cutoff, implicitSolventKappa)
        elif gbmodel == 'OBC2':
            if implicitSolventKappa > 0:
                gb = customgb.GBSAOBC2Force(solventDielectric, soluteDielectric, 'ACE', cutoff, implicitSolventKappa)
            else:
                gb = mm.GBSAOBCForce()
                gb.setSoluteDielectric(soluteDielectric)
                gb.setSolventDielectric(solventDielectric)
        elif gbmodel == 'GBn':
            gb = customgb.GBSAGBnForce(solventDielectric, soluteDielectric, 'ACE', cutoff, implicitSolventKappa)
        elif gbmodel == 'GBn2':
            gb = customgb.GBSAGBn2Force(solventDielectric, soluteDielectric, 'ACE', cutoff, implicitSolventKappa)
        else:
            raise Exception("Illegal value specified for implicit solvent model")
        for iAtom in range(prmtop.getNumAtoms()):
            if gbmodel == 'OBC2' and implicitSolventKappa == 0:
                gb.addParticle(charges[iAtom], gb_parms[iAtom][0], gb_parms[iAtom][1])
            elif gbmodel == 'GBn2':
                gb.addParticle([charges[iAtom], gb_parms[iAtom][0], gb_parms[iAtom][1],
                                gb_parms[iAtom][2], gb_parms[iAtom][3], gb_parms[iAtom][4]])
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
        # This applies the reaction field dielectric to the NonbondedForce
        # created above. Do not bind force to another name before this!
        force.setReactionFieldDielectric(1.0)

    # TODO: Add GBVI terms?

    return system

#=============================================================================================
# AMBER INPCRD loader classes
#=============================================================================================

class AmberAsciiRestart(object):
    """
    Class responsible for parsing Amber coordinates in the ASCII format.
    Automatically detects the presence of velocities or box parameters in the
    file.

    Parameters
    ----------
    filename : str
        Name of the restart file
    asNumpy : bool (False)
        Load the coordinates, velocities, and box as numpy ndarray objects

    Attributes
    ----------
    coordinates : natom x 3 array, Quantity
        Particle positions with units of length
    velocities : natom x 3 array, Quantity
        Particle velocities with units of length per time (None if velocities
        are not present in the inpcrd file)
    boxVectors : 3 x 3 array, Quantity
        Box vectors with units of length (None if no box is present in the
        inpcrd file)
    time : float, Quantity
        Simulation time (None if not present) with units of time
    title : str
        Title of the inpcrd file
    filename : str
        Name of the file we are parsing
    natom : int
        Number of atoms in the inpcrd file

    Raises
    ------
        `IOError' if the file does not exist
        `TypeError' if the format of the file is not recognized
        `ValueError' if not all fields are numbers (for example, if a field is
                     filled with ****'s)
        `IndexError' if the file is empty
        `ImportError' if numpy is requested but could not be imported
    Example
    -------
    >>> f = AmberAsciiRestart('alanine-dipeptide.inpcrd')
    >>> coordinates = f.coordinates
    """

    def __init__(self, filename, asNumpy=False):
        # Make sure numpy is available if requested
        if asNumpy and np is None:
            raise ImportError('asNumpy=True: numpy is not available')
        self._asNumpy = asNumpy
        self.filename = filename
        with open(filename, 'r') as f:
            lines = f.readlines()
            # Get rid of trailing blank lines
            while lines and not lines[-1].strip():
                lines.pop()
            self._parse(lines)

    def __str__(self):
        return self.filename

    def _parse(self, lines):
        """ Parses through the inpcrd file """
        global VELSCALE
        self.title = lines[0].strip()
        self.time = None

        try:
            words = lines[1].split()
            self.natom = int(words[0])
        except (IndexError, ValueError):
            raise TypeError('Unrecognized file type [%s]' % self.filename)
        
        if len(words) >= 2:
            self.time = float(words[1]) * units.picoseconds

        if len(lines) == int(ceil(self.natom / 2.0) + 2):
            hasbox = hasvels = False
            self.boxVectors = self.velocities = None
        elif self.natom in (1, 2) and len(lines) == 4:
            # This is the _only_ case where line counting does not work -- there
            # is either 1 or 2 atoms and there are 4 lines. The 1st 3 lines are
            # the title, natom/time, and coordinates. The 4th are almost always
            # velocities since Amber does not make it easy to make a periodic
            # system with only 2 atoms. If natom is 1, the 4th line is either a
            # velocity (3 #'s) or a box (6 #'s). If natom is 2, it is a bit
            # ambiguous. However, velocities (which are scaled by 20.445) have a
            # ~0% chance of being 60+, so we can pretty easily tell if the last
            # line has box dimensions and angles or velocities. I cannot
            # envision a _plausible_ scenario where the detection here will fail
            # in real life.
            line = lines[3]
            if self.natom == 1:
                tmp = [line[i:i+12] for i in range(0, 72, 12) if line[i:i+12]]
                if len(tmp) == 3:
                    hasvels = True
                    hasbox = False
                    self.boxVectors = False
                elif len(tmp) == 6:
                    hasbox = True
                    hasvels = False
                    self.velocities = None
                else:
                    raise TypeError('Unrecognized line in restart file %s' %
                                    self.filename)
            else:
                # Ambiguous case
                tmp = [float(line[i:i+12]) >= 60.0 for i in range(0, 72, 12)]
                if any(tmp):
                    hasbox = True
                    hasvels = False
                    self.velocities = False
                else:
                    hasvels = True
                    hasbox = False
                    self.boxVectors = False
        elif len(lines) == int(ceil(self.natom / 2.0) + 3):
            hasbox = True
            hasvels = False
            self.velocities = None
        elif len(lines) == int(2 * ceil(self.natom / 2.0) + 2):
            hasbox = False
            self.boxVectors = None
            hasvels = True
        elif len(lines) == int(2 * ceil(self.natom / 2.0) + 3):
            hasbox = hasvels = True
        else:
            raise TypeError('Badly formatted restart file. Has %d lines '
                            'for %d atoms.' % (len(self.lines), self.natom))

        if self._asNumpy:
            coordinates = np.zeros((self.natom, 3), np.float32)
            if hasvels:
                velocities = np.zeros((self.natom, 3), np.float32)
            if hasbox:
                boxVectors = np.zeros((3, 3), np.float32)
        else:
            coordinates = [[0.0, 0.0, 0.0] for i in range(self.natom)]
            if hasvels:
                velocities = [[0.0, 0.0, 0.0] for i in range(self.natom)]
            if hasbox:
                boxVectors = [[0.0, 0.0, 0.0] for i in range(3)]

        # Now it's time to parse.  Coordinates first
        startline = 2
        endline = startline + int(ceil(self.natom / 2.0))
        idx = 0
        for i in range(startline, endline):
            line = lines[i]
            coordinates[idx][0] = float(line[ 0:12])
            coordinates[idx][1] = float(line[12:24])
            coordinates[idx][2] = float(line[24:36])
            idx += 1
            if idx < self.natom:
                coordinates[idx][0] = float(line[36:48])
                coordinates[idx][1] = float(line[48:60])
                coordinates[idx][2] = float(line[60:72])
                idx += 1
        self.coordinates = units.Quantity(coordinates, units.angstroms)
        startline = endline
        # Now it's time to parse velocities if we have them
        if hasvels:
            endline = startline + int(ceil(self.natom / 2.0))
            idx = 0
            for i in range(startline, endline):
                line = lines[i]
                velocities[idx][0] = float(line[ 0:12]) * VELSCALE
                velocities[idx][1] = float(line[12:24]) * VELSCALE
                velocities[idx][2] = float(line[24:36]) * VELSCALE
                idx += 1
                if idx < self.natom:
                    velocities[idx][0] = float(line[36:48]) * VELSCALE
                    velocities[idx][1] = float(line[48:60]) * VELSCALE
                    velocities[idx][2] = float(line[60:72]) * VELSCALE
                    idx += 1
            startline = endline
            self.velocities = units.Quantity(velocities,
                                             units.angstroms/units.picoseconds)
        if hasbox:
            line = lines[startline]
            try:
                tmp = [float(line[i:i+12]) for i in range(0, 72, 12)]
            except (IndexError, ValueError):
                raise ValueError('Could not parse box line in %s' %
                                 self.filename)
            lengths = tmp[:3]
            angles = tmp[3:]
            _box_vectors_from_lengths_angles(lengths, angles, boxVectors)
            self.boxVectors = units.Quantity(boxVectors, units.angstroms)

class AmberNetcdfRestart(object):
    """
    Amber restart/inpcrd file in the NetCDF format (full double-precision
    coordinates, velocities, and unit cell parameters). Reads NetCDF restarts
    written by LEaP and pmemd/sander. Requires scipy to parse NetCDF files.

    Parameters
    ----------
    filename : str
        Name of the restart file
    asNumpy : bool (False)
        Load the coordinates, velocities, and box as numpy ndarray objects

    Attributes
    ----------
    coordinates : natom x 3 array, Quantity
        Particle positions with units of length
    velocities : natom x 3 array, Quantity
        Particle velocities with units of length per time (None if velocities
        are not present in the inpcrd file)
    boxVectors : 3 x 3 array, Quantity
        Box vectors with units of length (None if no box is present in the
        inpcrd file)
    time : float, Quantity
        Simulation time (None if not present) with units of time
    title : str
        Title of the inpcrd file
    filename : str
        Name of the file we are parsing
    natom : int
        Number of atoms in the inpcrd file

    Raises
    ------
        `IOError' if the file does not exist
        `TypeError' if the file is not a NetCDF v3 file
        `ImportError' if scipy is not available
    Example
    -------
    >>> f = AmberNetcdfRestart('alanine-dipeptide.ncrst')
    >>> coordinates = f.coordinates
    """
    def __init__(self, filename, asNumpy=False):
        try:
            from scipy.io.netcdf import NetCDFFile
        except ImportError:
            raise ImportError('scipy is necessary to parse NetCDF restarts')
        
        self.filename = filename
        self.velocities = self.boxVectors = self.time = None

        # Extract the information from the NetCDF file. We need to make copies
        # here because the NetCDF variables are mem-mapped, but is only mapped
        # to valid memory while the file handle is open. Since the context
        # manager GCs the ncfile handle, the memory for the original variables
        # is no longer valid. So copy those arrays while the handle is still
        # open. This is unnecessary in scipy v.0.12 and lower because NetCDFFile
        # accidentally leaks the file handle, but that was 'fixed' in 0.13. This
        # fix taken from MDTraj
        ncfile = NetCDFFile(filename, 'r')
        try:
            self.natom = ncfile.dimensions['atom']
            self.coordinates = np.array(ncfile.variables['coordinates'][:])
            if 'velocities' in ncfile.variables:
                vels = ncfile.variables['velocities']
                self.velocities = np.array(vels[:]) * vels.scale_factor
            if ('cell_lengths' in ncfile.variables and
                'cell_angles' in ncfile.variables):
                self.boxVectors = np.zeros((3,3), np.float32)
                _box_vectors_from_lengths_angles(
                        ncfile.variables['cell_lengths'][:],
                        ncfile.variables['cell_angles'][:],
                        self.boxVectors,
                )
            if 'time' in ncfile.variables:
                self.time = ncfile.variables['time'].getValue()
        finally:
            ncfile.close()

        # They are already numpy -- convert to list if we don't want numpy
        if not asNumpy:
            self.coordinates = self.coordinates.tolist()
            if self.velocities is not None:
                self.velocities = self.velocities.tolist()
            if self.boxVectors is not None:
                self.boxVectors = self.boxVectors.tolist()

        # Now add the units
        self.coordinates = units.Quantity(self.coordinates, units.angstroms)
        if self.velocities is not None:
            self.velocities = units.Quantity(self.velocities,
                                             units.angstroms/units.picoseconds)
        if self.boxVectors is not None:
            self.boxVectors = units.Quantity(self.boxVectors, units.angstroms)
        self.time = units.Quantity(self.time, units.picosecond)

def _box_vectors_from_lengths_angles(lengths, angles, boxVectors):
    """
    Converts lengths and angles into a series of box vectors and modifies
    boxVectors in-place (it must be a mutable sequence)

    Parameters
    ----------
    lengths : 3-element array of floats
        Lengths of the 3 periodic box vectors
    angles : 3-element array of floats
        Angles (in degrees) between the 3 periodic box vectors
    boxVectors : mutable 3x3 sequence
    """
    alpha = angles[0] * pi / 180.0
    beta = angles[1] * pi / 180.0
    gamma = angles[2] * pi / 180.0

    boxVectors[0][0] = lengths[0]

    boxVectors[1][0] = lengths[1] * cos(gamma)
    boxVectors[1][1] = lengths[1] * sin(gamma)

    boxVectors[2][0] = cx = lengths[2] * cos(beta)
    boxVectors[2][1] = cy = lengths[2] * (cos(alpha) - cos(beta) * cos(gamma))
    boxVectors[2][2] = sqrt(lengths[2]*lengths[2] - cx*cx - cy*cy)

    boxVectors[0][1] = boxVectors[0][2] = boxVectors[1][2] = 0.0

    # Now make sure any vector close to zero is zero exactly
    for i in range(3):
        for j in range(3):
            if abs(boxVectors[i][j]) < TINY:
                boxVectors[i][j] = 0.0

def readAmberCoordinates(filename, asNumpy=False):
    """
    Read atomic coordinates (and optionally, box vectors) from Amber formatted coordinate file.

    ARGUMENTS

    filename (string) - name of Amber coordinates file to be read in

    OPTIONAL ARGUMENTS

    asNumpy (boolean) - if True, results will be returned as Numpy arrays instead of lists of Vec3s

    RETURNS

    coordinates, velocities, boxVectors
        The velocities and boxVectors will be None if they are not found in the
        restart file

    EXAMPLES

    Read coordinates in vacuum.

    >>> directory = os.path.join(os.getenv('YANK_INSTALL_DIR'), 'test', 'systems', 'alanine-dipeptide-gbsa')
    >>> crd_filename = os.path.join(directory, 'alanine-dipeptide.inpcrd')
    >>> coordinates, velocities, box_vectors = readAmberCoordinates(crd_filename)

    Read coordinates in solvent.

    >>> directory = os.path.join(os.getenv('YANK_INSTALL_DIR'), 'test', 'systems', 'alanine-dipeptide-explicit')
    >>> crd_filename = os.path.join(directory, 'alanine-dipeptide.inpcrd')
    >>> coordinates, velocities, box_vectors = readAmberCoordinates(crd_filename)
    """

    try:
        crdfile = AmberNetcdfRestart(filename)
    except ImportError:
        # See if it's an ASCII file.  If so, no need to complain
        try:
            crdfile = AmberAsciiRestart(filename)
        except TypeError:
            raise TypeError('Problem parsing %s as an ASCII Amber restart file '
                            'and scipy could not be imported to try reading as '
                            'a NetCDF restart file.' % filename)
        except (IndexError, ValueError):
            raise TypeError('Could not parse Amber ASCII restart file %s' %
                            filename)
        except ImportError:
            raise ImportError('Could not find numpy; cannot use asNumpy=True')
    except TypeError:
        # We had scipy, but this is not a NetCDF v3 file. Try as ASCII now
        try:
            crdfile = AmberAsciiRestart(filename)
        except TypeError:
            raise TypeError('Problem parsing %s as an ASCII Amber restart file'
                            % filename)
        except (IndexError, ValueError):
            raise TypeError('Could not parse Amber ASCII restart file %s' %
                            filename)
        # Import error cannot happen, since we had scipy which has numpy as a
        # prereq. Do not catch that exception (only catch what you intend to
        # catch...)

    # We got here... one of the file types worked. Return the coordinates,
    # velocities, and boxVectors
    return crdfile.coordinates, crdfile.velocities, crdfile.boxVectors

#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

if __name__ == "__main__":
    import doctest
    doctest.testmod()
