"""
Provides a Python class for parsing a PSF file and setting up a system
structure for it within the OpenMM framework.

This file is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of Biological
Structures at Stanford, funded under the NIH Roadmap for Medical Research,
grant U54 GM072970. See https://simtk.org.  This code was originally part of
the ParmEd program and was ported for use with OpenMM.

Copyright (c) 2014 the Authors

Author: Jason M. Swails
Contributors:
Date: April 18, 2014

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
from __future__ import division

from functools import wraps
from math import pi, cos, sin, sqrt
import os
import simtk.openmm as mm
from simtk.openmm.vec3 import Vec3
import simtk.unit as u
from simtk.openmm.app import (forcefield as ff, Topology, element)
from simtk.openmm.app.amberprmtopfile import HCT, OBC1, OBC2, GBn, GBn2
from simtk.openmm.app.internal.customgbforces import (GBSAHCTForce,
                GBSAOBC1Force, GBSAOBC2Force, GBSAGBnForce, GBSAGBn2Force)
# CHARMM imports
from simtk.openmm.app.internal.charmm._charmmfile import CharmmFile
from simtk.openmm.app.internal.charmm.topologyobjects import (
                ResidueList, AtomList, TrackedList, Bond, Angle, Dihedral,
                Improper, AcceptorDonor, Group, Cmap, UreyBradley,
                NoUreyBradley)
from simtk.openmm.app.internal.charmm.exceptions import (
                CharmmPSFError, MoleculeError, CharmmPSFWarning,
                MissingParameter)
import warnings

TINY = 1e-8
WATNAMES = ('WAT', 'HOH', 'TIP3', 'TIP4', 'TIP5', 'SPCE', 'SPC')


def _catchindexerror(func):
    """
    Protects a function from raising an index error, and replace that exception
    with a CharmmPSFError instead
    """
    @wraps(func)
    def newfunc(*args, **kwargs):
        """ Catch the index error """
        try:
            return func(*args, **kwargs)
        except IndexError, e:
            raise CharmmPSFError('Array is too short: %s' % e)

    return newfunc

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CharmmPsfFile(object):
    """
    A chemical structure instantiated from CHARMM files.

    Example:
    >>> cs = CharmmPsfFile("testfiles/test.psf")
    
    This structure has numerous attributes that are lists of the elements of
    this structure, including atoms, bonds, torsions, etc. The attributes are
        - residue_list
        - atom_list
        - bond_list
        - angle_list
        - dihedral_list
        - improper_list
        - cmap_list
        - donor_list    # hbonds donors?
        - acceptor_list # hbond acceptors?
        - group_list    # list of nonbonded interaction groups

    Additional attribute is available if a CharmmParameterSet is loaded into
    this structure.
        
        - urey_bradley_list

    The lengths of each of these lists gives the pointers (e.g., natom, nres,
    etc.)

    Example:
    >>> cs = CharmmPsfFile("testfiles/test.psf")
    >>> len(cs.atom_list)
    33
    >>> len(cs.bond_list)
    32
    """
    # Define default force groups for all of the bonded terms. This allows them
    # to be turned on and off selectively. This is a way to implement per-term
    # energy decomposition to compare individual components

    BOND_FORCE_GROUP = 0
    ANGLE_FORCE_GROUP = 1
    DIHEDRAL_FORCE_GROUP = 2
    UREY_BRADLEY_FORCE_GROUP = 3
    IMPROPER_FORCE_GROUP = 4
    CMAP_FORCE_GROUP = 5
    NONBONDED_FORCE_GROUP = 6
    GB_FORCE_GROUP = 6
    
    @_catchindexerror
    def __init__(self, psf_name):
        """
        Opens and parses a PSF file, then instantiates a CharmmPsfFile instance
        from the data.
            
        Parameters:
            psf_name (str) : Name of the PSF file (it must exist)
        
        Exceptions Raised:
            IOError : If file "psf_name" does not exist
            CharmmPSFError: If any parsing errors are encountered
        """
        conv = CharmmPsfFile._convert
        # Make sure the file exists
        if not os.path.exists(psf_name):
            raise IOError('Could not find PSF file %s' % psf_name)
        # Open the PSF and read the first line. It must start with "PSF"
        psf = CharmmFile(psf_name, 'r')
        line = psf.readline()
        if not line.startswith('PSF'):
            raise CharmmPSFError('Unrecognized PSF file. First line is %s' %
                                 line.strip())
        # Store the flags
        psf_flags = line.split()[1:]
        # Next line is blank
        psf.readline()
        # The next line has one number -- the number of title lines
        ntitle = conv(psf.readline().strip(), int, 'title count')
        # store the title
        title = list()
        for i in range(ntitle):
            title.append(psf.readline().rstrip())
        # Skip the blank line
        psf.readline()
        # Next is the number of atoms
        natom = conv(psf.readline().strip(), int, 'natom')
        # Parse all of the atoms
        residue_list = ResidueList()
        atom_list = AtomList()
        for i in xrange(natom):
            words = psf.readline().split()
            atid = int(words[0])
            if atid != i + 1:
                raise CharmmPSFError('Nonsequential atoms detected!')
            system = words[1]
            resid = conv(words[2], int, 'residue number')
            resname = words[3]
            name = words[4]
            attype = words[5]
            # Try to convert the atom type to an integer a la CHARMM
            try:
                attype = int(attype)
            except ValueError:
                pass
            charge = conv(words[6], float, 'partial charge')
            mass = conv(words[7], float, 'atomic mass')
            props = words[8:]
            atom = residue_list.add_atom(system, resid, resname, name,
                                         attype, charge, mass, props)
            atom_list.append(atom)
        atom_list.assign_indexes()
        # Eat the next line
        psf.readline()
        # Now get the number of bonds
        nbond, holder = CharmmPsfFile._parse_psf_section(psf, int)
        bond_list = TrackedList()
        if len(holder) != nbond * 2:
            raise CharmmPSFError('Got %d indexes for %d bonds' %
                                 (len(holder), nbond))
        for i in range(nbond):
            id1 = holder[2*i  ] - 1
            id2 = holder[2*i+1] - 1
            bond_list.append(Bond(atom_list[id1], atom_list[id2]))
        bond_list.changed = False
        # Now get the number of angles and the angle list
        ntheta, holder = CharmmPsfFile._parse_psf_section(psf, int)
        angle_list = TrackedList()
        if len(holder) != ntheta * 3:
            raise CharmmPSFError('Got %d indexes for %d angles' %
                                 (len(holder), ntheta))
        for i in range(ntheta):
            id1 = holder[3*i  ] - 1
            id2 = holder[3*i+1] - 1
            id3 = holder[3*i+2] - 1
            angle_list.append(
                    Angle(atom_list[id1], atom_list[id2], atom_list[id3])
            )
        angle_list.changed = False
        # Now get the number of torsions and the torsion list
        nphi, holder = CharmmPsfFile._parse_psf_section(psf, int)
        dihedral_list = TrackedList()
        if len(holder) != nphi * 4:
            raise CharmmPSFError('Got %d indexes for %d torsions' %
                                 (len(holder), nphi))
        for i in range(nphi):
            id1 = holder[4*i  ] - 1
            id2 = holder[4*i+1] - 1
            id3 = holder[4*i+2] - 1
            id4 = holder[4*i+3] - 1
            dihedral_list.append(
                    Dihedral(atom_list[id1], atom_list[id2], atom_list[id3],
                             atom_list[id4])
            )
        dihedral_list.changed = False
        # Now get the number of improper torsions
        nimphi, holder = CharmmPsfFile._parse_psf_section(psf, int)
        improper_list = TrackedList()
        if len(holder) != nimphi * 4:
            raise CharmmPSFError('Got %d indexes for %d impropers' %
                                 (len(holder), nimphi))
        for i in range(nimphi):
            id1 = holder[4*i  ] - 1
            id2 = holder[4*i+1] - 1
            id3 = holder[4*i+2] - 1
            id4 = holder[4*i+3] - 1
            improper_list.append(
                    Improper(atom_list[id1], atom_list[id2], atom_list[id3],
                             atom_list[id4])
            )
        improper_list.changed = False
        # Now handle the donors (what is this used for??)
        ndon, holder = CharmmPsfFile._parse_psf_section(psf, int)
        donor_list = TrackedList()
        if len(holder) != ndon * 2:
            raise CharmmPSFError('Got %d indexes for %d donors' %
                                 (len(holder), ndon))
        for i in range(ndon):
            id1 = holder[2*i  ] - 1
            id2 = holder[2*i+1] - 1
            donor_list.append(AcceptorDonor(atom_list[id1], atom_list[id2]))
        donor_list.changed = False
        # Now handle the acceptors (what is this used for??)
        nacc, holder = CharmmPsfFile._parse_psf_section(psf, int)
        acceptor_list = TrackedList()
        if len(holder) != nacc * 2:
            raise CharmmPSFError('Got %d indexes for %d acceptors' %
                                 (len(holder), ndon))
        for i in range(nacc):
            id1 = holder[2*i  ] - 1
            id2 = holder[2*i+1] - 1
            acceptor_list.append(AcceptorDonor(atom_list[id1], atom_list[id2]))
        acceptor_list.changed = False
        # Now get the NNB section. Not sure what this section is for or what it
        # does...
        nnb, holder = CharmmPsfFile._parse_psf_section(psf, int)
        # Now get the group sections
        pointers, holder = CharmmPsfFile._parse_psf_section(psf, int)
        group_list = TrackedList()
        try:
            ngrp, nst2 = pointers
        except ValueError:
            raise CharmmPSFError('Could not unpack GROUP pointers')
        group_list.nst2 = nst2
        # Now handle the groups
        if len(holder) != ngrp * 3:
            raise CharmmPSFError('Got %d indexes for %d groups' %
                                 (len(holder), ngrp))
        for i in range(ngrp):
            i1 = holder[3*i  ]
            i2 = holder[3*i+1]
            i3 = holder[3*i+2]
            group_list.append(Group(i1, i2, i3))
        group_list.changed = False
        # The next section might be the number of molecules or it might be the
        # cross-term (cmap) section. The first thing we'll do is determine
        # molecularity based on the atom connectivity. If every PSF file was
        # guaranteed to be "correct", we could just compare the MOLNT
        # section with the one we compute here. However, CHARMM GUI appears
        # to assign MOLNT as a dummy section (with all 1's), so this
        # approach will not work. Instead, look at the value of the pointer
        # and the number of entries in the group. If the # of entries is
        # NATOM, assume we have MOLNT section. Warn if the MOLNT section is
        # 'wrong'...
        pointer, holder = CharmmPsfFile._parse_psf_section(psf, int)

        # Assign all of the atoms to molecules recursively
        set_molecules(atom_list)
        molecule_list = [atom.marked for atom in atom_list]
        if len(holder) == len(atom_list):
            if molecule_list != holder:
                warnings.warn('Detected PSF molecule section that is WRONG. '
                              'Resetting molecularity.', CharmmPSFWarning)
            # We have a CHARMM PSF file; now do NUMLP/NUMLPH sections
            words = psf.readline().split()
            numlp = conv(words[0], int, 'numlp')
            numlph = conv(words[1], int, 'numlph')
            if numlp != 0 or numlph != 0:
                raise NotImplemented('Cannot currently handle PSFs with lone '
                                     'pairs defined in the NUMLP/NUMLPH '
                                     'section.')
            psf.readline() # blank
            # Now we get to the cross-term section
            ncrterm, holder = CharmmPsfFile._parse_psf_section(psf, int)
        else:
            ncrterm = pointer
        # At this point, ncrterm and holder are both set to the CMAP list for
        # VMD and non-VMD PSFs.
        # Now get the cmaps
        cmap_list = TrackedList()
        if len(holder) != ncrterm * 8:
            raise CharmmPSFError('Got %d CMAP indexes for %d cmap terms' %
                                 (len(holder), ncrterm))
        for i in range(ncrterm):
            id1 = holder[8*i  ] - 1
            id2 = holder[8*i+1] - 1
            id3 = holder[8*i+2] - 1
            id4 = holder[8*i+3] - 1
            id5 = holder[8*i+4] - 1
            id6 = holder[8*i+5] - 1
            id7 = holder[8*i+6] - 1
            id8 = holder[8*i+7] - 1
            cmap_list.append(
                    Cmap(atom_list[id1], atom_list[id2], atom_list[id3],
                         atom_list[id4], atom_list[id5], atom_list[id6],
                         atom_list[id7], atom_list[id8])
            )
        cmap_list.changed = False

        # Now set the instance attributes
        self.residue_list = residue_list
        self.atom_list = atom_list
        self.bond_list = bond_list
        self.angle_list = angle_list
        self.dihedral_list = dihedral_list
        self.improper_list = improper_list
        self.cmap_list = cmap_list
        self.donor_list = donor_list
        self.acceptor_list = acceptor_list
        self.group_list = group_list
        self.title = title
        self.flags = psf_flags
        self.box_vectors = None

    @staticmethod
    def _convert(string, type, message):
        """
        Converts a string to a specific type, making sure to raise
        CharmmPSFError with the given message in the event of a failure.

        Parameters:
            - string (str) : Input string to process
            - type (type) : Type of data to convert to
            - message (str) : Error message to put in exception if failed
        """
        try:
            return type(string)
        except ValueError, e:
            print e
            raise CharmmPSFError('Could not convert %s' % message)

    @staticmethod
    def _parse_psf_section(psf, dtype):
        """
        This method parses a section of the PSF file

        Parameters:
            - psf (CharmmFile) : Open file that is pointing to the first line
                                 of the section that is to be parsed
            - dtype (type) : The data type to convert all of the data into
        
        Returns:
            (pointers, data)
            
            - pointers (int/tuple of ints) : If one pointer is set, pointers is
                    simply the integer that is value of that pointer. Otherwise
                    it is a tuple with every pointer value defined in the first
                    line
            - data (list) : A list of all data in the parsed section converted
                    to `dtype'
        """
        conv = CharmmPsfFile._convert
        words = psf.readline().split()
        if len(words) == 1:
            pointers = conv(words[0], int, 'pointer')
        else:
            pointers = tuple([conv(w, int, 'pointer') for w in words])
        line = psf.readline().strip()
        if not line:
            # This will correctly handle the NNB section (which has a spurious
            # blank line) as well as any sections that have 0 members.
            line = psf.readline().strip()
        data = []
        while line:
            words = line.split()
            data.extend([conv(w, dtype, 'PSF data') for w in words])
            line = psf.readline().strip()
        return pointers, data

    def writePsf(self, dest, vmd=False):
        """
        Writes a PSF file from the stored molecule

        Parameters:
            - dest (str or file-like) : The place to write the output PSF file.
                    If it has a "write" attribute, it will be used to print the
                    PSF file. Otherwise, it will be treated like a string and a
                    file will be opened, printed, then closed
            - vmd (bool) : If True, it will write out a PSF in the format that
                    VMD prints it in (i.e., no NUMLP/NUMLPH or MOLNT sections)
        Example:
            >>> cs = CharmmPsfFile('testfiles/test.psf')
            >>> cs.writePsf('testfiles/test2.psf')
        """
        # See if this is an extended format
        ext = 'EXT' in self.flags
        own_handle = False
        # Index the atoms and residues
        self.residue_list.assign_indexes()
        self.atom_list.assign_indexes()
        if not hasattr(dest, 'write'):
            own_handle = True
            dest = open(dest, 'w')

        # Assign the formats we need to write with
        if ext:
            atmfmt1 = ('%10d %-8s %-8i %-8s %-8s %4d %10.6f %13.4f' + 11*' ')
            atmfmt2 = ('%10d %-8s %-8i %-8s %-8s %-4s %10.6f %13.4f' + 11*' ')
            intfmt = '%10d' # For pointers
        else:
            atmfmt1 = ('%8d %-4s %-4i %-4s %-4s %4d %10.6f %13.4f' + 11*' ')
            atmfmt2 = ('%8d %-4s %-4i %-4s %-4s %-4s %10.6f %13.4f' + 11*' ')
            intfmt = '%8d' # For pointers

        # Now print the header then the title
        dest.write('PSF ' + ' '.join(self.flags) + '\n')
        dest.write('\n')
        dest.write(intfmt % len(self.title) + ' !NTITLE\n')
        dest.write('\n'.join(self.title) + '\n\n')
        # Now time for the atoms
        dest.write(intfmt % len(self.atom_list) + ' !NATOM\n')
        # atmfmt1 is for CHARMM format (i.e., atom types are integers)
        # atmfmt is for XPLOR format (i.e., atom types are strings)
        for i, atom in enumerate(self.atom_list):
            if isinstance(atom.attype, str):
                fmt = atmfmt2
            else:
                fmt = atmfmt1
            atmstr = fmt % (i+1, atom.system, atom.residue.resnum,
                            atom.residue.resname, atom.name, atom.attype,
                            atom.charge, atom.mass)
            dest.write(atmstr + '   '.join(atom.props) + '\n')
        dest.write('\n')
        # Bonds
        dest.write(intfmt % len(self.bond_list) + ' !NBOND: bonds\n')
        for i, bond in enumerate(self.bond_list):
            dest.write((intfmt*2) % (bond.atom1.idx+1, bond.atom2.idx+1))
            if i % 4 == 3: # Write 4 bonds per line
                dest.write('\n')
        # See if we need to terminate
        if len(self.bond_list) % 4 != 0 or len(self.bond_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Angles
        dest.write(intfmt % len(self.angle_list) + ' !NTHETA: angles\n')
        for i, angle in enumerate(self.angle_list):
            dest.write((intfmt*3) % (angle.atom1.idx+1, angle.atom2.idx+1,
                                     angle.atom3.idx+1)
            )
            if i % 3 == 2: # Write 3 angles per line
                dest.write('\n')
        # See if we need to terminate
        if len(self.angle_list) % 3 != 0 or len(self.angle_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Dihedrals
        dest.write(intfmt % len(self.dihedral_list) + ' !NPHI: dihedrals\n')
        for i, dih in enumerate(self.dihedral_list):
            dest.write((intfmt*4) % (dih.atom1.idx+1, dih.atom2.idx+1,
                                     dih.atom3.idx+1, dih.atom4.idx+1)
            )
            if i % 2 == 1: # Write 2 dihedrals per line
                dest.write('\n')
        # See if we need to terminate
        if len(self.dihedral_list) % 2 != 0 or len(self.dihedral_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Impropers
        dest.write(intfmt % len(self.improper_list) + ' !NIMPHI: impropers\n')
        for i, imp in enumerate(self.improper_list):
            dest.write((intfmt*4) % (imp.atom1.idx+1, imp.atom2.idx+1,
                                     imp.atom3.idx+1, imp.atom4.idx+1)
            )
            if i % 2 == 1: # Write 2 dihedrals per line
                dest.write('\n')
        # See if we need to terminate
        if len(self.improper_list) % 2 != 0 or len(self.improper_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Donor section
        dest.write(intfmt % len(self.donor_list) + ' !NDON: donors\n')
        for i, don in enumerate(self.donor_list):
            dest.write((intfmt*2) % (don.atom1.idx+1, don.atom2.idx+1))
            if i % 4 == 3: # 4 donors per line
                dest.write('\n')
        if len(self.donor_list) % 4 != 0 or len(self.donor_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # Acceptor section
        dest.write(intfmt % len(self.acceptor_list) + ' !NACC: acceptors\n')
        for i, acc in enumerate(self.acceptor_list):
            dest.write((intfmt*2) % (acc.atom1.idx+1, acc.atom2.idx+1))
            if i % 4 == 3: # 4 donors per line
                dest.write('\n')
        if len(self.acceptor_list) % 4 != 0 or len(self.acceptor_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # NNB section ??
        dest.write(intfmt % 0 + ' !NNB\n\n')
        for i in range(len(self.atom_list)):
            dest.write(intfmt % 0)
            if i % 8 == 7: # Write 8 0's per line
                dest.write('\n')
        if len(self.atom_list) % 8 != 0: dest.write('\n')
        dest.write('\n')
        # Group section
        dest.write((intfmt*2) % (len(self.group_list), self.group_list.nst2))
        dest.write(' !NGRP NST2\n')
        for i, gp in enumerate(self.group_list):
            dest.write((intfmt*3) % (gp.bs, gp.type, gp.move))
            if i % 3 == 2: dest.write('\n')
        if len(self.group_list) % 3 != 0 or len(self.group_list) == 0:
            dest.write('\n')
        dest.write('\n')
        # The next two sections are never found in VMD prmtops...
        if not vmd:
            # Molecule section; first set molecularity
            set_molecules(self.atom_list)
            mollist = [a.marked for a in self.atom_list]
            dest.write(intfmt % max(mollist) + ' !MOLNT\n')
            for i, atom in enumerate(self.atom_list):
                dest.write(intfmt % atom.marked)
                if i % 8 == 7: dest.write('\n')
            if len(self.atom_list) % 8 != 0: dest.write('\n')
            dest.write('\n')
            # NUMLP/NUMLPH section
            dest.write((intfmt*2) % (0, 0) + ' !NUMLP NUMLPH\n')
            dest.write('\n')
        # CMAP section
        dest.write(intfmt % len(self.cmap_list) + ' !NCRTERM: cross-terms\n')
        for i, cmap in enumerate(self.cmap_list):
            dest.write((intfmt*4) % (cmap.atom1.idx+1, cmap.atom2.idx+1,
                                     cmap.atom3.idx+1, cmap.atom4.idx+1)
            )
            if cmap.consecutive:
                dest.write((intfmt*4) % (cmap.atom2.idx+1, cmap.atom3.idx+1,
                                         cmap.atom4.idx+1, cmap.atom5.idx+1)
                )
            else:
                dest.write((intfmt*4) % (cmap.atom5.idx+1, cmap.atom6.idx+1,
                                         cmap.atom7.idx+1, cmap.atom8.idx+1)
                )
            dest.write('\n')
        # Done!
        # If we opened our own handle, close it
        if own_handle:
            dest.close()
        
    def loadParameters(self, parmset):
        """
        Loads parameters from a parameter set that was loaded via CHARMM RTF,
        PAR, and STR files.

        Parameters:
            - parmset (CharmmParameterSet) : List of all parameters

        Notes:
            - If any parameters that are necessary cannot be found, a
              MissingParameter exception is raised.

            - If any dihedral or improper parameters cannot be found, I will try
              inserting wildcards (at either end for dihedrals and as the two
              central atoms in impropers) and see if that matches.  Wild-cards
              will apply ONLY if specific parameters cannot be found.

            - This method will expand the dihedral_list attribute by adding a
              separate Dihedral object for each term for types that have a
              multi-term expansion
        """
        # First load the atom types
        types_are_int = False
        for atom in self.atom_list:
            try:
                if isinstance(atom.attype, int):
                    atype = parmset.atom_types_int[atom.attype]
                    types_are_int = True # if we have to change back
                else:
                    atype = parmset.atom_types_str[atom.attype]
            except KeyError:
                raise MissingParameter('Could not find atom type for %s' %
                                       atom.attype)
            atom.type = atype
            # Change to string attype to look up the rest of the parameters
            atom.type_to_str()

        # Next load all of the bonds
        for bond in self.bond_list:
            # Construct the key
            key = (min(bond.atom1.attype, bond.atom2.attype),
                   max(bond.atom1.attype, bond.atom2.attype))
            try:
                bond.bond_type = parmset.bond_types[key]
            except KeyError:
                raise MissingParameter('Missing bond type for %r' % bond)
        # Next load all of the angles. If a Urey-Bradley term is defined for
        # this angle, also build the urey_bradley and urey_bradley_type lists
        self.urey_bradley_list = TrackedList()
        for ang in self.angle_list:
            # Construct the key
            key = (min(ang.atom1.attype, ang.atom3.attype), ang.atom2.attype,
                   max(ang.atom1.attype, ang.atom3.attype))
            try:
                ang.angle_type = parmset.angle_types[key]
                ubt = parmset.urey_bradley_types[key]
                if ubt is not NoUreyBradley:
                    ub = UreyBradley(ang.atom1, ang.atom3, ubt)
                    self.urey_bradley_list.append(ub)
            except KeyError:
                raise MissingParameter('Missing angle type for %r' % ang)
        # Next load all of the dihedrals. This is a little trickier since we
        # need to back up the existing dihedral list and replace it with a
        # longer one that has only one Fourier term per Dihedral instance.
        dihedral_list = self.dihedral_list
        self.dihedral_list = TrackedList()
        for dih in dihedral_list:
            # Store the atoms
            a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
            at1, at2, at3, at4 = a1.attype, a2.attype, a3.attype, a4.attype
            # First see if the exact dihedral is specified
            key = min((at1,at2,at3,at4), (at4,at3,at2,at1))
            if not key in parmset.dihedral_types:
                # Check for wild-cards
                key = min(('X',at2,at3,'X'), ('X',at3,at2,'X'))
                if not key in parmset.dihedral_types:
                    raise MissingParameter('No dihedral parameters found for '
                                           '%r' % dih)
            dtlist = parmset.dihedral_types[key]
            for i, dt in enumerate(dtlist):
                self.dihedral_list.append(Dihedral(a1, a2, a3, a4, dt))
                # See if we include the end-group interactions for this
                # dihedral. We do IFF it is the last or only dihedral term and
                # it is NOT in the angle/bond partners
                if i != len(dtlist) - 1:
                    self.dihedral_list[-1].end_groups_active = False
                elif a1 in a4.bond_partners or a1 in a4.angle_partners:
                    self.dihedral_list[-1].end_groups_active = False
        # Now do the impropers
        for imp in self.improper_list:
            # Store the atoms
            a1, a2, a3, a4 = imp.atom1, imp.atom2, imp.atom3, imp.atom4
            at1, at2, at3, at4 = a1.attype, a2.attype, a3.attype, a4.attype
            key = tuple(sorted([at1, at2, at3, at4]))
            if not key in parmset.improper_types:
                # Check for wild-cards
                for anchor in (at2, at3, at4):
                    key = tuple(sorted([at1, anchor, 'X', 'X']))
                    if key in parmset.improper_types:
                        break # This is the right key
            try:
                imp.improper_type = parmset.improper_types[key]
            except KeyError:
                raise MissingParameter('No improper parameters found for %r' %
                                       imp)
        # Now do the cmaps. These will not have wild-cards
        for cmap in self.cmap_list:
            # Store the atoms for easy reference
            if cmap.consecutive:
                a1, a2, a3, a4 = cmap.atom1, cmap.atom2, cmap.atom3, cmap.atom4
                a5, a6, a7, a8 = cmap.atom2, cmap.atom3, cmap.atom4, cmap.atom5
            else:
                a1, a2, a3, a4 = cmap.atom1, cmap.atom2, cmap.atom3, cmap.atom4
                a5, a6, a7, a8 = cmap.atom5, cmap.atom6, cmap.atom7, cmap.atom8
            at1, at2, at3, at4 = a1.attype, a2.attype, a3.attype, a4.attype
            at5, at6, at7, at8 = a5.attype, a6.attype, a7.attype, a8.attype
            # Construct the keys
            k1 = list(min((at1,at2,at3,at4), (at4,at3,at2,at1)))
            k2 = list(min((at5,at6,at7,at8), (at8,at7,at6,at5)))
            key = tuple(k1 + k2)
            try:
                cmap.cmap_type = parmset.cmap_types[key]
            except KeyError:
                raise MissingParameter('No CMAP parameters found for %r' % cmap)
        # If the types started out as integers, change them back
        if types_are_int:
            for atom in self.atom_list: atom.type_to_int()

    def setCoordinates(self, positions, velocities=None):
        """
        This method loads the coordinates and velocity information from an
        external object or passed data.

        Parameters:
            - positions (list of floats) : A 3-N length iterable with all of the
                coordinates in the order [x1, y1, z1, x2, y2, z2, ...].
            - velocities (list of floats) : If not None, is the velocity
                equivalent of the positions
        """
        if len(positions) / 3 != len(self.atom_list):
            raise ValueError('Coordinates given for %s atoms, but %d atoms '
                             'exist in this structure.' %
                             (len(positions)/3, len(self.atom_list)))
        # Now assign all of the atoms positions
        for i, atom in enumerate(self.atom_list):
            atom.xx = positions[3*i  ]
            atom.xy = positions[3*i+1]
            atom.xz = positions[3*i+2]

        # Do velocities if given
        if velocities is not None:
            if len(velocities) / 3 != len(self.atom_list):
                raise ValueError('Velocities given for %s atoms, but %d atoms '
                                 'exist in this structure.' %
                                 (len(velocities)/3, len(self.atom_list)))
            for i, atom in enumerate(self.atom_list):
                atom.vx = velocities[3*i  ]
                atom.vy = velocities[3*i+1]
                atom.vz = velocities[3*i+2]
            self.velocities = velocities

        self.positions = positions

    def setBox(self, a, b, c, alpha=90.0*u.degrees, beta=90.0*u.degrees,
               gamma=90.0*u.degrees):
        """
        Sets the periodic box boundary conditions.

        Parameters:
            - a, b, c (floats) : Lengths of the periodic cell
            - alpha, beta, gamma (floats, optional) : Angles between the
                periodic cells.
        """
        self.box_vectors = _box_vectors_from_lengths_angles(a, b, c,
                                                            alpha, beta, gamma)
        # If we already have a _topology instance, then we have possibly changed
        # the existence of box information (whether or not this is a periodic
        # system), so delete any cached reference to a topology so it's
        # regenerated with updated information
        if hasattr(self, '_topology'):
            del self._topology

    @property
    def topology(self):
        """ Create an OpenMM Topology object from the stored bonded network """
        try:
            return self._topology
        except AttributeError:
            # If none exists, we need to create it
            pass
        # Cache the topology for easy returning later
        self._topology = topology = Topology()
        
        last_chain = None
        last_residue = None
        # Add each chain (separate 'system's) and residue
        for atom in self.atom_list:
            if atom.system != last_chain:
                chain = topology.addChain()
                last_chain = atom.system
            if atom.residue.idx != last_residue:
                last_residue = atom.residue.idx
                residue = topology.addResidue(atom.residue.resname, chain)
            if atom.type is not None:
                # This is the most reliable way of determining the element
                atomic_num = atom.type.atomic_number
                elem = element.Element.getByAtomicNumber(atomic_num)
            else:
                # Figure it out from the mass
                elem = element.Element.getByMass(atom.mass)
            topology.addAtom(atom.name, elem, residue)

        # Add all of the bonds
        atoms = list(topology.atoms())
        # Assign atom indexes to make sure they're current
        self.atom_list.assign_indexes()
        for bond in self.bond_list:
            topology.addBond(atoms[bond.atom1.idx], atoms[bond.atom2.idx])

        # Add the periodic box if there is one
        if self.box_vectors is not None:
            topology.setUnitCellDimensions(self.boxLengths)

        return topology

    def _get_gb_params(self, gb_model=HCT):
        """ Gets the GB parameters. Need this method to special-case GB neck """
        screen = [0 for atom in self.atom_list]
        if gb_model is GBn:
            radii = _bondi_radii(self.atom_list)
            screen = [0.5 for atom in self.atom_list]
            for i, atom in enumerate(self.atom_list):
                if atom.type.atomic_number == 6:
                    screen[i] = 0.48435382330
                elif atom.type.atomic_number == 1:
                    screen[i] = 1.09085413633
                elif atom.type.atomic_number == 7:
                    screen[i] = 0.700147318409
                elif atom.type.atomic_number == 8:
                    screen[i] = 1.06557401132
                elif atom.type.atomic_number == 16:
                    screen[i] = 0.602256336067
        elif gb_model is GBn2:
            radii = _mbondi3_radii(self.atom_list)
            # Add non-optimized values as defaults
            alpha = [1.0 for i in self.atom_list]
            beta = [0.8 for i in self.atom_list]
            gamma = [4.85 for i in self.atom_list]
            screen = [0.5 for i in self.atom_list]
            for i, atom in enumerate(self.atom_list):
                if atom.type.atomic_number == 6:
                    screen[i] = 1.058554
                    alpha[i] = 0.733756
                    beta[i] = 0.506378
                    gamma[i] = 0.205844
                elif atom.type.atomic_number == 1:
                    screen[i] = 1.425952
                    alpha[i] = 0.788440
                    beta[i] = 0.798699
                    gamma[i] = 0.437334
                elif atom.type.atomic_number == 7:
                    screen[i] = 0.733599
                    alpha[i] = 0.503364
                    beta[i] = 0.316828
                    gamma[i] = 0.192915
                elif atom.type.atomic_number == 8:
                    screen[i] = 1.061039
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
                elif atom.type.atomic_number == 16:
                    screen[i] = -0.703469
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
        else:
            # Set the default screening parameters
            for i, atom in enumerate(self.atom_list):
                if atom.type.atomic_number == 1:
                    screen[i] = 0.85
                elif atom.type.atomic_number == 6:
                    screen[i] = 0.72
                elif atom.type.atomic_number == 7:
                    screen[i] = 0.79
                elif atom.type.atomic_number == 8:
                    screen[i] = 0.85
                elif atom.type.atomic_number == 9:
                    screen[i] = 0.88
                elif atom.type.atomic_number == 15:
                    screen[i] = 0.86
                elif atom.type.atomic_number == 16:
                    screen[i] = 0.96
                else:
                    screen[i] = 0.8
            # Determine which radii set we need
            if gb_model is OBC1 or gb_model is OBC2:
                radii = _mbondi2_radii(self.atom_list)
            elif gb_model is HCT:
                radii = _mbondi_radii(self.atom_list)

        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        radii = [x * length_conv for x in radii]

        if gb_model is GBn2:
            return zip(radii, screen, alpha, beta, gamma)
        return zip(radii, screen)

    def createSystem(self, params, nonbondedMethod=ff.NoCutoff,
                     nonbondedCutoff=1.0*u.nanometer,
                     switchDistance=0.0*u.nanometer,
                     constraints=None,
                     rigidWater=True,
                     implicitSolvent=None,
                     implicitSolventKappa=None,
                     implicitSolventSaltConc=0.0*u.moles/u.liter,
                     temperature=298.15*u.kelvin,
                     soluteDielectric=1.0,
                     solventDielectric=78.5,
                     removeCMMotion=True,
                     hydrogenMass=None,
                     ewaldErrorTolerance=0.0005,
                     flexibleConstraints=True,
                     verbose=False):
        """
        Construct an OpenMM System representing the topology described by the
        prmtop file. You MUST have loaded a parameter set into this PSF before
        calling createSystem. If not, AttributeError will be raised. ValueError
        is raised for illegal input.

        Parameters:
         -  params (CharmmParameterSet) The parameter set to use to parametrize
               this molecule
         -  nonbondedMethod (object=NoCutoff) The method to use for nonbonded
               interactions. Allowed values are NoCutoff, CutoffNonPeriodic,
               CutoffPeriodic, Ewald, or PME.
         -  nonbondedCutoff (distance=1*nanometer) The cutoff distance to use
               for nonbonded interactions.
         -  switchDistance (distance=0*nanometer) The distance at which the
               switching function is active for nonbonded interactions. If the
               switchDistance evaluates to boolean False (if it is 0), no
               switching function will be used. Illegal values will raise a
               ValueError
         -  constraints (object=None) Specifies which bonds or angles should be
               implemented with constraints. Allowed values are None, HBonds,
               AllBonds, or HAngles.
         -  rigidWater (boolean=True) If true, water molecules will be fully
               rigid regardless of the value passed for the constraints argument
         -  implicitSolvent (object=None) If not None, the implicit solvent
               model to use. Allowed values are HCT, OBC1, OBC2, or GBn
         -  implicitSolventKappa (float=None): Debye screening parameter to
               model salt concentrations in GB solvent.
         -  implicitSolventSaltConc (float=0.0*u.moles/u.liter): Salt
               concentration for GB simulations. Converted to Debye length
               `kappa'
         -  temperature (float=298.15*u.kelvin): Temperature used in the salt
               concentration-to-kappa conversion for GB salt concentration term
         -  soluteDielectric (float=1.0) The solute dielectric constant to use
               in the implicit solvent model.
         -  solventDielectric (float=78.5) The solvent dielectric constant to
               use in the implicit solvent model.
         -  removeCMMotion (boolean=True) If true, a CMMotionRemover will be
               added to the System.
         -  hydrogenMass (mass=None) The mass to use for hydrogen atoms bound to
               heavy atoms. Any mass added to a hydrogen is subtracted from the
               heavy atom to keep their total mass the same.
         -  ewaldErrorTolerance (float=0.0005) The error tolerance to use if the
               nonbonded method is Ewald or PME.
         -  flexibleConstraints (bool=True) Are our constraints flexible or not?
         -  verbose (bool=False) Optionally prints out a running progress report
        """
        # back up the dihedral list
        dihedral_list = self.dihedral_list
        # Load the parameter set
        self.loadParameters(params.condense())
        hasbox = self.topology.getUnitCellDimensions() is not None
        # Set the cutoff distance in nanometers
        cutoff = None
        if nonbondedMethod is not ff.NoCutoff:
            cutoff = nonbondedCutoff
            # Remove units from cutoff
            if u.is_quantity(cutoff):
                cutoff = cutoff.value_in_unit(u.nanometers)

        if nonbondedMethod not in (ff.NoCutoff, ff.CutoffNonPeriodic,
                                   ff.CutoffPeriodic, ff.Ewald, ff.PME):
            raise ValueError('Illegal value for nonbonded method')
        if not hasbox and nonbondedMethod in (ff.CutoffPeriodic,
                                              ff.Ewald, ff.PME):
            raise ValueError('Illegal nonbonded method for a '
                             'non-periodic system')
        if implicitSolvent not in (HCT, OBC1, OBC2, GBn, GBn2, None):
            raise ValueError('Illegal implicit solvent model choice.')
        if not constraints in (None, ff.HAngles, ff.HBonds, ff.AllBonds):
            raise ValueError('Illegal constraints choice')
      
        # Define conversion factors
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        _chmfrc = u.kilocalorie_per_mole/(u.angstrom*u.angstrom)
        _openmmfrc = u.kilojoule_per_mole/(u.nanometer*u.nanometer)
        bond_frc_conv = _chmfrc.conversion_factor_to(_openmmfrc)
        _chmfrc = u.kilocalorie_per_mole/(u.radians*u.radians)
        _openmmfrc = u.kilojoule_per_mole/(u.radians*u.radians)
        angle_frc_conv = _chmfrc.conversion_factor_to(_openmmfrc)
        dihe_frc_conv = u.kilocalorie_per_mole.conversion_factor_to(
                            u.kilojoule_per_mole)
        ene_conv = dihe_frc_conv
      
        # Create the system
        system = mm.System()
        if verbose: print('Adding particles...')
        for atom in self.atom_list:
            system.addParticle(atom.mass)
        # Set up the constraints
        if verbose and (constraints is not None and not rigidWater):
            print('Adding constraints...')
        if constraints in (ff.HBonds, ff.AllBonds, ff.HAngles):
            for bond in self.bond_list:
                if (bond.atom1.type.atomic_number != 1 and
                    bond.atom2.type.atomic_number != 1):
                    continue
                system.addConstraint(bond.atom1.idx, bond.atom2.idx,
                                     bond.bond_type.req*length_conv)
        if constraints in (ff.AllBonds, ff.HAngles):
            for bond in self.bond_list:
                if (bond.atom1.type.atomic_number == 1 or
                    bond.atom2.type.atomic_number == 1):
                    continue
                system.addConstraint(bond.atom1.idx, bond.atom2.idx,
                                     bond.bond_type.req*length_conv)
        if rigidWater and constraints is None:
            for bond in self.bond_list:
                if (bond.atom1.type.atomic_number != 1 and
                    bond.atom2.type.atomic_number != 1):
                    continue
                if (bond.atom1.residue.resname in WATNAMES and
                    bond.atom2.residue.resname in WATNAMES):
                    system.addConstraint(bond.atom1.idx, bond.atom2.idx,
                                         bond.bond_type.req*length_conv)
        # Add Bond forces
        if verbose: print('Adding bonds...')
        force = mm.HarmonicBondForce()
        force.setForceGroup(self.BOND_FORCE_GROUP)
        # See which, if any, energy terms we omit
        omitall = not flexibleConstraints and constraints is ff.AllBonds
        omith = omitall or (flexibleConstraints and constraints in
                            (ff.HBonds, ff.AllBonds, ff.HAngles))
        for bond in self.bond_list:
            if omitall: continue
            if omith and (bond.atom1.type.atomic_number == 1 or
                          bond.atom2.type.atomic_number == 1):
                continue
            force.addBond(bond.atom1.idx, bond.atom2.idx,
                          bond.bond_type.req*length_conv,
                          2*bond.bond_type.k*bond_frc_conv)
        system.addForce(force)
        # Add Angle forces
        if verbose: print('Adding angles...')
        force = mm.HarmonicAngleForce()
        force.setForceGroup(self.ANGLE_FORCE_GROUP)
        if constraints is ff.HAngles:
            num_constrained_bonds = system.getNumConstraints()
            atom_constraints = [[]] * system.getNumParticles()
            for i in range(num_constrained_bonds):
                c = system.getConstraintParameters(i)
                dist = c[2].value_in_unit(u.nanometer)
                atom_constraints[c[0]].append((c[1], dist))
                atom_constraints[c[1]].append((c[0], dist))
        for angle in self.angle_list:
            # Only constrain angles including hydrogen here
            if (angle.atom1.type.atomic_number != 1 and
                angle.atom2.type.atomic_number != 1 and
                angle.atom3.type.atomic_number != 1):
                continue
            if constraints is ff.HAngles:
                a1 = angle.atom1.atomic_number
                a2 = angle.atom2.atomic_number
                a3 = angle.atom3.atomic_number
                nh = int(a1==1) + int(a2==1) + int(a3==1)
                constrained = (nh >= 2 or (nh == 1 and a2 == 8))
            else:
                constrained = False # no constraints
            if constrained:
                l1 = l2 = None
                for bond in angle.atom2.bonds:
                    if bond.atom1 is angle.atom1 or bond.atom2 is angle.atom1:
                        l1 = bond.bond_type.req * length_conv
                    elif bond.atom1 is angle.atom3 or bond.atom2 is angle.atom3:
                        l2 = bond.bond_type.req * length_conv
                # Compute the distance between the atoms and add a constraint
                length = sqrt(l1*l1 + l2*l2 - 2*l1*l2*
                              cos(angle.angle_type.theteq))
                system.addConstraint(bond.atom1.idx, bond.atom2.idx, length)
            if flexibleConstraints or not constrained:
                force.addAngle(angle.atom1.idx, angle.atom2.idx,
                               angle.atom3.idx, angle.angle_type.theteq*pi/180,
                               2*angle.angle_type.k*angle_frc_conv)
        for angle in self.angle_list:
            # Already did the angles with hydrogen above. So skip those here
            if (angle.atom1.type.atomic_number == 1 or
                angle.atom2.type.atomic_number == 1 or
                angle.atom3.type.atomic_number == 1):
                continue
            force.addAngle(angle.atom1.idx, angle.atom2.idx,
                           angle.atom3.idx, angle.angle_type.theteq*pi/180,
                           2*angle.angle_type.k*angle_frc_conv)
        system.addForce(force)

        # Add the urey-bradley terms
        if verbose: print('Adding Urey-Bradley terms')
        force = mm.HarmonicBondForce()
        force.setForceGroup(self.UREY_BRADLEY_FORCE_GROUP)
        for ub in self.urey_bradley_list:
            force.addBond(ub.atom1.idx, ub.atom2.idx,
                          ub.ub_type.req*length_conv,
                          2*ub.ub_type.k*bond_frc_conv)
        system.addForce(force)

        # Add dihedral forces
        if verbose: print('Adding torsions...')
        force = mm.PeriodicTorsionForce()
        force.setForceGroup(self.DIHEDRAL_FORCE_GROUP)
        for tor in self.dihedral_list:
            force.addTorsion(tor.atom1.idx, tor.atom2.idx, tor.atom3.idx,
                             tor.atom4.idx, tor.dihedral_type.per,
                             tor.dihedral_type.phase*pi/180,
                             tor.dihedral_type.phi_k*dihe_frc_conv)
        system.addForce(force)

        if verbose: print('Adding impropers...')
        # Ick. OpenMM does not have an improper torsion class. Need to
        # construct one from CustomTorsionForce
        force = mm.CustomTorsionForce('k*(theta-theta0)^2')
        force.addPerTorsionParameter('k')
        force.addPerTorsionParameter('theta0')
        force.setForceGroup(self.IMPROPER_FORCE_GROUP)
        for imp in self.improper_list:
            force.addTorsion(imp.atom1.idx, imp.atom2.idx,
                             imp.atom3.idx, imp.atom4.idx,
                             (imp.improper_type.k*dihe_frc_conv,
                              imp.improper_type.phieq*pi/180)
            )
        system.addForce(force)

        if hasattr(self, 'cmap_list'):
            if verbose: print('Adding CMAP coupled torsions...')
            force = mm.CMAPTorsionForce()
            force.setForceGroup(self.CMAP_FORCE_GROUP)
            # First get the list of cmap maps we're going to use. Just store the
            # IDs so we have simple integer comparisons to do later
            cmap_type_list = []
            cmap_map = dict()
            for cmap in self.cmap_list:
                if not id(cmap.cmap_type) in cmap_type_list:
                    ct = cmap.cmap_type
                    cmap_type_list.append(id(ct))
                    # Our torsion correction maps need to go from 0 to 360
                    # degrees
                    grid = ct.grid.switch_range().T
                    m = force.addMap(ct.resolution, [x*ene_conv for x in grid])
                    cmap_map[id(ct)] = m
            # Now add in all of the cmaps
            for cmap in self.cmap_list:
                if cmap.consecutive:
                    id1, id2 = cmap.atom1.idx, cmap.atom2.idx
                    id3, id4 = cmap.atom3.idx, cmap.atom4.idx
                    id5, id6 = cmap.atom2.idx, cmap.atom3.idx
                    id7, id8 = cmap.atom4.idx, cmap.atom5.idx
                else:
                    id1, id2 = cmap.atom1.idx, cmap.atom2.idx
                    id3, id4 = cmap.atom3.idx, cmap.atom4.idx
                    id5, id6 = cmap.atom5.idx, cmap.atom6.idx
                    id7, id8 = cmap.atom7.idx, cmap.atom8.idx
                force.addTorsion(cmap_map[id(cmap.cmap_type)],
                                 id1, id2, id3, id4, id5, id6, id7, id8)
            system.addForce(force)
        # Add nonbonded terms now
        if verbose: print('Adding nonbonded interactions...')
        force = mm.NonbondedForce()
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        if not hasbox: # non-periodic
            if nonbondedMethod is ff.NoCutoff:
                force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
            elif nonbondedMethod is ff.CutoffNonPeriodic:
                if cutoff is None:
                    raise ValueError('No cutoff value specified')
                force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
                force.setCutoffDistance(cutoff)
            else:
                raise ValueError('Illegal nonbonded method for non-periodic '
                                 'system')

            # See if we need to use a switching function
            if switchDistance and nonbondedMethod is not ff.NoCutoff:
                # make sure it's legal
                if switchDistance >= nonbondedCutoff:
                    raise ValueError('switchDistance is too large compared '
                                     'to the cutoff!')
                    force.setUseSwitchingFunction(True)
                    force.setSwitchingDistance(switchDistance)

        else: # periodic
            # Set up box vectors (from inpcrd if available, or fall back to
            # prmtop definitions
            system.setDefaultPeriodicBoxVectors(*self.box_vectors)

            # Set cutoff
            if cutoff is None:
                # Compute cutoff automatically
                box = self.boxLengths
                min_box_width = min((box[0]/u.nanometers,
                                     box[1]/u.nanometers,
                                     box[2]/u.nanometers))
                CLEARANCE_FACTOR = 0.97
                cutoff = u.Quantity((min_box_width*CLEARANCE_FACTOR)/2.0,
                                    u.nanometers)
            if nonbondedMethod is not ff.NoCutoff:
                force.setCutoffDistance(cutoff)

            # Set nonbonded method.
            if nonbondedMethod is ff.NoCutoff:
                force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
            elif nonbondedMethod is ff.CutoffNonPeriodic:
                force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
            elif nonbondedMethod is ff.CutoffPeriodic:
                force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            elif nonbondedMethod is ff.Ewald:
                force.setNonbondedMethod(mm.NonbondedForce.Ewald)
            elif nonbondedMethod is ff.PME:
                force.setNonbondedMethod(mm.NonbondedForce.PME)
            else:
                raise ValueError('Cutoff method is not understood')

            # See if we need to use a switching function
            if switchDistance and nonbondedMethod is not ff.NoCutoff:
                # make sure it's legal
                if switchDistance >= nonbondedCutoff:
                    raise ValueError('switchDistance is too large compared '
                                     'to the cutoff!')
                    force.setUseSwitchingFunction(True)
                    force.setSwitchingDistance(switchDistance)

            if ewaldErrorTolerance is not None:
                force.setEwaldErrorTolerance(ewaldErrorTolerance)

        # Add per-particle nonbonded parameters (LJ params)
        sigma_scale = 2**(-1/6) * 2
        for i, atm in enumerate(self.atom_list):
            force.addParticle(atm.charge, sigma_scale*atm.type.rmin*length_conv,
                              abs(atm.type.epsilon*ene_conv))

        # Add 1-4 interactions
        excluded_atom_pairs = set() # save these pairs so we don't zero them out
        sigma_scale = 2**(-1/6)
        for tor in self.dihedral_list:
            # First check to see if atoms 1 and 4 are already excluded because
            # they are 1-2 or 1-3 pairs (would happen in 6-member rings or
            # fewer). Then check that they're not already added as exclusions
            if tor.atom1 in tor.atom4.bond_partners: continue
            if tor.atom1 in tor.atom4.angle_partners: continue
            key = min((tor.atom1.idx, tor.atom4.idx),
                      (tor.atom4.idx, tor.atom1.idx))
            if key in excluded_atom_pairs: continue # multiterm...
            charge_prod = (tor.atom1.charge * tor.atom4.charge)
            epsilon = (sqrt(abs(tor.atom1.type.epsilon_14) * ene_conv *
                            abs(tor.atom4.type.epsilon_14) * ene_conv))
            sigma = (tor.atom1.type.rmin_14 + tor.atom4.type.rmin_14) * (
                     length_conv * sigma_scale)
            force.addException(tor.atom1.idx, tor.atom4.idx,
                               charge_prod, sigma, epsilon)
            excluded_atom_pairs.add(
                    min((tor.atom1.idx, tor.atom4.idx),
                        (tor.atom4.idx, tor.atom1.idx))
            )

        # Add excluded atoms
        for atom in self.atom_list:
            # Exclude all bonds and angles
            for atom2 in atom.bond_partners:
                if atom2.idx > atom.idx:
                    force.addException(atom.idx, atom2.idx, 0.0, 0.1, 0.0)
            for atom2 in atom.angle_partners:
                if atom2.idx > atom.idx:
                    force.addException(atom.idx, atom2.idx, 0.0, 0.1, 0.0)
            for atom2 in atom.dihedral_partners:
                if atom2.idx <= atom.idx: continue
                if ((atom.idx, atom2.idx) in excluded_atom_pairs):
                    continue
                force.addException(atom.idx, atom2.idx, 0.0, 0.1, 0.0)
        system.addForce(force)

        # Add GB model if we're doing one
        if implicitSolvent is not None:
            if verbose: print('Adding GB parameters...')
            gb_parms = self._get_gb_params(implicitSolvent)

            # If implicitSolventKappa is None, compute it from salt
            # concentration
            if implicitSolventKappa is None:
                if u.is_quantity(implicitSolventSaltConc):
                    sc = implicitSolventSaltConc.value_in_unit(u.moles/u.liter)
                    implicitSolventSaltConc = sc
                if u.is_quantity(temperature):
                    temperature = temperature.value_in_unit(u.kelvin)
                # The constant is 1 / sqrt( epsilon_0 * kB / (2 * NA * q^2 *
                # 1000) ) where NA is avogadro's number, epsilon_0 is the
                # permittivity of free space, q is the elementary charge (this
                # number matches sander/pmemd's kappa conversion factor)
                implicitSolventKappa = 50.33355 * sqrt(implicitSolventSaltConc /
                                                solventDielectric / temperature)
                # Multiply by 0.73 to account for ion exclusions, and multiply
                # by 10 to convert to 1/nm from 1/angstroms
                implicitSolventKappa *= 7.3
            elif implicitSolvent is None:
                implicitSolventKappa = 0.0

            if u.is_quantity(implicitSolventKappa):
                implicitSolventKappa = implicitSolventKappa.value_in_unit(
                                            (1.0/u.nanometer).unit)
            if implicitSolvent is HCT:
                gb = GBSAHCTForce(solventDielectric, soluteDielectric, None,
                                  cutoff, kappa=implicitSolventKappa)
            elif implicitSolvent is OBC1:
                gb = GBSAOBC1Force(solventDielectric, soluteDielectric, None,
                                   cutoff, kappa=implicitSolventKappa)
            elif implicitSolvent is OBC2:
                gb = GBSAOBC2Force(solventDielectric, soluteDielectric, None,
                                   cutoff, kappa=implicitSolventKappa)
            elif implicitSolvent is GBn:
                gb = GBSAGBnForce(solventDielectric, soluteDielectric, None,
                                  cutoff, kappa=implicitSolventKappa)
            elif implicitSolvent is GBn2:
                gb = GBSAGBn2Force(solventDielectric, soluteDielectric, None,
                                   cutoff, kappa=implicitSolventKappa)
            for i, atom in enumerate(self.atom_list):
                gb.addParticle([atom.charge] + list(gb_parms[i]))
            # Set cutoff method
            if nonbondedMethod is ff.NoCutoff:
                gb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
            elif nonbondedMethod is ff.CutoffNonPeriodic:
                gb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
                gb.setCutoffDistance(cutoff)
            elif nonbondedMethod is ff.CutoffPeriodic:
                gb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
                gb.setCutoffDistance(cutoff)
            else:
                raise ValueError('Illegal nonbonded method for use with GBSA')
            gb.setForceGroup(self.GB_FORCE_GROUP)
            system.addForce(gb)
            force.setReactionFieldDielectric(1.0) # applies to NonbondedForce

        # See if we repartition the hydrogen masses
        if hydrogenMass is not None:
            for bond in self.bond_list:
                # Only take the ones with at least one hydrogen
                if (bond.atom1.type.atomic_number != 1 and
                    bond.atom2.type.atomic_number != 1):
                    continue
                atom1, atom2 = bond.atom1, bond.atom2
                if atom1.type.atomic_number == 1:
                    atom1, atom2 = atom2, atom1 # now atom2 is hydrogen for sure
                if atom1.type.atomic_number != 1:
                    transfer_mass = hydrogenMass - atom2.mass
                    new_mass1 = (system.getParticleMass(atom1.idx) -
                                 transfer_mass)
                    system.setParticleMass(atom2.idx, hydrogenMass)
                    system.setParticleMass(atom1.idx, new_mass1)
        # See if we want to remove COM motion
        if removeCMMotion:
            system.addForce(mm.CMMotionRemover())

        # Cache our system for easy access
        self._system = system

        # Restore the dihedral list to allow reparametrization later
        self.dihedral_list = dihedral_list

        return system

    @property
    def system(self):
        """
        Return the cached system class -- it needs to be initialized via
        "createSystem" first!
        """
        try:
            return self._system
        except AttributeError:
            raise AttributeError('You must initialize the system with '
                                 'createSystem before accessing the cached '
                                 'object.')

    @property
    def positions(self):
        """
        Return the cached positions or create new ones from the atoms
        """
        try:
            if len(self._positions) == len(self.atom_list):
                return self._positions
        except AttributeError:
            pass

        self._positions = tuple([Vec3(a.xx, a.xy, a.xz)
                               for a in self.atom_list]) * u.angstroms
        return self._positions

    @positions.setter
    def positions(self, stuff):
        """
        Replace the cached positions and the positions of each atom. If no units
        are applied to "stuff", it is assumed to be Angstroms.
        """
        if not u.is_quantity(stuff):
            # Assume this is Angstroms
            stuff *= u.angstroms

        # If we got a 1-D array, reshape it into an natom list of Vec3's
        if len(stuff) == len(self.atom_list) * 3:
            stuff = [Vec3(stuff[i*3], stuff[i*3+1], stuff[i*3+2])
                     for i in range(len(self.atom_list))]
        self._positions = stuff
        for atom, pos in zip(self.atom_list, stuff):
            atom.xx, atom.xy, atom.xz = pos.value_in_unit(u.angstrom)

    @property
    def velocities(self):
        """ Same as for positions, but for velocities """
        try:
            if len(self._velocities) == len(self.atom_list):
                return self._velocities
        except AttributeError:
            pass

        self._velocities = tuple([Vec3(a.vx, a.vy, a.vz)
                    for a in self.atom_list]) * (u.angstroms/u.picosecond) 
        return self._velocities

    @property
    def boxLengths(self):
        """ Return tuple of 3 units """
        if self.box_vectors is not None:
            return (self.box_vectors[0][0], self.box_vectors[0][1],
                    self.box_vectors[0][2])
        return None

    @boxLengths.setter
    def boxLengths(self, stuff):
        raise RuntimeError('Use setBox to set a box with lengths and angles '
                           'or set the boxVectors attribute with box vectors')
    
    @property
    def boxVectors(self):
        """ Return the box vectors """
        return self.box_vectors

    @boxVectors.setter
    def boxVectors(self, stuff):
        """ Sets the box vectors """
        self.box_vectors = stuff

    def deleteCmap(self):
        """ Deletes the CMAP terms from the CHARMM PSF """
        self.cmap_list = TrackedList()

def _box_vectors_from_lengths_angles(a, b, c, alpha, beta, gamma):
    """
    This method takes the lengths and angles from a unit cell and creates unit
    cell vectors.

    Parameters:
        - a (unit, dimension length): Length of the first vector
        - b (unit, dimension length): Length of the second vector
        - c (unit, dimension length): Length of the third vector
        - alpha (float): Angle between b and c in degrees
        - beta (float): Angle between a and c in degrees
        - gamma (float): Angle between a and b in degrees

    Returns:
        Tuple of box vectors (as Vec3 instances)
   """
    if not (u.is_quantity(a) and u.is_quantity(b) and u.is_quantity(c)):
        raise TypeError('a, b, and c must be units of dimension length')
    if u.is_quantity(alpha): alpha = alpha.value_in_unit(u.degree)
    if u.is_quantity(beta): beta = beta.value_in_unit(u.degree)
    if u.is_quantity(gamma): gamma = gamma.value_in_unit(u.degree)
    a = a.value_in_unit(u.angstrom)
    b = b.value_in_unit(u.angstrom)
    c = c.value_in_unit(u.angstrom)

    if alpha <= 2 * pi and beta <= 2 * pi and gamma <= 2 * pi:
        raise ValueError('box angles must be given in degrees')

    alpha *= pi / 180
    beta *= pi / 180
    gamma *= pi / 180

    av = Vec3(a, 0.0, 0.0) * u.angstrom
    bx = b * cos(gamma)
    by = b * sin(gamma)
    bz = 0.0
    cx = c * cos(beta)
    cy = c * (cos(alpha) - cos(beta) * cos(gamma))
    cz = sqrt(c * c - cx * cx - cy * cy)
   
    # Make sure any components that are close to zero are set to zero exactly
    if abs(bx) < TINY: bx = 0.0
    if abs(by) < TINY: by = 0.0
    if abs(cx) < TINY: cx = 0.0
    if abs(cy) < TINY: cy = 0.0
    if abs(cz) < TINY: cz = 0.0

    bv = Vec3(bx, by, bz) * u.angstrom
    cv = Vec3(cx, cy, cz) * u.angstrom

    return (av, bv, cv)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def set_molecules(atom_list):
    """
    Correctly sets the molecularity of the system based on connectivity
    """
    from sys import setrecursionlimit, getrecursionlimit
    # Since we use a recursive function here, we make sure that the recursion
    # limit is large enough to handle the maximum possible recursion depth we'll
    # need (NATOM). We don't want to shrink it, though, since we use list
    # comprehensions in list constructors in some places that have an implicit
    # (shallow) recursion, therefore, reducing the recursion limit too much here
    # could raise a recursion depth exceeded exception during a _Type/Atom/XList
    # creation. Therefore, set the recursion limit to the greater of the current
    # limit or the number of atoms
    setrecursionlimit(max(len(atom_list), getrecursionlimit()))

    # Unmark all atoms so we can track which molecule each goes into
    atom_list.unmark()

    # The molecule "ownership" list
    owner = []
    # The way I do this is via a recursive algorithm, in which
    # the "set_owner" method is called for each bonded partner an atom
    # has, which in turn calls set_owner for each of its partners and 
    # so on until everything has been assigned.
    molecule_number = 1 # which molecule number we are on
    for i in range(len(atom_list)):
        # If this atom has not yet been "owned", make it the next molecule
        # However, we only increment which molecule number we're on if 
        # we actually assigned a new molecule (obviously)
        if not atom_list[i].marked:
            tmp = [i]
            _set_owner(atom_list, tmp, i, molecule_number)
            # Make sure the atom indexes are sorted
            tmp.sort()
            owner.append(tmp)
            molecule_number += 1
    return owner

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _set_owner(atom_list, owner_array, atm, mol_id):
    """ Recursively sets ownership of given atom and all bonded partners """
    atom_list[atm].marked = mol_id
    for partner in atom_list[atm].bond_partners:
        if not partner.marked:
            owner_array.append(partner.idx)
            _set_owner(atom_list, owner_array, partner.idx, mol_id)
        elif partner.marked != mol_id:
            raise MoleculeError('Atom %d in multiple molecules' % 
                                partner.idx)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Routines that set the necessary radii lists based on a list of atoms with
# their connectivities

def _bondi_radii(atom_list):
    """ Sets the bondi radii """
    radii = [0.0 for atom in atom_list]
    for i, atom in enumerate(atom_list):
        if atom.type.atomic_number == 6:
            radii[i] = 1.7
        elif atom.type.atomic_number == 1:
            radii[i] = 1.2
        elif atom.type.atomic_number == 7:
            radii[i] = 1.55
        elif atom.type.atomic_number == 8:
            radii[i] = 1.5
        elif atom.type.atomic_number == 9:
            radii[i] = 1.5
        elif atom.type.atomic_number == 14:
            radii[i] = 2.1
        elif atom.type.atomic_number == 15:
            radii[i] = 1.85
        elif atom.type.atomic_number == 16:
            radii[i] = 1.8
        elif atom.type.atomic_number == 17:
            radii[i] = 1.5
        else:
            radii[i] = 1.5
    return radii  # converted to nanometers above

def _mbondi_radii(atom_list):
    """ Sets the mbondi radii """
    radii = [0.0 for atom in atom_list]
    for i, atom in enumerate(atom_list):
        # Radius of H atom depends on element it is bonded to
        if atom.type.atomic_number == 1:
            bondeds = list(atom.bond_partners)
            if bondeds[0].type.atomic_number in (6, 7): # C or N
                radii[i] = 1.3
            elif bondeds[0].type.atomic_number in (8, 16): # O or S
                radii[i] = 0.8
            else:
                radii[i] = 1.2
        # Radius of C atom depends on what type it is
        elif atom.type.atomic_number == 6:
            radii[i] = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.type.atomic_number == 7:
            radii[i] = 1.55
        elif atom.type.atomic_number == 8:
            radii[i] = 1.5
        elif atom.type.atomic_number == 9:
            radii[i] = 1.5
        elif atom.type.atomic_number == 14:
            radii[i] = 2.1
        elif atom.type.atomic_number == 15:
            radii[i] = 1.85
        elif atom.type.atomic_number == 16:
            radii[i] = 1.8
        elif atom.type.atomic_number == 17:
            radii[i] = 1.5
        else:
            radii[i] = 1.5
    return radii  # converted to nanometers above

def _mbondi2_radii(atom_list):
    """ Sets the mbondi2 radii """
    radii = [0.0 for atom in atom_list]
    for i, atom in enumerate(atom_list):
        # Radius of H atom depends on element it is bonded to
        if atom.type.atomic_number == 1:
            if atom.bond_partners[0].type.atomic_number == 7:
                radii[i] = 1.3
            else:
                radii[i] = 1.2
        # Radius of C atom depends on what type it is
        elif atom.type.atomic_number == 6:
            radii[i] = 1.7
        # All other elements have fixed radii for all types/partners
        elif atom.type.atomic_number == 7:
            radii[i] = 1.55
        elif atom.type.atomic_number == 8:
            radii[i] = 1.5
        elif atom.type.atomic_number == 9:
            radii[i] = 1.5
        elif atom.type.atomic_number == 14:
            radii[i] = 2.1
        elif atom.type.atomic_number == 15:
            radii[i] = 1.85
        elif atom.type.atomic_number == 16:
            radii[i] = 1.8
        elif atom.type.atomic_number == 17:
            radii[i] = 1.5
        else:
            radii[i] = 1.5
    return radii  # Converted to nanometers above

def _mbondi3_radii(atom_list):
    """ Sets the mbondi3 radii """
    radii = _mbondi2_radii(atom_list)
    for i, atom in enumerate(atom_list):
        # Adjust OE (GLU), OD (ASP) and HH/HE (ARG)
        if atom.residue.resname in ('GLU', 'ASP'):
            if atom.name.startswith('OE') or atom.name.startswith('OD'):
                radii[i] = 1.4
        elif atom.residue.resname == 'ARG':
            if atom.name.startswith('HH') or atom.name.startswith('HE'):
                radii[i] = 1.17
        # Adjust carboxylate O radii on C-termini
        if atom.name == 'OXT':
            radii[i] = 1.4
    return radii  # Converted to nanometers above

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == '__main__':
    import doctest
    doctest.testmod()
