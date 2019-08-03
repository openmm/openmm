"""
This module contains classes for parsing and processing CHARMM parameter,
topology, and stream files. It only extracts atom properties from the
topology files and extracts all parameters from the parameter files

This file is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of Biological
Structures at Stanford, funded under the NIH Roadmap for Medical Research,
grant U54 GM072970. See https://simtk.org.  This code was originally part of
the ParmEd program and was ported for use with OpenMM.

Copyright (c) 2014 the Authors

Author: Jason M. Swails
Contributors: Jing Huang
Date: Sep. 17, 2014

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
import os
from simtk.openmm.app.internal.charmm._charmmfile import (
            CharmmFile, CharmmStreamFile)
from simtk.openmm.app.internal.charmm.topologyobjects import (
            AtomType, BondType, AngleType, DihedralType, ImproperType, CmapType,
            UreyBradleyType, NoUreyBradley)
from simtk.openmm.app.internal.charmm.exceptions import CharmmFileError
from simtk.openmm.app.element import Element, get_by_symbol
import simtk.unit as u
import warnings

class CharmmParameterSet(object):
    """
    Stores a parameter set defined by CHARMM files. It stores the equivalent of
    the information found in the MASS section of the CHARMM topology file
    (TOP/RTF) and all of the information in the parameter files (PAR)

    Parameters
    ----------
    filenames : List of topology, parameter, and stream files to load into the parameter set.
        The following file type suffixes are recognized. Unrecognized file types raise a TypeError
          * .rtf, .top -- Residue topology file
          * .par, .prm -- Parameter file
          * .str -- Stream file
          * .inp -- If "par" is in the file name, it is a parameter file. If
                    "top" is in the file name, it is a topology file. Otherwise,
                    raise TypeError

    Attributes
    ----------
    All type lists are dictionaries whose keys are tuples (with however
    many elements are needed to define that type of parameter). The types
    that can be in any order are SORTED.

    - atom_types_str
    - atom_types_int
    - atom_types_tuple
    - bond_types
    - angle_types
    - urey_bradley_types
    - dihedral_types
    - improper_types
    - cmap_types
    - nbfix_types

    The dihedral types can be multiterm, so the values for each dict key is
    actually a list of DihedralType instances. The atom_types are dicts that
    match the name (str), number (int), or (name, number) tuple (tuple) to
    the atom type. The tuple is guaranteed to be the most robust, although
    when only the integer or string is available the other dictionaries are
    helpful

    Examples
    --------
    >>> params = CharmmParameterSet('charmm22.top', 'charmm22.par', 'file.str')
    """

    @staticmethod
    def _convert(data, type, msg=''):
        """
        Converts a data type to a desired type, raising CharmmFileError if it
        fails
        """
        try:
            return type(data)
        except ValueError:
            raise CharmmFileError('Could not convert %s to %s' % (msg, type))

    def __init__(self, *args, **kwargs):
        # Instantiate the list types
        self.atom_types_str = dict()
        self.atom_types_int = dict()
        self.atom_types_tuple = dict()
        self.bond_types = dict()
        self.angle_types = dict()
        self.urey_bradley_types = dict()
        self.dihedral_types = dict()
        self.improper_types = dict()
        self.cmap_types = dict()
        self.nbfix_types = dict()
        self.nbthole_types = dict()
        self.parametersets = []

        # Load all of the files
        tops, pars, strs = [], [], []
        for arg in args:
            if arg.endswith('.rtf') or arg.endswith('.top'):
                tops.append(arg)
            elif arg.endswith('.par') or arg.endswith('.prm'):
                pars.append(arg)
            elif arg.endswith('.str'):
                strs.append(arg)
            elif arg.endswith('.inp'):
                # Only consider the file name (since the directory is likely
                # "toppar" and will screw up file type detection)
                fname = os.path.split(arg)[1]
                if 'par' in fname:
                    pars.append(arg)
                elif 'top' in fname:
                    tops.append(arg)
                else:
                    raise TypeError('Unrecognized file type: %s' % arg)
            else:
                raise TypeError('Unrecognized file type: %s' % arg)

        permissive=kwargs.pop("permissive", False)
        if len(kwargs):
            raise TypeError('Unrecognised named argument')

        for top in tops: self.readTopologyFile(top)
        for par in pars: self.readParameterFile(par, permissive=permissive)
        for strf in strs: self.readStreamFile(strf)

    @classmethod
    def loadSet(cls, tfile=None, pfile=None, sfiles=[], permissive=False):
        """
        Instantiates a CharmmParameterSet from a Topology file and a Parameter
        file (or just a Parameter file if it has all information)

        Parameters
        -----------
        tfile : str
            Name of the Topology (RTF/TOP) file
        pfile : str
            Name of the Parameter (PAR) file
        sfiles : list of str
            List or tuple of stream (STR) file names.
        permissive : bool=False
            Accept non-bonbded parameters for undefined atom types

        Returns
        -------
        CharmmParameterSet
            New CharmmParameterSet populated with the parameters found in the
            provided files.

        Notes
        -----
        The RTF file is read first (if provided), followed by the PAR file,
        followed by the list of stream files (in the order they are
        provided). Parameters in each stream file will overwrite those that
        came before (or simply append to the existing set if they are
        different)
        """
        inst = cls()
        if tfile is not None:
            inst.readTopologyFile(tfile)
        if pfile is not None:
            inst.readParameterFile(pfile, permissive=permissive)
        if isinstance(sfiles, str):
            # The API docstring requests a list, but allow for users to pass a
            # string with a single filename instead
            inst.readStreamFile(sfiles)
        elif sfiles is not None:
            for sfile in sfiles:
                inst.readStreamFile(sfile)
        return inst

    def readParameterFile(self, pfile, permissive=False):
        """Reads all of the parameters from a parameter file. Versions 36 and later
        of the CHARMM force field files have an ATOMS section defining all of
        the atom types.  Older versions need to load this information from the
        RTF/TOP files.

        Parameters
        ----------
        pfile : str
            Name of the CHARMM PARameter file to read
        permissive : bool
            Accept non-bonbded parameters for undefined atom types (default:
            False).

        Notes
        -----
        The atom types must all be loaded by the end of this routine. Either
        supply a PAR file with atom definitions in them or read in a RTF/TOP
        file first. Failure to do so will result in a raised RuntimeError.
        """
        conv = CharmmParameterSet._convert
        if isinstance(pfile, str):
            own_handle = True
            f = CharmmFile(pfile)
        else:
            own_handle = False
            f = pfile
        # What section are we parsing?
        section = None
        # The current cmap we are building (these span multiple lines)
        current_cmap = None
        current_cmap_data = []
        current_cmap_res = 0
        nonbonded_types = dict() # Holder
        parameterset = None
        read_first_nonbonded = False
        for line in f:
            line = line.strip()
            if not line:
                # This is a blank line
                continue
            if parameterset is None and line.strip().startswith('*>>'):
                parameterset = line.strip()[1:78]
                continue
            # Set section if this is a section header
            if line.startswith('ATOMS'):
                section = 'ATOMS'
                continue
            if line.startswith('BONDS'):
                section = 'BONDS'
                continue
            if line.startswith('ANGLES'):
                section = 'ANGLES'
                continue
            if line.startswith('DIHEDRALS'):
                section = 'DIHEDRALS'
                continue
            if line.startswith('IMPROPER'):
                section = 'IMPROPER'
                continue
            if line.startswith('CMAP'):
                section = 'CMAP'
                continue
            if line.startswith('NONBONDED'):
                read_first_nonbonded = False
                section = 'NONBONDED'
                continue
            if line.startswith('NBFIX'):
                section = 'NBFIX'
                continue
            if line.startswith('THOLE'):
                section = 'NBTHOLE'
                continue
            if line.startswith('HBOND'):
                section = None
                continue
            # It seems like files? sections? can be terminated with 'END'
            if line.startswith('END'): # should this be case-insensitive?
                section = None
                continue
            # If we have no section, skip
            if section is None: continue
            # Now handle each section specifically
            if section == 'ATOMS':
                if not line.startswith('MASS'): continue # Should this happen?
                words = line.split()
                try:
                    idx = conv(words[1], int, 'atom type')
                    name = words[2]
                    mass = conv(words[3], float, 'atom mass')
                except IndexError:
                    raise CharmmFileError('Could not parse MASS section.')
                # The parameter file might or might not have an element name
                try:
                    elem = words[4]
                    atomic_number = get_by_symbol(elem).atomic_number
                except (IndexError, KeyError):
                    # Figure it out from the mass
                    masselem = Element.getByMass(mass)
                    if masselem is None:
                        atomic_number = 0 # Extra point or something
                    else:
                        atomic_number = masselem.atomic_number
                atype = AtomType(name=name, number=idx, mass=mass,
                                 atomic_number=atomic_number)
                self.atom_types_str[atype.name] = atype
                self.atom_types_int[atype.number] = atype
                self.atom_types_tuple[(atype.name, atype.number)] = atype
                continue
            if section == 'BONDS':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    k = conv(words[2], float, 'bond force constant')
                    req = conv(words[3], float, 'bond equilibrium dist')
                except IndexError:
                    raise CharmmFileError('Could not parse bonds.')
                key = (min(type1, type2), max(type1, type2))
                self.bond_types[key] = BondType(k, req)
                continue
            if section == 'ANGLES':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    k = conv(words[3], float, 'angle force constant')
                    theteq = conv(words[4], float, 'angle equilibrium value')
                except IndexError:
                    raise CharmmFileError('Could not parse angles.')
                key = (min(type1, type3), type2, max(type1, type3))
                self.angle_types[key] = AngleType(k, theteq)
                # See if we have a urey-bradley
                try:
                    ubk = conv(words[5], float, 'Urey-Bradley force constant')
                    ubeq = conv(words[6], float, 'Urey-Bradley equil. value')
                    ubtype = UreyBradleyType(ubk, ubeq)
                except IndexError:
                    ubtype = NoUreyBradley
                self.urey_bradley_types[key] = ubtype
                continue
            if section == 'DIHEDRALS':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    type4 = words[3]
                    k = conv(words[4], float, 'dihedral force constant')
                    n = conv(words[5], float, 'dihedral periodicity')
                    phase = conv(words[6], float, 'dihedral phase')
                except IndexError:
                    raise CharmmFileError('Could not parse dihedrals.')
                # Torsion can be in either direction. Sort by end groups first,
                # then sort by middle 2
                if type1 < type4:
                    key = (type1, type2, type3, type4)
                elif type1 > type4:
                    key = (type4, type3, type2, type1)
                else:
                    # OK, we need to sort by the middle atoms now
                    if type2 < type3:
                        key = (type1, type2, type3, type4)
                    else:
                        key = (type4, type3, type2, type1)
                # See if this is a second (or more) term of the dihedral group
                # that's already present.
                dihedral = DihedralType(k, n, phase)
                if key in self.dihedral_types:
                    # See if the existing dihedral type list has a term with
                    # the same periodicity -- If so, replace it
                    replaced = False
                    for i, dtype in enumerate(self.dihedral_types[key]):
                        if dtype.per == dihedral.per:
                            # Replace. Warn if they are different
                            if dtype != dihedral:
                                warnings.warn('Replacing dihedral %r with %r' %
                                              (dtype, dihedral))
                            self.dihedral_types[key]
                            replaced = True
                            break
                    if not replaced:
                        self.dihedral_types[key].append(dihedral)
                else: # key not present
                    self.dihedral_types[key] = [dihedral]
                continue
            if section == 'IMPROPER':
                words = line.split()
                try:
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    type4 = words[3]
                    k = conv(words[4], float, 'improper force constant')
                    theteq = conv(words[5], float, 'improper equil. value')
                except IndexError:
                    raise CharmmFileError('Could not parse dihedrals.')
                # If we have a 7th column, that is the real psi0 (and the 6th
                # is just a dummy 0)
                try:
                    tmp = conv(words[6], float, 'improper equil. value')
                    theteq = tmp
                except IndexError:
                    pass # Do nothing
                # Improper types seem not to have the central atom defined in
                # the first place, so just have the key a fully sorted list. We
                # still depend on the PSF having properly ordered improper atoms
                key = tuple(sorted([type1, type2, type3, type4]))
                self.improper_types[key] = ImproperType(k, theteq)
                continue
            if section == 'CMAP':
                # This is the most complicated part, since cmap parameters span
                # many lines. We won't do much error catching here.
                words = line.split()
                try:
                    holder = [float(w) for w in words]
                    current_cmap_data.extend(holder)
                except ValueError:
                    # We assume this is a definition of a new CMAP, so
                    # terminate the last CMAP if applicable
                    if current_cmap is not None:
                        # We have a map to terminate
                        ty = CmapType(current_cmap_res, current_cmap_data)
                        self.cmap_types[current_cmap] = ty
                    try:
                        type1 = words[0]
                        type2 = words[1]
                        type3 = words[2]
                        type4 = words[3]
                        type5 = words[4]
                        type6 = words[5]
                        type7 = words[6]
                        type8 = words[7]
                        res = conv(words[8], int, 'CMAP resolution')
                    except IndexError:
                        raise CharmmFileError('Could not parse CMAP data.')
                    # order the torsions independently
                    k1 = [type1, type2, type3, type4]
                    k2 = [type4, type3, type2, type1]
                    key1 = min(k1, k2)
                    k1 = [type5, type6, type7, type8]
                    k2 = [type8, type7, type6, type5]
                    key2 = min(k1, k2)
                    current_cmap = tuple(key1 + key2)
                    current_cmap_res = res
                    current_cmap_data = []
                continue
            if section == 'NONBONDED':
                # Now get the nonbonded values
                words = line.split()
                try:
                    atype = words[0]
                    # 1st column is ignored
                    epsilon = conv(words[2], float, 'vdW epsilon term')
                    rmin = conv(words[3], float, 'vdW Rmin/2 term')
                except IndexError:
                    # If we haven't read our first nonbonded term yet, we may
                    # just be parsing the settings that should be used. So
                    # soldier on
                    if not read_first_nonbonded: continue
                    raise CharmmFileError('Could not parse nonbonded terms.')
                except CharmmFileError as e:
                    if not read_first_nonbonded: continue
                    raise CharmmFileError(str(e))
                else:
                    # OK, we've read our first nonbonded section for sure now
                    read_first_nonbonded = True
                # See if we have 1-4 parameters
                try:
                    # 4th column is ignored
                    eps14 = conv(words[5], float, '1-4 vdW epsilon term')
                    rmin14 = conv(words[6], float, '1-4 vdW Rmin/2 term')
                except IndexError:
                    eps14 = rmin14 = None
                nonbonded_types[atype] = [epsilon, rmin, eps14, rmin14]
                continue
            if section == 'NBFIX':
                words = line.split()
                try:
                    at1 = words[0]
                    at2 = words[1]
                    emin = abs(conv(words[2], float, 'NBFIX Emin'))
                    rmin = conv(words[3], float, 'NBFIX Rmin')
                    try:
                        emin14 = abs(conv(words[4], float, 'NBFIX Emin 1-4'))
                        rmin14 = conv(words[5], float, 'NBFIX Rmin 1-4')
                    except IndexError:
                        emin14 = rmin14 = None
                    try:
                        self.atom_types_str[at1].add_nbfix(at2, rmin, emin,
                                                           rmin14, emin14)
                        self.atom_types_str[at2].add_nbfix(at1, rmin, emin,
                                                           rmin14, emin14)
                    except KeyError:
                        # Some stream files define NBFIX terms with an atom that
                        # is defined in another toppar file that does not
                        # necessarily have to be loaded. As a result, not every
                        # NBFIX found here will necessarily need to be applied.
                        # If we can't find a particular atom type, don't bother
                        # adding that nbfix and press on
                        pass
                except IndexError:
                    raise CharmmFileError('Could not parse NBFIX terms.')
                self.nbfix_types[(min(at1, at2), max(at1, at2))] = (emin, rmin)
                continue
            # Here parse the possible nbthole section
            if section == 'NBTHOLE':
                words = line.split()
                try:
                    at1 = words[0]
                    at2 = words[1]
                    nbt = abs(conv(words[2], float, 'NBTHOLE a'))
                    try:
                        self.atom_types_str[at1].add_nbthole(at2, nbt)
                        self.atom_types_str[at2].add_nbthole(at1, nbt)
                    except KeyError:
                        pass
                except IndexError:
                    raise CharmmFileError('Could not parse NBTHOLE terms.')
                self.nbthole_types[(min(at1, at2), max(at1, at2))] = (nbt)
        # If there were any CMAP terms stored in the parameter set, the last one
        # defined will not have been added to the set. Add it now.
        if current_cmap is not None:
            ty = CmapType(current_cmap_res, current_cmap_data)
            self.cmap_types[current_cmap] = ty

        # If in permissive mode create an atomtype for every type used in
        # the nonbonded parameters. This is a work-around for when all that's
        # available is a CHARMM22 inp file, which has no ATOM/MASS fields

        if permissive:
            try:
               idx = max(self.atom_types_int.keys())+1000
            except ValueError:
               idx = 10000
            for key in nonbonded_types:
                if not key in self.atom_types_str:
                    atype =AtomType(name=key, number=idx, mass= float('NaN'), atomic_number= 1 )
                    self.atom_types_str[key] = atype
                    self.atom_types_int[idx] = atype
                    idx=idx+1

        # Now we're done. Load the nonbonded types into the relevant AtomType
        # instances. In order for this to work, all keys in nonbonded_types
        # must be in the self.atom_types_str dict. Raise a RuntimeError if this
        # is not satisfied
        try:
            for key in nonbonded_types:
                self.atom_types_str[key].set_lj_params(*nonbonded_types[key])
        except KeyError:
            raise RuntimeError('Atom type %s not present in AtomType list' %
                               key)

        if parameterset is not None: self.parametersets.append(parameterset)
        if own_handle: f.close()

    def readTopologyFile(self, tfile):
        """Reads _only_ the atom type definitions from a topology file. This is
        unnecessary for versions 36 and later of the CHARMM force field.

        Parameters
        ----------
        tfile : str
            : Name of the CHARMM TOPology file to read

        Notes
        -----
        The CHARMM TOPology file is also called a Residue Topology File
        """
        conv = CharmmParameterSet._convert
        if isinstance(tfile, str):
            own_handle = True
            f = CharmmFile(tfile)
        else:
            own_handle = False
            f = tfile
        for line in f:
            line = line.strip()
            if line[:4] != 'MASS': continue
            words = line.split()
            try:
                idx = conv(words[1], int, 'atom type')
                name = words[2]
                mass = conv(words[3], float, 'atom mass')
            except IndexError:
                raise CharmmFileError('Could not parse MASS section of %s' %
                                      tfile)
            # The parameter file might or might not have an element name
            try:
                elem = words[4]
                atomic_number = get_by_symbol(elem).atomic_number
            except (IndexError, KeyError):
                # Figure it out from the mass
                masselem = Element.getByMass(mass)
                if masselem is None:
                    atomic_number = 0 # Extra point or something
                else:
                    atomic_number = masselem.atomic_number
            atype = AtomType(name=name, number=idx, mass=mass,
                             atomic_number=atomic_number)
            self.atom_types_str[atype.name] = atype
            self.atom_types_int[atype.number] = atype
            self.atom_types_tuple[(atype.name, atype.number)] = atype
        if own_handle: f.close()

    def readStreamFile(self, sfile):
        """Reads RTF and PAR sections from a stream file and dispatches the
        sections to readTopologyFile or readParameterFile

        Parameters
        ----------
        sfile : str or CharmmStreamFile
            Stream file to parse
        """
        if isinstance(sfile, CharmmStreamFile):
            f = sfile
        else:
            f = CharmmStreamFile(sfile)

        title, section = f.next_section()
        while title is not None and section is not None:
            words = title.lower().split()
            if words[1] == 'rtf':
                # This is a Residue Topology File section.
                self.readTopologyFile(section)
            elif words[1].startswith('para'):
                # This is a Parameter file section
                self.readParameterFile(section)
            title, section = f.next_section()

    def condense(self):
        """
        This function goes through each of the parameter type dicts and
        eliminates duplicate types. After calling this function, every unique
        bond, angle, dihedral, improper, or cmap type will pair with EVERY key
        in the type mapping dictionaries that points to the equivalent type

        Example
        -------
        >>> params = CharmmParameterSet('charmm.prm').condense()
        """
        # First scan through all of the bond types
        self._condense_types(self.bond_types)
        self._condense_types(self.angle_types)
        self._condense_types(self.urey_bradley_types)
        self._condense_types(self.improper_types)
        self._condense_types(self.cmap_types)
        # Dihedrals have to be handled separately, since each key corresponds to
        # a list of (potentially multiterm) dihedral terms. Since all terms in a
        # multiterm dihedral have to have a DIFFERENT periodicity, we don't have
        # to condense _within_ a single list of torsions assigned to the same
        # key (they're guaranteed to be different)
        keylist = list(self.dihedral_types.keys())
        for i in range(len(keylist) - 1):
            key1 = keylist[i]
            for dihedral in self.dihedral_types[key1]:
                for j in range(i+1, len(keylist)):
                    key2 = keylist[j]
                    for jj, dihedral2 in enumerate(self.dihedral_types[key2]):
                        if dihedral2 == dihedral:
                            self.dihedral_types[key2][jj] = dihedral
        return self

    @staticmethod
    def _condense_types(typedict):
        """
        Loops through the given dict and condenses all types.

        Parameters
        ----------
        typedict
            Type dictionary to condense
        """
        keylist = list(typedict.keys())
        for i in range(len(keylist) - 1):
            key1 = keylist[i]
            for j in range(i+1, len(keylist)):
                key2 = keylist[j]
                if typedict[key1] == typedict[key2]:
                    typedict[key2] = typedict[key1]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
