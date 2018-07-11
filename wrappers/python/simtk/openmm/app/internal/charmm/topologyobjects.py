"""
This module contains data structures useful for parsing in CHARMM files and
constructing chemical structures from those files

This file is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of Biological
Structures at Stanford, funded under the NIH Roadmap for Medical Research,
grant U54 GM072970. See https://simtk.org.  This code was originally part of
the ParmEd program and was ported for use with OpenMM.

Copyright (c) 2014 the Authors

Author: Jason M. Swails
Contributors:
Date: August 15, 2014

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
from simtk.openmm.app.internal.charmm.exceptions import (
                SplitResidueWarning, BondError, ResidueError, CmapError,
                MissingParameter)
import simtk.unit as u
import warnings

TINY = 1e-8

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _tracking(fcn):
    """ Decorator to indicate the list has changed """
    def new_fcn(self, *args):
        self.changed = True
        return fcn(self, *args)
    return new_fcn

class TrackedList(list):
    """
    This creates a list type that allows you to see if anything has changed
    """
    def __init__(self, arg=[]):
        self.changed = False
        list.__init__(self, arg)

    __delitem__ = _tracking(list.__delitem__)
    append = _tracking(list.append)
    extend = _tracking(list.extend)
    __setitem__ = _tracking(list.__setitem__)

# Python 3 does not have __delslice__, but make sure we override it for Python 2
if hasattr(TrackedList, '__delslice__'):
    TrackedList.__delslice__ = _tracking(TrackedList.__delslice__)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomType(object):
    """
    Atom types can either be compared by indexes or names. Can be assigned with
    a string, integer, (string is not automatically cast) or with both. Create
    new atom types with the "add" constructor to make sure the registry is
    filled with only unique types

    Parameters
    ----------
    name : str
        The name of the atom type
    number : int
        The integer index of the atom type
    mass : float
        The mass of the atom type
    atomic_number : int
        The atomic number of the element of the atom type

    Attributes
    ----------
    name : str
        The name of the atom type
    number : int
        The integer index of the atom type
    _member_number : int, private)
        The order in which this atom type was 'added' this is used to make
        sure that atom types added last have priority in assignment in the
        generated hash tables
    nbfix : dict
        Dictionary that maps nbfix terms with other atom types. Dict entries
        are (rmin, epsilon) -- precombined values for that particular atom pair
    nbthole : dict
            Dictionary that maps nbthole terms with other atom types. Dict entries
            are the value of a -- screening factor for that pair

    Examples
    --------
    >>> at = AtomType('HA', 1, 1.008, 1)
    >>> at.name, at.number
    ('HA', 1)
    >>> at2 = AtomType('CA', 2, 12.01, 6)
    >>> at2.name, at2.number
    ('CA', 2)
    >>> print "%s: %d" % (str(at), int(at))
    HA: 1
    >>> print at == WildCard
    True
    >>> print at2 == WildCard
    True
    """

    def __init__(self, name, number, mass, atomic_number):
        if number is None and name is not None:
            # If we were given an integer, assign it to number. Otherwise,
            # assign it to the name
            if isinstance(name, int):
                self.number = name
                self.name = None
            else:
                self.name = name
                self.number = None
        else:
            self.name = name
            self.number = int(number)
        self.mass = mass * u.daltons
        self.atomic_number = atomic_number
        # We have no LJ parameters as of yet
        self.epsilon = self.rmin = self.epsilon_14 = self.rmin_14 = None
        # Store each NBFIX term as a dict with the atom type string matching to
        # a 2-element tuple that is rmin, epsilon
        self.nbfix = dict()
        # Likewise, store each NBTHOLE term as a dict 
        self.nbthole = dict()

    def __eq__(self, other):
        """
        Compares based on available properties (name and number, just name,
        or just number)
        """
        if other is WildCard: return True # all atom types match wild cards
        if isinstance(other, AtomType):
            return self.name == other.name and self.number == other.number
        if isinstance(other, basestring):
            return self.name == other
        if isinstance(other, int):
            return self.number == other
        return other == (self.number, self.name)

    def set_lj_params(self, eps, rmin, eps14=None, rmin14=None):
        """ Sets Lennard-Jones parameters on this atom type """
        if eps14 is None:
            eps14 = eps
        if rmin14 is None:
            rmin14 = rmin
        self.epsilon = eps
        self.rmin = rmin
        self.epsilon_14 = eps14
        self.rmin_14 = rmin14

    def __int__(self):
        """ The integer representation of an AtomType is its index """
        return self.number

    def add_nbfix(self, typename, rmin, epsilon, rmin14, epsilon14):
        """ Adds a new NBFIX exclusion for this atom """
        if rmin14 is None: rmin14 = rmin
        if epsilon14 is None: epsilon14 = epsilon
        self.nbfix[typename] = (rmin, epsilon, rmin14, epsilon14)

    def add_nbthole(self, typename, nbt):
        """ Adds a new NBTHOLE screening factor for this atom """
        self.nbthole[typename] = (nbt)

    def __str__(self):
        return self.name

    # Comparisons are all based on number
    def __gt__(self, other):
        return self._member_number > other._member_number
    def __lt__(self, other):
        return self._member_number < other._member_number
    def __ge__(self, other):
        return self._member_number > other._member_number or self == other
    def __le__(self, other):
        return self._member_number < other._member_number or self == other

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class WildCard(AtomType):
    """
    This is a wild-card atom type that matches ANY atom type. It is often used
    in torsion parameters for the end-groups. Some properties:

        - It is a singleton, so seeing if an AtomType is a WildCard should be
          done using the "is" binary operator

        - It is defined as "equal" to every other object via the "==" operator

        - It is less than everything

        - It is greater than nothing
    """

    def __init__(self):
        self._member_number = -1
        self.name = 'X'
        self.number = None
        self.mass = None

    # Define comparison operators
    def __eq__(self, other): return True
    def __lt__(self, other): return True
    def __gt__(self, other): return False
    def __le__(self, other): return True
    def __ge__(self, other): return True

WildCard = WildCard() # Turn it into a singleton

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Atom(object):
    """ An atom in a structure.

    Parameters
    ----------
    system : str
        Name of the system this atom belongs to
    name : str
        name of the atom
    type : str or int
        Type of the atom
    charge : float
        Partial atomic charge (elementary charge units)
    mass : float
        Atomic mass (amu)
    props : list
        Other properties from the PSF

    Attributes
    ----------
    attype : str
        Name of the atom type
    system : str
        The system name associated with this atom
    name : str
        Name of the atom (str)
    charge : float
        Partial atomic charge
    mass : float
        Mass of the atom (amu)
    idx : int
        index of the atom in the system, starting from 0
    props : list
        List of extraneous properties parsed from a PSF
    type : AtomType
        If assigned, has additional properties like the non-bonded LJ
        parameters. If None, it has not yet been assigned

    Possible Attributes
    -------------------
    bond_partners : set of Atoms
        List of all atoms I am bonded to
    angle_partners set of Atoms
        List of all atoms I am angled to
    dihedral_partners : et of Atoms
        List of all atoms I am dihedraled to
    bonds : list of Bonds
        All bonds to which I belong
    angles : list of Angles
        All angles to which I belong
    dihedrals : list of Dihedrals
        All dihedrals to which I belong
    impropers : list of Impropers
        All impropers to which I belong
    cmaps : list of Cmaps
        All correction maps to which I belong
    """
    def __init__(self, system, name, attype, charge, mass, props=None):
        self.name = name
        self.attype = attype
        self.type = None
        self.charge = charge
        self.mass = mass * u.daltons
        self.idx = -1
        self.props = props
        self.system = system
        self.marked = 0 # For recursive molecule determination
        self._bond_partners = set()
        self._angle_partners = set()
        self._dihedral_partners = set()
        self.bonds = []
        self.angles = []
        self.urey_bradleys = []
        self.dihedrals = []
        self.impropers = []
        self.cmaps = []

    def bond_to(self, other):
        """
        Register this atom as bonded partners. Cannot bond to itself. If that
        is attempted, a BondError is raised
        """
        if self is other:
            raise BondError('Cannot bond atom to itself')
        self._bond_partners.add(other)
        other._bond_partners.add(self)

    def angle_to(self, other):
        """
        Register this atom as angle partners. Cannot angle to itself. If that
        is attempted, a BondError is raised
        """
        if self is other:
            raise BondError('Cannot angle atom to itself')
        self._angle_partners.add(other)
        other._angle_partners.add(self)

    def dihedral_to(self, other):
        """
        Register this atom as dihedral partners. Cannot dihedral to itself. If
        that is attempted, a BondError is raised
        """
        if self is other:
            raise BondError('Cannot dihedral atom to itself')
        self._dihedral_partners.add(other)
        other._dihedral_partners.add(self)

    @property
    def bond_partners(self):
        return sorted(list(self._bond_partners))

    @property
    def angle_partners(self):
        return sorted(list(self._angle_partners - self._bond_partners))

    @property
    def dihedral_partners(self):
        return sorted(list(self._dihedral_partners - self._angle_partners -
                           self._bond_partners))

    def type_to_int(self):
        """
        Changes the type to an integer, matching CHARMM conventions. This can
        only be done if a type mapping has been loaded (i.e., if the type
        attribute is not None).

        If the type identifier is not already an integer and no mapping is
        available, MissingParameter is raised.
        """
        if isinstance(self.attype, int):
            return
        if self.type is None:
            raise MissingParameter('No type mapping loaded. Cannot change '
                                   'type identifier to integer for %s' %
                                   self.attype)
        self.attype = self.type.number

    def type_to_str(self):
        """
        Changes the type to a string, matching XPLOR conventions. This can only
        be done if a type mapping has been loaded (i.e., if the type attribute
        is not None).

        If the type identifier is not already a string and no mapping is
        available, MissingParameter is raised.
        """
        if isinstance(self.attype, str):
            return
        if self.type is None:
            raise MissingParameter('No type mapping loaded. Cannot change '
                                   'type identifier to string for %s' %
                                   self.attype)
        self.attype = self.type.name

    def __repr__(self):
        retstr = '<Atom %d' % self.idx
        if hasattr(self, 'residue'):
            retstr += '; %d %s [%s: %s]>' % (
                        self.residue.idx, self.residue.resname, self.name,
                        self.attype)
        else:
            retstr += '; %s> ' % (self.name)
        return retstr

    def __lt__(self, other):
        return self.idx < other.idx
    def __gt__(self, other):
        return self.idx > other.idx
    def __le__(self, other):
        return not self > other
    def __ge__(self, other):
        return not self < other

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomList(TrackedList):
    """ A list of Atom instances.  """

    def unmark(self):
        for atom in self: atom.marked = 0

    def assign_indexes(self):
        self._index_us()

    def _index_us(self):
        for i, atom in enumerate(self): atom.idx = i

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Residue(object):
    """ Residue class """

    def __init__(self, resname, idx, inscode=''):
        self.resname = resname
        self.idx = idx
        self.atoms = []
        self.resnum = None # Numbered based on SYSTEM name
        self.system = None
        self.inscode = inscode

    def add_atom(self, atom):
        if self.system is None:
            self.system = atom.system
        else:
            if self.system != atom.system:
                raise ResidueError('Added atom has a different system than '
                                   'the other atoms in the residue!')
        atom.residue = self
        self.atoms.append(atom)

    def delete_atom(self, atom):
        """
        If an atom is present in this residue, delete it from the list of
        atoms
        """
        self.atoms = [a for a in self.atoms if a is not atom]

    # Implement some container methods over the list of atoms
    def __contains__(self, thing):
        """ True if an atom is present in this residue """
        return thing in self.atoms

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, idx):
        return self.atoms[idx]

    # Equality
    def __eq__(self, thing):
        if isinstance(thing, Residue):
            return self.resname == thing.resname and self.idx == thing.idx
        if isinstance(thing, tuple) or isinstance(thing, list):
            # Must be resnum, resname
            return thing == (self.resname, self.idx, self.inscode)
        return False # No other type can be equal.

    def __ne__(self, thing):
        return not self.__eq__(thing)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueList(list):
    """
    A list of residues in a chemical system. Atoms are added to this list and
    spawn new residues when the residue number or name changes.
    """
    def __init__(self):
        list.__init__(self)
        self._last_residue = None

    def assign_indexes(self):
        """ Indexes all residues starting from 1 """
        last_system = None
        resnum = 1
        for i, res in enumerate(self):
            res.idx = i + 1
            if res.system != last_system:
                res.resnum = resnum = 1
                last_system = res.system
            else:
                res.resnum = resnum
            resnum += 1

    def add_atom(self, system, resnum, resname, name,
                 attype, charge, mass, inscode, props=None):
        """Adds an atom to the list of residues. If the residue is not the same as
        the last residue that was created, a new Residue is created and added to
        this list

        Parameters
        ----------
        system : str
            The system this atom belongs to
        resnum : int
            Residue number
        resname : str
            Name of the residue
        name : str
            Name of the atom
        attype : int or str
            Type of the atom
        charge : float
            Partial atomic charge of the atom
        mass : float
            Mass (amu) of the atom
        inscode : str
            Insertion code, if it is specified
        props : list
            Other properties from the PSF

        Returns
        -------
            The Atom instance created and added to the list of residues
        """
        lr = self._last_residue
        if lr is None:
            res = self._last_residue = Residue(resname, resnum, inscode)
            list.append(self, res)
        elif (lr.resname != resname or lr.idx != resname or
                lr.inscode != inscode  or system != lr.system):
            res = self._last_residue = Residue(resname, resnum, inscode)
            res.system = system
            list.append(self, res)
        else:
            res = lr
        atom = Atom(system, name, attype, float(charge), float(mass), props)
        res.add_atom(atom)
        return atom

    def append(self, thing):
        raise NotImplementedError('Use "add_atom" to build a residue list')

    extend = append

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Bond(object):
    """A bond object that links 2 atoms

    Parameters
    ----------
    atom1 : Atom
        First atom included in the bond
    atom2 : Atom
        Second atom included in the bond
    bond_type : BondType
        Type for the bond (None if unknown)
    """
    def __init__(self, atom1, atom2, bond_type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_type = bond_type
        self.atom1.bond_to(atom2)
        # Add this bond to the atoms
        self.atom1.bonds.append(self)
        self.atom2.bonds.append(self)

    def __repr__(self):
        return '<Bond; %r -- %r; type=%r>' % (self.atom1, self.atom2,
                                              self.bond_type)

    def __contains__(self, thing):
        """ See if an atom is part of this bond """
        return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Angle(object):
    """An angle object that links 3 atoms

    Parameters
    ----------
    atom1 : Atom
        First atom included in the angle
    atom2 : Atom
        Central atom in the valence angle
    atom3 : Atom
        Third atom in the valence angle
    angle_type : AngleType
        Type for the angle (None if unknown)
    """
    def __init__(self, atom1, atom2, atom3, angle_type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.angle_type = angle_type
        self.atom1.angle_to(atom3)
        # Add this angle to the atoms
        self.atom1.angles.append(self)
        self.atom2.angles.append(self)
        self.atom3.angles.append(self)

    def __repr__(self):
        return '<Angle; %r-%r-%r; type=%r>' % (self.atom1, self.atom2,
                                               self.atom3, self.angle_type)

    def __contains__(self, thing):
        """ See if a bond or an atom is in this angle """
        if isinstance(thing, Bond):
            return self.atom2 in thing and (self.atom1 in thing or
                                            self.atom3 in thing)
        # Otherwise assume it's an atom
        return self.atom1 is thing or self.atom2 is thing or self.atom3 is thing

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradley(object):
    """A harmonic restraint between two atoms separated by 2 valence bonds
    (i.e., involved in a valence angle with each other

    Parameters
    ----------
    atom1 : Atom
        The first atom included in the Urey-Bradley term
    atom2 : Atom
        The second atom included in the Urey-Bradley term
    ub_type : UreyBradleyType
        The type for the Urey-Bradley term (None if unknown)
    """
    def __init__(self, atom1, atom2, ub_type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.ub_type = ub_type
        # Add this Urey-Bradley to the atoms
        atom1.urey_bradleys.append(self)
        atom2.urey_bradleys.append(self)

    def __repr__(self):
        return '<UreyBradley; %r-%r; type=%r>' % (self.atom1, self.atom2,
                                                  self.ub_type)

    def __contains__(self, thing):
        # See the comments in chemistry/amber/topologyobjects.py for what's
        # being done here (under the same method of the UreyBradley type there)
        if isinstance(thing, Bond):
            if not thing.atom1 in self:
                if not thing.atom2 in self:
                    return False
                end1 = thing.atom2
                cent = thing.atom1
            else:
                end1 = thing.atom1
                cent = thing.atom2
            for atm in cent.bonds():
                if atm is end1: continue
                if atm in self:
                    return True
            return False
        # thing must be an atom
        return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Dihedral(object):
    """A torsion angle object that links 4 atoms

    Parameters
    ----------
    atom1 : Atom
        First atom included in the torsion
    atom2 : Atom
        Second atom included in the torsion
    atom3 : Atom
        Third atom included in the torsion
    atom4 : Atom
        Fourth atom included in the torsion
    dihedral_type : DihedralType
        Type for the torsion (None if unknown)
    """
    def __init__(self, atom1, atom2, atom3, atom4, dihedral_type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.dihedral_type = dihedral_type
        self.atom1.dihedral_to(atom4)
        # Add this torsion to the atoms
        self.atom1.dihedrals.append(self)
        self.atom2.dihedrals.append(self)
        self.atom3.dihedrals.append(self)
        self.atom4.dihedrals.append(self)
        # Add a marker for indicating if this dihedral is not the final term in
        # a multi-term expansion or if the atoms at the end are also 1-2 or 1-3
        # pairs (this can happen for 5- and 6-member rings, respectively)
        self.end_groups_active = True

    def __repr__(self):
        return '<Dihedral; %r-%r-%r-%r; type=%r>' % (
                    self.atom1, self.atom2, self.atom3, self.atom4,
                    self.dihedral_type)

    def __contains__(self, thing):
        """ See if a bond or an atom is in this torsion """
        if isinstance(thing, Bond):
            if self.atom1 in thing and self.atom2 in thing:
                return True
            if self.atom2 in thing and self.atom3 in thing:
                return True
            if self.atom3 in thing and self.atom4 in thing:
                return True
            return False
        # Otherwise assume it's an atom
        return (thing is self.atom1 or thing is self.atom2 or
                thing is self.atom3 or thing is self.atom4)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Improper(object):
    """An improper torsion object. The third atom is bonded to each other atom

    Parameters
    ----------
    atom1 : Atom
        First atom included in the torsion
    atom2 : Atom
        Second atom included in the torsion
    atom3 : Atom
        Third atom included in the torsion
    atom4 : Atom
        Fourth atom included in the torsion
    improper_type : ImproperType
        Type for the improper (None if unknown)
    """
    def __init__(self, atom1, atom2, atom3, atom4, improper_type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.improper_type = improper_type
        # Nothing to register -- all bonds should already be formed
        # Add this improper to the atoms
        self.atom1.impropers.append(self)
        self.atom2.impropers.append(self)
        self.atom3.impropers.append(self)
        self.atom4.impropers.append(self)

    def __repr__(self):
        return '<Improper; %r-%r-%r-%r; type=%r>' % (
                    self.atom1, self.atom2, self.atom3, self.atom4,
                    self.improper_type)

    def __contains__(self, thing):
        """
        See if a bond or an atom is in this improper
        An improper is defined as shown below

                                A3
                                |
                                |
                       A4 ----- A1 ----- A2

        So the bonds will either be between atom1 and any other atom
        """
        if isinstance(thing, Bond):
            if self.atom2 in thing and self.atom1 in thing:
                return True
            if self.atom3 in thing and self.atom1 in thing:
                return True
            if self.atom4 in thing and self.atom1 in thing:
                return True
            return False
        # Otherwise assume it's an atom
        return (thing is self.atom1 or thing is self.atom2 or
                thing is self.atom3 or thing is self.atom4)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AcceptorDonor(object):
    """Just a holder for donors and acceptors in CHARMM speak

    Parameters
    ----------
    atom1 : Atom
        First atom in the donor/acceptor group
    atom2 : Atom
        Second atom in the donor/acceptor group
    """
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2

    def __repr__(self):
        return '<AcceptorDonor; %r %r>' % (self.atom1, self.atom2)

    def __contains__(self, thing):
        """ See if the atom is in this donor/acceptor """
        return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Group(object):
    """An 'interacting' group defined by the PSF.

    Parameters
    ----------
    bs : int
        ??
    type : int
        The group type
    move : int
        If the group moves ??

    Disclaimer: I really don't know what these numbers mean. I'm speculating
    based on the source code of 'chamber', and this section is simply ignored
    there.
    """
    def __init__(self, bs, type, move):
        self.bs = bs
        self.type = type
        self.move = move

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Cmap(object):
    """
    A correction-map (two coupled torsions). Defined between 8 atoms (2
    consecutive correction maps). "Consecutive torsions" (i.e., those definable
    by 5 atoms) will only be recognized if the two torsions have the same order

    Parameters
    ----------
    atom1 : Atom
        1st atom of first dihedral
    atom2 : Atom
        2nd atom of first dihedral
    atom3 : Atom
        3rd atom of first dihedral
    atom4 : Atom
        4th atom of first dihedral
    atom5 : Atom
        1st atom of second dihedral
    atom6 : Atom
        2nd atom of second dihedral
    atom7 : Atom
        3rd atom of second dihedral
    atom8 : Atom
        4th atom of second dihedral
    cmap_type : CmapType
        Cmap type for this cmap (None if unknown)

    Attributes
    ----------
    consecutive : bool
        Are the dihedrals consecutive?

    if consecutive:
        atom1 : Atom
            1st atom of 1st dihedral
        atom2 L Atom
            2nd atom of 1st dihedral && 1st atom of 2nd dihedral
        atom3 : Atom
            3rd atom of 1st dihedral && 2nd atom of 2nd dihedral
        atom4 : Atom
            4th atom of 1st dihedral && 3rd atom of 2nd dihedral
        atom5 : Atom
            4th atom of 2nd dihedral

     if not consecutive:
         atom1 : Atom
            1st atom of first dihedral
         atom2 : Atom
            2nd atom of first dihedral
         atom3 : Atom
            3rd atom of first dihedral
         atom4 : Atom
            4th atom of first dihedral
         atom5 : Atom
            1st atom of second dihedral
         atom6 : Atom
            2nd atom of second dihedral
         atom7 : Atom
            3rd atom of second dihedral
         atom8 : Atom
            4th atom of second dihedral
    """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, atom6, atom7,
                 atom8, cmap_type=None):
        self.consecutive = False
        if atom2 is atom5 and atom3 is atom6 and atom4 is atom7:
            self.consecutive = True
        if self.consecutive:
            self.atom1 = atom1
            self.atom2 = atom2
            self.atom3 = atom3
            self.atom4 = atom4
            self.atom5 = atom8
            # Add this cmap to the atoms
            self.atom1.cmaps.append(self)
            self.atom2.cmaps.append(self)
            self.atom3.cmaps.append(self)
            self.atom4.cmaps.append(self)
            self.atom5.cmaps.append(self)
        else:
            self.atom1 = atom1
            self.atom2 = atom2
            self.atom3 = atom3
            self.atom4 = atom4
            self.atom5 = atom5
            self.atom6 = atom6
            self.atom7 = atom7
            self.atom8 = atom8
            # Add this cmap to the atoms
            self.atom1.cmaps.append(self)
            self.atom2.cmaps.append(self)
            self.atom3.cmaps.append(self)
            self.atom4.cmaps.append(self)
            self.atom5.cmaps.append(self)
            self.atom6.cmaps.append(self)
            self.atom7.cmaps.append(self)
            self.atom8.cmaps.append(self)
        self.cmap_type = cmap_type

    def __repr__(self):
        if self.consecutive:
            return '<Cmap; %r-%r-%r-%r-%r; cmap_type=%r>' % (
                        self.atom1, self.atom2, self.atom3, self.atom4,
                        self.atom5, self.cmap_type)
        else:
            return '<Cmap; %r-%r-%r-%r && %r-%r-%r-%r; cmap_type=%r>' % (
                        self.atom1, self.atom2, self.atom3, self.atom4,
                        self.atom5, self.atom6, self.atom7, self.atom8,
                        self.cmap_type)

    def __contains__(self, thing):
        """ See if a Bond or Atom is inside this torsion """
        if isinstance(thing, Bond):
            if self.consecutive:
                return self._consecutive_has_bond(thing)
            else:
                return self._nonconsecutive_has_bond(thing)
        # Must be an atom
        if self.consecutive:
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3 or thing is self.atom4 or
                    thing is self.atom5)
        return (thing is self.atom1 or thing is self.atom2 or
                thing is self.atom3 or thing is self.atom4 or
                thing is self.atom5 or thing is self.atom6 or
                thing is self.atom7 or thing is self.atom8)

    def _consecutive_has_bond(self, thing):
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom3 in thing and self.atom4 in thing) or
                (self.atom4 in thing and self.atom5 in thing))

    def _nonconsecutive_has_bond(self, thing):
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom3 in thing and self.atom4 in thing) or
                (self.atom5 in thing and self.atom6 in thing) or
                (self.atom6 in thing and self.atom7 in thing) or
                (self.atom7 in thing and self.atom8 in thing))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondType(object):
    """
    A bond type with an equilibrium length (Angstroms) and force constant
    (kcal/mol/Angstrom^2)

    Parameters
    ----------
    k : float
        : Force constant (kcal/mol/A^2)
    req : float
        : Equilibrium distance
    """
    def __init__(self, k, req):
        self.k = k
        self.req = req

    def __eq__(self, other):
        return self.k == other.k and self.req == other.req

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleType(object):
    """An angle type with an equilibrium angle (degrees) and force constant
    (kcal/mol/radians^2)

    Parameters
    ----------
    k : float
        Force constant (kcal/mol/radians^2)
    theteq : float
        Equilibrium angle value (degrees)
    """
    def __init__(self, k, theteq):
        self.k = k
        self.theteq = theteq

    def __eq__(self, other):
        return self.k == other.k and self.theteq == other.theteq

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradleyType(BondType):
    """
    A harmonic spring connecting two atoms separated by two valence bonds (a
    valence angle). It is functionally equivalent to a Bond and is actually
    implemented as a (unaltered) Bond subclass. See BondType documentation.
    """

# Not all angles have Urey-Bradley terms attached to them. This is a singleton
# that indicates that there is NO U-B term for this particular type
NoUreyBradley = UreyBradleyType(None, None)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralType(object):
    """A torsion angle type with a force constant (kcal/mol), periodicity
    (int), and phase (degrees)

    Parameters
    ----------
    phi_k : float
        Force constant (kcal/mol)
    per : int
        Periodicity
    phase : float
        Phase of the torsion
    """
    def __init__(self, phi_k, per, phase):
        self.phi_k = float(phi_k)
        self.per = int(per)
        self.phase = float(phase)

    def __repr__(self):
        return "<DihedralType: k=%s; phase=%s; per=%s>" % (self.phi_k,
                self.phase, self.per)

    def __eq__(self, other):
        return (self.phi_k == other.phi_k and self.per == other.per and
                self.phase == other.phase)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ImproperType(object):
    """An improper torsion angle type with a force constant (kcal/mol) and
    equilibrium angle (degrees)

    Parameters
    ----------
    k : float
        : Force constant (kcal/mol)
    phieq : int
        : Equilibrium angle (degrees)
    """
    def __init__(self, k, phieq):
        self.k = k
        self.phieq = phieq

    def __eq__(self, other):
        return self.k == other.k and self.phieq == other.phieq

    def __repr__(self):
        return '<ImproperType; k=%s; phieq=%s>' % (self.k, self.phieq)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CmapType(object):
    """
    Contains a correction map interpolation grid

    Parameters
    ----------
    resolution : int
        Number of interpolation points for each dihedral
    grid : list of floats
        resolution x resolution list of energy values (kcal/mol) for the
        angles with the 2nd angle changing fastest.

    The grid object is converted to a _CmapGrid instance which can be treated
    like a normal list, but also has the ability to quickly return a transpose
    (so the 1st angle changes fastest) and to switch the range from -180 -- 180
    to 0 -- 360 (and back again). This is particularly helpful because CHARMM
    defines CMAP tables from -180 -- 180 whereas OpenMM expects them from
    0 -- 360 (with the 1st angle changing fastest!!)
    """

    def __init__(self, resolution, grid):
        self.resolution = resolution
        self.grid = _CmapGrid(resolution, grid)
        if len(grid) != self.resolution * self.resolution:
            raise CmapError('CMAP grid does not match expected resolution')

    def __eq__(self, other):
        return (self.resolution == other.resolution and
                all([abs(i - j) < TINY for i, j in zip(self.grid, other.grid)]))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Take the CmapGrid class from the Chamber prmtop topology objects
class _CmapGrid(object):
    """
    A grid object for storing Correction map data. Data can be accessed in one
    of two ways; either with 1 or 2 indexes. If 2 indexes are given, the index
    into the flattened array is i*resolution+j. Indexing starts from 0.

    The _CmapGrid usually has ranges for the two angles from -180 to 180. Some
    places will expect the range to be 0-360 degrees (e.g., OpenMM). The
    switch_range method returns a _CmapGrid with this range. This method will
    also work to go backwards.

    Example:
    >>> g = _CmapGrid(2, [0, 1, 2, 3])
    >>> print g[0], g[0,0]
    0 0
    >>> print g[1], g[0,1]
    1 1
    >>> print g[1,0]
    2
    >>> g[1,1] = 10
    >>> print g[3]
    10
    >>> print g.switch_range()
    [10.0000, 2.0000
     1.0000, 0.0000]
    """

    def __init__(self, resolution, data=None):
        self.resolution = resolution
        if data is None:
            self._data = [0 for i in range(self.resolution*self.resolution)]
        else:
            self._data = data

    @property
    def transpose(self):
        """ Returns the transpose of the grid """
        try:
            return self._transpose
        except AttributeError:
            pass
        _transpose = []
        size = len(self._data)
        for i in range(self.resolution):
            piece = [self[j] for j in range(i, size, self.resolution)]
            _transpose += piece
        self._transpose = _CmapGrid(self.resolution, _transpose)
        return self._transpose

    T = transpose

    def __getitem__(self, idx):
        if isinstance(idx, tuple):
            return self._data[self.resolution*idx[0]+idx[1]]
        return self._data[idx]

    def __setitem__(self, idx, val):
        if isinstance(idx, tuple):
            self._data[self.resolution*idx[0]+idx[1]] = val
        else:
            self._data[idx] = val

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self):
        return '<_CmapGrid: %dx%d>' % (self.resolution, self.resolution)

    def __str__(self):
        retstr = '[%.4f,' % self._data[0]
        fmt = ' %.4f'
        for i, val in enumerate(self):
            if i == 0: continue
            retstr += fmt % val
            if (i+1) % self.resolution == 0 and i != len(self._data) - 1:
                retstr += '\n'
            elif i != len(self) - 1:
                retstr += ','
        return retstr + ']'

    def __eq__(self, other):
        try:
            if self.resolution != other.resolution:
                return False
            for x, y in zip(self, other):
                if abs(x - y) > TINY:
                    return False
            return True
        except AttributeError:
            return TypeError('Bad type comparison with _CmapGrid')

    def switch_range(self):
        """
        Returns a grid object whose range is 0 to 360 degrees in both dimensions
        instead of -180 to 180 degrees (or -180 to 180 degrees if the range is
        already 0 to 360 degrees)
        """
        res = self.resolution
        mid = res // 2
        newgrid = _CmapGrid(res)
        for i in range(res):
            for j in range(res):
                # Start from the middle
                newgrid[i, j] = self[(i+mid)%res, (j+mid)%res]
        return newgrid

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == '__main__':
    import doctest
    doctest.testmod()
