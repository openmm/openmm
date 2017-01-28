"""
desmonddmsfile.py: Load Desmond dms files

Portions copyright (c) 2013 Stanford University and the Authors
Authors: Robert McGibbon
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
from __future__ import absolute_import

import os
import math

from simtk import openmm as mm
from simtk.openmm.app import forcefield as ff
from simtk.openmm.app import Element, Topology, PDBFile
from simtk.openmm.app.element import hydrogen
from simtk.unit import (nanometer, angstrom, dalton, radian,
                        kilocalorie_per_mole, kilojoule_per_mole,
                        degree, elementary_charge)


class DesmondDMSFile(object):
    """DesmondDMSFile parses a Desmond DMS (desmond molecular system) and
    constructs a topology and (optionally) an OpenMM System from it
    """

    def __init__(self, file):
        """Load a DMS file

        Parameters
        ----------
        file : string
            the name of the file to load
        """
        # sqlite3 is included in the standard lib, but at python
        # compile time, you can disable support (I think), so it's
        # not *guarenteed* to be available. Doing the import here
        # means we only raise an ImportError if people try to use
        # this class, so the module can be safely imported
        import sqlite3

        self._open = False
        self._tables = None
        if not  os.path.exists(str(file)):
            raise IOError("No such file or directory: '%s'" % str(file))
        self._conn = sqlite3.connect(file)
        self._open = True
        self._readSchemas()

        if len(self._tables) == 0:
            raise IOError('DMS file was not loaded sucessfully. No tables found')
        if 'nbtype' not in self._tables['particle']:
            raise ValueError('No nonbonded parameters associated with this '
                             'DMS file. You can add a forcefield with the '
                             'viparr command line tool distributed with desmond')

        # build the provenance string
        provenance = []
        q = """SELECT id, user, timestamp, version, workdir, cmdline, executable
        FROM provenance"""
        #for id, user, timestamp, version, workdir, cmdline, executable in self._conn.execute(q):
        for row in self._conn.execute('SELECT * FROM provenance'):
            rowdict = dict(zip(self._tables['provenance'], row))
            provenance.append('%(id)d) %(timestamp)s: %(user)s\n  version: %(version)s\n  '
                              'cmdline: %(cmdline)s\n  executable: %(executable)s\n' % rowdict)
        self.provenance = ''.join(provenance)

        # Build the topology
        self.topology, self.positions = self._createTopology()
        self._topologyAtoms = list(self.topology.atoms())
        self._atomBonds = [{} for x in range(len(self._topologyAtoms))]
        self._angleConstraints = [{} for x in range(len(self._topologyAtoms))]

    def getPositions(self):
        """Get the positions of each atom in the system
        """
        return self.positions

    def getTopology(self):
        """Get the topology of the system
        """
        return self.topology

    def getProvenance(self):
        """Get the provenance string of this system
        """
        return self.provenance

    def _createTopology(self):
        """Build the topology of the system
        """
        top = Topology()
        positions = []

        boxVectors = []
        for x, y, z in self._conn.execute('SELECT x, y, z FROM global_cell'):
            boxVectors.append(mm.Vec3(x, y, z))
        unitCellDimensions = [boxVectors[0][0], boxVectors[1][1], boxVectors[2][2]]
        top.setUnitCellDimensions(unitCellDimensions*angstrom)

        atoms = {}
        lastChain = None
        lastResId = None
        c = top.addChain()
        q = """SELECT id, name, anum, resname, resid, chain, x, y, z
        FROM particle"""
        for (atomId, atomName, atomNumber, resName, resId, chain, x, y, z) in self._conn.execute(q):
            newChain = False
            if chain != lastChain:
                lastChain = chain
                c = top.addChain()
                newChain = True
            if resId != lastResId or newChain:
                lastResId = resId
                if resName in PDBFile._residueNameReplacements:
                    resName = PDBFile._residueNameReplacements[resName]
                r = top.addResidue(resName, c)
                if resName in PDBFile._atomNameReplacements:
                    atomReplacements = PDBFile._atomNameReplacements[resName]
                else:
                    atomReplacements = {}

            if atomNumber == 0 and atomName.startswith('Vrt'):
                elem = None
            else:
                elem = Element.getByAtomicNumber(atomNumber)

            if atomName in atomReplacements:
                atomName = atomReplacements[atomName]

            atoms[atomId] = top.addAtom(atomName, elem, r)
            positions.append(mm.Vec3(x, y, z))

        for p0, p1 in self._conn.execute('SELECT p0, p1 FROM bond'):
            top.addBond(atoms[p0], atoms[p1])

        positions = positions*angstrom
        return top, positions

    def createSystem(self, nonbondedMethod=ff.NoCutoff, nonbondedCutoff=1.0*nanometer,
                     ewaldErrorTolerance=0.0005, removeCMMotion=True, hydrogenMass=None):
        """Construct an OpenMM System representing the topology described by this
        DMS file

        Parameters
        ----------
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, or LJPME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        ewaldErrorTolerance : float=0.0005
            The error tolerance to use if nonbondedMethod is Ewald, PME, or LJPME.
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep their
            total mass the same.
        """
        self._checkForUnsupportedTerms()
        sys = mm.System()

        # Buld the box dimensions
        sys = mm.System()
        boxSize = self.topology.getUnitCellDimensions()
        if boxSize is not None:
            sys.setDefaultPeriodicBoxVectors((boxSize[0], 0, 0), (0, boxSize[1], 0), (0, 0, boxSize[2]))
        elif nonbondedMethod in (ff.CutoffPeriodic, ff.Ewald, ff.PME, ff.LJPME):
            raise ValueError('Illegal nonbonded method for a non-periodic system')

        # Create all of the particles
        for mass in self._conn.execute('SELECT mass from particle'):
            sys.addParticle(mass[0]*dalton)

        # Add all of the forces
        self._addBondsToSystem(sys)
        self._addAnglesToSystem(sys)
        self._addConstraintsToSystem(sys)
        self._addPeriodicTorsionsToSystem(sys)
        self._addImproperHarmonicTorsionsToSystem(sys)
        self._addCMAPToSystem(sys)
        self._addVirtualSitesToSystem(sys)
        nb = self._addNonbondedForceToSystem(sys)

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

        # Adjust masses.
        if hydrogenMass is not None:
            for atom1, atom2 in self.topology.bonds():
                if atom1.element == hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == hydrogen and atom1.element not in (hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Add a CMMotionRemover.
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        return sys

    def _addBondsToSystem(self, sys):
        """Create the harmonic bonds
        """
        bonds = mm.HarmonicBondForce()
        sys.addForce(bonds)

        q = """SELECT p0, p1, r0, fc, constrained
        FROM stretch_harm_term INNER JOIN stretch_harm_param
        ON stretch_harm_term.param=stretch_harm_param.id"""
        for p0, p1, r0, fc, constrained in self._conn.execute(q):
            if constrained:
                sys.addConstraint(p0, p1, r0*angstrom)
            else:
                # Desmond writes the harmonic bond force without 1/2
                # so we need to to double the force constant
                bonds.addBond(p0, p1, r0*angstrom, 2*fc*kilocalorie_per_mole/angstrom**2)

            # Record information that will be needed for constraining angles.
            self._atomBonds[p0][p1] = r0*angstrom
            self._atomBonds[p1][p0] = r0*angstrom

        return bonds

    def _addAnglesToSystem(self, sys):
        """Create the harmonic angles
        """
        angles = mm.HarmonicAngleForce()
        sys.addForce(angles)
        degToRad = math.pi/180

        q = """SELECT p0, p1, p2, theta0, fc, constrained
        FROM angle_harm_term INNER JOIN angle_harm_param
        ON angle_harm_term.param=angle_harm_param.id"""
        for p0, p1, p2, theta0, fc, constrained in self._conn.execute(q):
            if constrained:
                l1 = self._atomBonds[p1][p0]
                l2 = self._atomBonds[p1][p2]
                length = (l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta0*degToRad)).sqrt()
                sys.addConstraint(p0, p2, length)
                self._angleConstraints[p1][p0] = p2
                self._angleConstraints[p1][p2] = p0
            else:
                # Desmond writes the harmonic angle force without 1/2
                # so we need to to double the force constant
                angles.addAngle(p0, p1, p2, theta0*degToRad, 2*fc*kilocalorie_per_mole/radian**2)

        return angles

    def _addConstraintsToSystem(self, sys):
        """Add constraints to system. Normally these should already be
        added by the bonds table, but we want to make sure that there's
        no extra information in the constraints table that we're not
        including in the system"""
        for term_table in [n for n in self._tables.keys() if n.startswith('constraint_a') and n.endswith('term')]:
            param_table = term_table.replace('term', 'param')
            q = """SELECT p0, p1, r1
            FROM %(term)s INNER JOIN %(param)s
            ON %(term)s.param=%(param)s.id""" % \
                {'term': term_table, 'param': param_table}
            for p0, p1, r1 in self._conn.execute(q):
                if not p1 in self._atomBonds[p0]:
                    sys.addConstraint(p0, p1, r1*angstrom)
                    self._atomBonds[p0][p1] = r1*angstrom
                    self._atomBonds[p1][p0] = r1*angstrom

        if 'constraint_hoh_term' in self._tables:
            degToRad = math.pi/180
            q = """SELECT p0, p1, p2, r1, r2, theta
            FROM constraint_hoh_term INNER JOIN constraint_hoh_param
            ON constraint_hoh_term.param=constraint_hoh_param.id"""
            for p0, p1, p2, r1, r2, theta in self._conn.execute(q):
                # Here, p0 is the heavy atom and p1 and p2 are the H1 and H2
                # wihth O-H1 and O-H2 distances r1 and r2
                if not (self._angleConstraints[p0].get(p1, None) == p2):
                    length = (r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta*degToRad)).sqrt()
                    sys.addConstraint(p1, p2, length)

    def _addPeriodicTorsionsToSystem(self, sys):
        """Create the torsion terms
        """
        periodic = mm.PeriodicTorsionForce()
        sys.addForce(periodic)

        q = """SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
        FROM dihedral_trig_term INNER JOIN dihedral_trig_param
        ON dihedral_trig_term.param=dihedral_trig_param.id"""
        for p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 in self._conn.execute(q):
            for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                if fc == 0:
                    continue
                periodic.addTorsion(p0, p1, p2, p3, order, phi0*degree, fc*kilocalorie_per_mole)


    def _addImproperHarmonicTorsionsToSystem(self, sys):
        """Create the improper harmonic torsion terms
        """
        if not self._hasTable('improper_harm_term'):
            return

        harmonicTorsion = mm.CustomTorsionForce('k*(theta-theta0)^2')
        harmonicTorsion.addPerTorsionParameter('theta0')
        harmonicTorsion.addPerTorsionParameter('k')
        sys.addForce(harmonicTorsion)

        q = """SELECT p0, p1, p2, p3, phi0, fc
        FROM improper_harm_term INNER JOIN improper_harm_param
        ON improper_harm_term.param=improper_harm_param.id"""
        for p0, p1, p2, p3, phi0, fc in self._conn.execute(q):
            harmonicTorsion.addTorsion(p0, p1, p2, p3, [phi0*degree, fc*kilocalorie_per_mole])

    def _addCMAPToSystem(self, sys):
        """Create the CMAP terms
        """
        if not self._hasTable('torsiontorsion_cmap_term'):
            return

        # Create CMAP torsion terms
        cmap = mm.CMAPTorsionForce()
        sys.addForce(cmap)
        cmap_indices = {}

        for name in [k for k in self._tables.keys() if k.startswith('cmap')]:
            size2 = self._conn.execute('SELECT COUNT(*) FROM %s' % name).fetchone()[0]
            fsize = math.sqrt(size2)
            if fsize != int(fsize):
                raise ValueError('Non-square CMAPs are not supported')
            size = int(fsize)

            map = [0 for i in range(size2)]
            for phi, psi, energy in self._conn.execute("SELECT phi, psi, energy FROM %s" % name):
                i = int((phi % 360) / (360.0 / size))
                j = int((psi % 360) / (360.0 / size))
                map[i+size*j] = energy
            index = cmap.addMap(size, map*kilocalorie_per_mole)
            cmap_indices[name] = index

        q = """SELECT p0, p1, p2, p3, p4, p5, p6, p7, cmapid
        FROM torsiontorsion_cmap_term INNER JOIN torsiontorsion_cmap_param
        ON torsiontorsion_cmap_term.param=torsiontorsion_cmap_param.id"""
        for p0, p1, p2, p3, p4, p5, p6, p7, cmapid in self._conn.execute(q):
            cmap.addTorsion(cmap_indices[cmapid], p0, p1, p2, p3, p4, p5, p6, p7)

    def _addNonbondedForceToSystem(self, sys):
        """Create the nonbonded force
        """
        nb = mm.NonbondedForce()
        sys.addForce(nb)

        q = """SELECT charge, sigma, epsilon
        FROM particle INNER JOIN nonbonded_param
        ON particle.nbtype=nonbonded_param.id"""
        for charge, sigma, epsilon in self._conn.execute(q):
            nb.addParticle(charge, sigma*angstrom, epsilon*kilocalorie_per_mole)

        for p0, p1 in self._conn.execute('SELECT p0, p1 FROM exclusion'):
            nb.addException(p0, p1, 0.0, 1.0, 0.0)

        q = """SELECT p0, p1, aij, bij, qij
        FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param
        ON pair_12_6_es_term.param=pair_12_6_es_param.id;"""
        for p0, p1, a_ij, b_ij, q_ij in self._conn.execute(q):
            a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
            b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
            q_ij = q_ij*elementary_charge**2

            if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                new_epsilon = 0
                new_sigma = 1
            else:
                new_epsilon =  b_ij**2/(4*a_ij)
                new_sigma = (a_ij / b_ij)**(1.0/6.0)
            nb.addException(p0, p1, q_ij, new_sigma, new_epsilon, True)

        n_total = self._conn.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
        n_in_exclusions = self._conn.execute("""SELECT COUNT(*)
        FROM exclusion INNER JOIN pair_12_6_es_term
        ON exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1""").fetchone()
        if not n_total == n_in_exclusions:
            raise NotImplementedError('All pair_12_6_es_terms must have a corresponding exclusion')

        return nb

    def _addVirtualSitesToSystem(self, sys):
        """Create any virtual sites in the systempy
        """
        if not any(t.startswith('virtual_') for t in self._tables.keys()):
            return

        if 'virtual_lc2_term' in self._tables:
            q = """SELECT p0, p1, p2, c1
            FROM virtual_lc2_term INNER JOIN virtual_lc2_param
            ON virtual_lc2_term.param=virtual_lc2_param.id"""
            for p0, p1, p2, c1 in self._conn.execute(q):
                vsite = mm.TwoParticleAverageSite(p1, p2, (1-c1), c1)
                sys.setVirtualSite(p0, vsite)

        if 'virtual_lc3_term' in self._tables:
            q = """SELECT p0, p1, p2, p3, c1, c2
            FROM virtual_lc3_term INNER JOIN virtual_lc3_param
            ON virtual_lc3_term.param=virtual_lc3_param.id"""
            for p0, p1, p2, p3, c1, c2 in self._conn.execute(q):
                vsite = mm.ThreeParticleAverageSite(p1, p2, p3, (1-c1-c2), c1, c2)
                sys.setVirtualSite(p0, vsite)

        if 'virtual_out3_term' in self._tables:
            q = """SELECT p0, p1, p2, p3, c1, c2, c3
            FROM virtual_out3_term INNER JOIN virtual_out3_param
            ON virtual_out3_term.param=virtual_out3_param.id"""
            for p0, p1, p2, p3, c1, c2, c3 in self._conn.execute(q):
                vsite = mm.OutOfPlaneSite(p1, p2, p3, c1, c2, c3)
                sys.setVirtualSite(p0, vsite)

        if 'virtual_fdat3_term' in self._tables:
            raise NotImplementedError('OpenMM does not currently support '
                                      'fdat3-style virtual sites')


    def _hasTable(self, table_name):
        """Does our DMS file contain this table?
        """
        return table_name in self._tables

    def _readSchemas(self):
        """Read the schemas of each of the tables in the dms file, populating
        the `_tables` instance attribute
        """
        tables = {}
        for table in self._conn.execute("SELECT name FROM sqlite_master WHERE type='table'"):
            names = []
            for e in self._conn.execute('PRAGMA table_info(%s)' % table):
                names.append(str(e[1]))
            tables[str(table[0])] = names
        self._tables = tables

    def _checkForUnsupportedTerms(self):
        """Check the file for forcefield terms that are not currenty supported,
        raising a NotImplementedError
        """
        if 'posre_harm_term' in self._tables:
            raise NotImplementedError('Position restraints are not implemented.')
        flat_bottom_potential_terms = ['stretch_fbhw_term', 'angle_fbhw_term',
                                       'improper_fbhw_term', 'posre_fbhw_term']
        if any((t in self._tables) for t in flat_bottom_potential_terms):
            raise NotImplementedError('Flat bottom potential terms '
                                      'are not implemeneted')

        nbinfo = dict(zip(self._tables['nonbonded_info'],
                          self._conn.execute('SELECT * FROM nonbonded_info').fetchone()))

        if nbinfo['vdw_funct'] != u'vdw_12_6':
            raise NotImplementedError('Only Leonard-Jones van der Waals '
                                      'interactions are currently supported')
        if nbinfo['vdw_rule'] != u'arithmetic/geometric':
            raise NotImplementedError('Only Lorentz-Berthelot nonbonded '
                                      'combining rules are currently supported')

        if 'nonbonded_combined_param' in self._tables:
            raise NotImplementedError('nonbonded_combined_param interactions '
                                      'are not currently supported')

        if 'alchemical_particle' in self._tables:
            raise NotImplementedError('Alchemical particles are not supported')
        if 'alchemical_stretch_harm' in self._tables:
            raise NotImplementedError('Alchemical bonds are not supported')

        if 'polar_term' in self._tables:
            if self._conn.execute("SELECT COUNT(*) FROM polar_term").fetchone()[0] != 0:
                raise NotImplementedError('Drude particles are not currently supported')

    def close(self):
        """Close the SQL connection
        """
        if self._open:
            self._conn.close()

    def __del__(self):
        self.close()
