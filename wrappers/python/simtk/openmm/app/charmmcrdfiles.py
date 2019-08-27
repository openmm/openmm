"""
Provides a class for parsing CHARMM-style coordinate files, namely CHARMM .crd
(coordinate) files and CHARMM .rst (restart) file. Uses CharmmFile class in
_charmmfile.py for reading files

This file is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of Biological
Structures at Stanford, funded under the NIH Roadmap for Medical Research,
grant U54 GM072970. See https://simtk.org.  This code was originally part of
the ParmEd program and was ported for use with OpenMM.

Copyright (c) 2014 the Authors

Author: Jason Deckman
Contributors: Jason M. Swails
Date: June 6, 2014

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
from __future__ import print_function

from simtk.openmm.app.internal.charmm.exceptions import CharmmFileError
import simtk.unit as u
from simtk.openmm.vec3 import Vec3

CHARMMLEN = 22
TIMESCALE = 4.888821E-14 * 1e12 # AKMA time units to picoseconds
ONE_TIMESCALE = 1 / TIMESCALE

class CharmmCrdFile(object):
    """
    Reads and parses a CHARMM coordinate file (.crd) into its components,
    namely the coordinates, CHARMM atom types, resname, etc.

    Attributes
    ----------
    natom : int
        Number of atoms in the system
    resname : list
        Names of all residues
    positions : list
        All cartesian coordinates [x1, y1, z1, x2, ...]

    Examples
    --------
    >>> chm = CharmmCrdFile('testfiles/1tnm.crd')
    >>> print '%d atoms; %d coords' % (chm.natom, len(chm.positions))
    1414 atoms; 1414 coords
    """

    def __init__(self, fname):
        self.atomno = []                   # Atom number
        self.resno = []                    # Residue number
        self.resname = []                  # Residue name
        self.attype = []                   # Atom type
        self.positions = []                # 3N atomic coordinates
        self.title = []                    # .crd file title block
        self.segid = []                    # Segment ID
        self.weighting = []                # Atom weighting

        self.natom = 0                     # Number of atoms specified in file
        self._parse(fname)

    def _parse(self, fname):

        with open(fname, 'r') as crdfile:
            line = crdfile.readline()
    
            while len(line.strip()) == 0:      # Skip whitespace, as a precaution
                line = crdfile.readline()
    
            intitle = True
    
            while intitle:
                self.title.append(line.strip())
                line = crdfile.readline()
                if len(line.strip()) == 0:
                    intitle = False
                elif line.strip()[0] != '*':
                    intitle = False
                else:
                    intitle = True
    
            while len(line.strip()) == 0:      # Skip whitespace
                line = crdfile.readline()
    
            try:
                self.natom = int(line.strip().split()[0])
    
                for _ in range(self.natom):
                    line = crdfile.readline().strip().split()
                    self.atomno.append(int(line[0]))
                    self.resno.append(int(line[1]))
                    self.resname.append(line[2])
                    self.attype.append(line[3])
                    pos = Vec3(float(line[4]), float(line[5]), float(line[6]))
                    self.positions.append(pos)
                    self.segid.append(line[7])
                    self.weighting.append(float(line[9]))
    
                if self.natom != len(self.positions):
                    raise CharmmFileError("Error parsing CHARMM .crd file: %d "
                                          "atoms requires %d positions (not %d)" %
                                          (self.natom, self.natom,
                                           len(self.positions))
                    )
    
            except (ValueError, IndexError):
                raise CharmmFileError('Error parsing CHARMM coordinate file')

        # Apply units to the positions now. Do it this way to allow for
        # (possible) numpy functionality in the future.
        self.positions = u.Quantity(self.positions, u.angstroms)

class CharmmRstFile(object):
    """
    Reads and parses data, velocities and coordinates from a CHARMM restart
    file (.rst) of file name 'fname' into class attributes

    Attributes
    ----------
    natom : int
        Number of atoms in the system
    resname : list
        Names of all residues
    positions : list
        All cartesian coordinates [x1, y1, z1, x2, ...]
    positionsold : list
        Old cartesian coordinates
    velocities : list
        List of all cartesian velocities

    Examples
    --------
    >>> chm = CharmmRstFile('testfiles/sample-charmm.rst')
    >>> print chm.header[0]
    REST    37     1
    >>> natom, nc, nco = chm.natom, len(chm.positions), len(chm.positionsold)
    >>> nv = len(chm.velocities)
    >>> print '%d atoms; %d crds; %d old crds; %d vels' % (natom, nc, nco, nv)
    256 atoms; 256 crds; 256 old crds; 256 vels
    """

    def __init__(self, fname):
        self.header = []
        self.title = []
        self.enrgstat = []
        self.positionsold = []
        self.positions = []
        self.velocities = []

        self.ff_version = 0
        self.natom = 0
        self.npriv = 0
        self.nstep = 0
        self.nsavc = 0
        self.nsavv = 0
        self.jhstrt = 0

        self._parse(fname)

    def _parse(self, fname):

        crdfile = open(fname, 'r')
        readingHeader = True

        while readingHeader:
            line = crdfile.readline()
            if not len(line):
                raise CharmmFileError('Premature end of file')
            line = line.strip()
            words = line.split()
            if len(line) != 0:
                if words[0] == 'ENERGIES' or words[0] == '!ENERGIES':
                    readingHeader = False
                else:
                    self.header.append(line.strip())
            else:
                self.header.append(line.strip())

        for row in range(len(self.header)):
            if len(self.header[row].strip()) != 0:
                line = self.header[row].strip().split()
                if line[0][0:5] == 'NATOM' or line[0][0:6] == '!NATOM':
                    try:
                        line = self.header[row+1].strip().split()
                        self.natom = int(line[0])
                        self.npriv = int(line[1])     # num. previous steps
                        self.nstep = int(line[2])     # num. steps in file
                        self.nsavc = int(line[3])     # coord save frequency
                        self.nsavv = int(line[4])     # velocities "
                        self.jhstrt = int(line[5])    # Num total steps?
                        break

                    except (ValueError, IndexError) as e:
                        raise CharmmFileError('Problem parsing CHARMM restart')

        self._scan(crdfile, '!XOLD')
        self._get_formatted_crds(crdfile, self.positionsold)

        self._scan(crdfile, '!VX')
        self._get_formatted_crds(crdfile, self.velocities)

        self._scan(crdfile, '!X')
        self._get_formatted_crds(crdfile, self.positions)

        # Convert velocities to angstroms/ps
        self.velocities = [v * ONE_TIMESCALE for v in self.velocities]

        # Add units to positions and velocities
        self.positions = u.Quantity(self.positions, u.angstroms)
        self.positionsold = u.Quantity(self.positionsold, u.angstroms)
        self.velocities = u.Quantity(self.velocities, u.angstroms/u.picoseconds)

    def _scan(self, handle, str, r=0): # read lines in file until str is found
        scanning = True

        if(r): handle.seek(0)

        while scanning:
            line = handle.readline()
            if not line:
                raise CharmmFileError('Premature end of file')

            if len(line.strip()) != 0:
                if line.strip().split()[0][0:len(str)] == str:
                    scanning = False


    def _get_formatted_crds(self, crdfile, crds):
        for row in range(self.natom):
            line = crdfile.readline()

            if not line:
                raise CharmmFileError('Premature end of file')

            if len(line) < 3*CHARMMLEN:
                raise CharmmFileError("Less than 3 coordinates present in "
                                      "coordinate row or positions may be "
                                      "truncated.")

            line = line.replace('D','E')     # CHARMM uses 'D' for exponentials

            # CHARMM uses fixed format (len = CHARMMLEN = 22) for crds in .rst's

            c = Vec3(float(line[0:CHARMMLEN]), float(line[CHARMMLEN:2*CHARMMLEN]),
                     float(line[2*CHARMMLEN:3*CHARMMLEN]))
            crds.append(c)


    def printcoords(self, crds):
        for crd in range(len(crds)):
            print(crds[crd], end=' ')
            if not (crd+1) % 3:
                print('\n', end=' ')

if __name__ == '__main__':
    import doctest
    doctest.testmod()
