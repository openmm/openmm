"""
dcdfile.py: Used for writing DCD files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman
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
__author__ = "Peter Eastman"
__version__ = "1.0"

import array
import os
import time
import struct
import math
from simtk.unit import picoseconds, nanometers, angstroms, is_quantity, norm
from simtk.openmm import Vec3
from simtk.openmm.app.internal.unitcell import computeLengthsAndAngles

class DCDFile(object):
    """DCDFile provides methods for creating DCD files.

    DCD is a file format for storing simulation trajectories.  It is supported by many programs, such
    as CHARMM, NAMD, and X-PLOR.  Note, however, that different programs produce subtly different
    versions of the format.  This class generates the CHARMM version.  Also note that there is no
    standard byte ordering (big-endian or little-endian) for this format.  This class always generates
    files with little-endian ordering.

    To use this class, create a DCDFile object, then call writeModel() once for each model in the file."""

    def __init__(self, file, topology, dt, firstStep=0, interval=1, append=False):
        """Create a DCD file and write out the header, or open an existing file to append.

        Parameters
        ----------
        file : file
            A file to write to
        topology : Topology
            The Topology defining the molecular system being written
        dt : time
            The time step used in the trajectory
        firstStep : int=0
            The index of the first step in the trajectory
        interval : int=1
            The frequency (measured in time steps) at which states are written
            to the trajectory
        append : bool=False
            If True, open an existing DCD file to append to.  If False, create a new file.
        """
        self._file = file
        self._topology = topology
        self._firstStep = firstStep
        self._interval = interval
        self._modelCount = 0
        if is_quantity(dt):
            dt = dt.value_in_unit(picoseconds)
        dt /= 0.04888821
        self._dt = dt
        boxFlag = 0
        if topology.getUnitCellDimensions() is not None:
            boxFlag = 1
        if append:
            file.seek(8, os.SEEK_SET)
            self._modelCount = struct.unpack('<i', file.read(4))[0]
            file.seek(268, os.SEEK_SET)
            numAtoms = struct.unpack('<i', file.read(4))[0]
            if numAtoms != len(list(topology.atoms())):
                raise ValueError('Cannot append to a DCD file that contains a different number of atoms')
        else:
            header = struct.pack('<i4c9if', 84, b'C', b'O', b'R', b'D', 0, firstStep, interval, 0, 0, 0, 0, 0, 0, dt)
            header += struct.pack('<13i', boxFlag, 0, 0, 0, 0, 0, 0, 0, 0, 24, 84, 164, 2)
            header += struct.pack('<80s', b'Created by OpenMM')
            header += struct.pack('<80s', b'Created '+time.asctime(time.localtime(time.time())).encode('ascii'))
            header += struct.pack('<4i', 164, 4, len(list(topology.atoms())), 4)
            file.write(header)

    def writeModel(self, positions, unitCellDimensions=None, periodicBoxVectors=None):
        """Write out a model to the DCD file.

        The periodic box can be specified either by the unit cell dimensions
        (for a rectangular box), or the full set of box vectors (for an
        arbitrary triclinic box).  If neither is specified, the box vectors
        specified in the Topology will be used. Regardless of the value
        specified, no dimensions will be written if the Topology does not
        represent a periodic system.

        Parameters
        ----------
        positions : list
            The list of atomic positions to write
        unitCellDimensions : Vec3=None
            The dimensions of the crystallographic unit cell.
        periodicBoxVectors : tuple of Vec3=None
            The vectors defining the periodic box.
        """
        if len(list(self._topology.atoms())) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
        if is_quantity(positions):
            positions = positions.value_in_unit(nanometers)
        if any(math.isnan(norm(pos)) for pos in positions):
            raise ValueError('Particle position is NaN')
        if any(math.isinf(norm(pos)) for pos in positions):
            raise ValueError('Particle position is infinite')
        file = self._file

        self._modelCount += 1
        if self._interval > 1 and self._firstStep+self._modelCount*self._interval > 1<<31:
            # This will exceed the range of a 32 bit integer.  To avoid crashing or producing a corrupt file,
            # update the header to say the trajectory consisted of a smaller number of larger steps (so the
            # total trajectory length remains correct).
            self._firstStep //= self._interval
            self._dt *= self._interval
            self._interval = 1
            file.seek(0, os.SEEK_SET)
            file.write(struct.pack('<i4c9if', 84, b'C', b'O', b'R', b'D', 0, self._firstStep, self._interval, 0, 0, 0, 0, 0, 0, self._dt))

        # Update the header.

        file.seek(8, os.SEEK_SET)
        file.write(struct.pack('<i', self._modelCount))
        file.seek(20, os.SEEK_SET)
        file.write(struct.pack('<i', self._firstStep+self._modelCount*self._interval))

        # Write the data.

        file.seek(0, os.SEEK_END)
        boxVectors = self._topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            if periodicBoxVectors is not None:
                boxVectors = periodicBoxVectors
            elif unitCellDimensions is not None:
                if is_quantity(unitCellDimensions):
                    unitCellDimensions = unitCellDimensions.value_in_unit(nanometers)
                boxVectors = (Vec3(unitCellDimensions[0], 0, 0), Vec3(0, unitCellDimensions[1], 0), Vec3(0, 0, unitCellDimensions[2]))*nanometers
            (a_length, b_length, c_length, alpha, beta, gamma) = computeLengthsAndAngles(boxVectors)
            a_length = a_length * 10.  # computeLengthsAndAngles returns unitless nanometers, but need angstroms here.
            b_length = b_length * 10.  # computeLengthsAndAngles returns unitless nanometers, but need angstroms here.
            c_length = c_length * 10.  # computeLengthsAndAngles returns unitless nanometers, but need angstroms here.
            angle1 = math.sin(math.pi/2-gamma)
            angle2 = math.sin(math.pi/2-beta)
            angle3 = math.sin(math.pi/2-alpha)
            file.write(struct.pack('<i6di', 48, a_length, angle1, b_length, angle2, angle3, c_length, 48))
        length = struct.pack('<i', 4*len(positions))
        for i in range(3):
            file.write(length)
            data = array.array('f', (10*x[i] for x in positions))
            data.tofile(file)
            file.write(length)
        try:
            file.flush()
        except AttributeError:
            pass
