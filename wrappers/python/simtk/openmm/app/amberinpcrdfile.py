"""
armberinpcrdfile.py: Used for loading AMBER inpcrd files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2014 Stanford University and the Authors.
Authors: Peter Eastman
Contributors: Jason Swails

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

from functools import wraps
from simtk.openmm.app.internal import amber_file_parser
from simtk.unit import Quantity, nanometers, picoseconds
import warnings
try:
    import numpy as np
except:
    np = None

def numpy_protector(func):
    """
    Decorator to emit useful error messages if users try to request numpy
    processing if numpy is not available. Raises ImportError if numpy could not
    be found
    """
    @wraps(func)
    def wrapper(self, asNumpy=False):
        if asNumpy and np is None:
            raise ImportError('Could not import numpy. Cannot set asNumpy=True')
        return func(self, asNumpy=asNumpy)
    return wrapper

class AmberInpcrdFile(object):
    """AmberInpcrdFile parses an AMBER inpcrd file and loads the data stored in it."""

    def __init__(self, file, loadVelocities=None, loadBoxVectors=None):
        """Load an inpcrd file.

        An inpcrd file contains atom positions and, optionally, velocities and
        periodic box dimensions.

        Parameters
        ----------
        file : str
             The name of the file to load
        loadVelocities : bool
             Deprecated. Velocities are loaded automatically if present
        loadBoxVectors : bool
            Deprecated. Box vectors are loaded automatically if present
        """
        self.file = file
        if loadVelocities is not None or loadBoxVectors is not None:
            warnings.warn('loadVelocities and loadBoxVectors have been '
                          'deprecated. velocities and box information '
                          'is loaded automatically if the inpcrd file contains '
                          'them.', DeprecationWarning)
        results = amber_file_parser.readAmberCoordinates(file)
        self.positions, self.velocities, self.boxVectors = results
        # Cached numpy arrays
        self._numpyPositions = None
        self._numpyVelocities = None
        self._numpyBoxVectors = None

    @numpy_protector
    def getPositions(self, asNumpy=False):
        """Get the atomic positions.

        Parameters
        ----------
        asNumpy : bool=False
            if true, the values are returned as a numpy array instead of a list
            of Vec3s
        """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = Quantity(np.array(self.positions.value_in_unit(nanometers)), nanometers)
            return self._numpyPositions
        return self.positions

    @numpy_protector
    def getVelocities(self, asNumpy=False):
        """Get the atomic velocities.

        Parameters
        ----------
        asNumpy : bool=False
            if true, the vectors are returned as numpy arrays instead of Vec3s
        """
        if self.velocities is None:
            raise AttributeError('velocities not found in %s' % self.file)
        if asNumpy:
            if self._numpyVelocities is None:
                self._numpyVelocities = Quantity(np.array(self.velocities.value_in_unit(nanometers/picoseconds)), nanometers/picoseconds)
            return self._numpyVelocities
        return self.velocities

    @numpy_protector
    def getBoxVectors(self, asNumpy=False):
        """Get the periodic box vectors.

        Parameters
        ----------
        asNumpy : bool=False
            if true, the values are returned as a numpy array instead of a list
            of Vec3s
        """
        if self.boxVectors is None:
            raise AttributeError('Box information not found in %s' % self.file)
        if asNumpy:
            if self._numpyBoxVectors is None:
                self._numpyBoxVectors = []
                self._numpyBoxVectors.append(Quantity(np.array(self.boxVectors[0].value_in_unit(nanometers)), nanometers))
                self._numpyBoxVectors.append(Quantity(np.array(self.boxVectors[1].value_in_unit(nanometers)), nanometers))
                self._numpyBoxVectors.append(Quantity(np.array(self.boxVectors[2].value_in_unit(nanometers)), nanometers))
            return self._numpyBoxVectors
        return self.boxVectors

