"""
armberinpcrdfile.py: Used for loading AMBER inpcrd files.

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
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm.app.internal import amber_file_parser
from simtk.unit import Quantity, nanometers, picoseconds
try:
    import numpy
except:
    pass

class AmberInpcrdFile(object):
    """AmberInpcrdFile parses an AMBER inpcrd file and loads the data stored in it."""
    
    def __init__(self, file, loadVelocities=False, loadBoxVectors=False):
        """Load an inpcrd file.
        
        An inpcrd file contains atom positions and, optionally, velocities and periodic box dimensions.
        Unfortunately, it is sometimes impossible to determine from the file itself exactly what data
        it contains.  You therefore must specify in advance what data to load.  It is stored into this
        object's "positions", "velocities", and "boxVectors" fields.
        
        Parameters:
         - file (string) the name of the file to load
         - loadVelocities (boolean=False) whether to load velocities from the file
         - loadBoxVectors (boolean=False) whether to load the periodic box vectors
        """
        results = amber_file_parser.readAmberCoordinates(file, read_velocities=loadVelocities, read_box=loadBoxVectors)
        if loadVelocities:
            ## The atom positions read from the inpcrd file
            self.positions = results[0]
            if loadBoxVectors:
                ## The periodic box vectors read from the inpcrd file
                self.boxVectors = results[1]
                ## The atom velocities read from the inpcrd file
                self.velocities = results[2]
            else:
                self.velocities = results[1]
        elif loadBoxVectors:
            self.positions = results[0]
            self.boxVectors = results[1]
        else:
            self.positions = results
        self._numpyPositions = None
        if loadVelocities:
            self._numpyVelocities = None
        if loadBoxVectors:
            self._numpyBoxVectors = None

    def getPositions(self, asNumpy=False):
        """Get the atomic positions.
        
        Parameters:
         - asNumpy (boolean=False) if true, the values are returned as a numpy array instead of a list of Vec3s
         """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = Quantity(numpy.array(self.positions.value_in_unit(nanometers)), nanometers)
            return self._numpyPositions
        return self.positions
    
    def getVelocities(self, asNumpy=False):
        """Get the atomic velocities.
        
        Parameters:
         - asNumpy (boolean=False) if true, the vectors are returned as numpy arrays instead of Vec3s
         """
        if asNumpy:
            if self._numpyVelocities is None:
                self._numpyVelocities = Quantity(numpy.array(self.velocities.value_in_unit(nanometers/picoseconds)), nanometers/picoseconds)
            return self._numpyVelocities
        return self.velocities
    
    def getBoxVectors(self, asNumpy=False):
        """Get the periodic box vectors.
        
        Parameters:
         - asNumpy (boolean=False) if true, the values are returned as a numpy array instead of a list of Vec3s
         """
        if asNumpy:
            if self._numpyBoxVectors is None:
                self._numpyBoxVectors = []
                self._numpyBoxVectors.append(Quantity(numpy.array(self.boxVectors[0].value_in_unit(nanometers)), nanometers))
                self._numpyBoxVectors.append(Quantity(numpy.array(self.boxVectors[1].value_in_unit(nanometers)), nanometers))
                self._numpyBoxVectors.append(Quantity(numpy.array(self.boxVectors[2].value_in_unit(nanometers)), nanometers))
            return self._numpyBoxVectors
        return self.boxVectors
        
