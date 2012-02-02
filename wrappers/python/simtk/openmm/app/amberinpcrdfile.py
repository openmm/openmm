"""
armberinpcrdfile.py: Used for loading AMBER inpcrd files.
"""
__author__ = "Peter Eastman"
__version__ = "1.0"

from simtk.openmm.app.internal import amber_file_parser
from simtk.unit import nanometers, picoseconds
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
            self.positions = results[0]
            if loadBoxVectors:
                self.boxVectors = results[1]
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
                self._numpyPositions = numpy.array(self.positions.value_in_unit(nanometers))*nanometers
            return self._numpyPositions
        return self.positions
    
    def getVelocities(self, asNumpy=False):
        """Get the atomic velocities.
        
        Parameters:
         - asNumpy (boolean=False) if true, the vectors are returned as numpy arrays instead of Vec3s
         """
        if asNumpy:
            if self._numpyVelocities is None:
                self._numpyVelocities = numpy.array(self.velocities.value_in_unit(nanometers/picoseconds))*nanometers/picoseconds
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
                self._numpyBoxVectors.append(numpy.array(self.boxVectors[0].value_in_unit(nanometers))*nanometers)
                self._numpyBoxVectors.append(numpy.array(self.boxVectors[1].value_in_unit(nanometers))*nanometers)
                self._numpyBoxVectors.append(numpy.array(self.boxVectors[2].value_in_unit(nanometers))*nanometers)
            return self._numpyBoxVectors
        return self.boxVectors
        