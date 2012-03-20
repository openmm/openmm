"""
dcdfile.py: Used for writing DCD files.
"""
__author__ = "Peter Eastman"
__version__ = "1.0"

import array
import os
import time
import struct
from simtk.unit import picoseconds, nanometers, angstroms, is_quantity
    
class DCDFile(object):
    """DCDFile provides methods for creating DCD files.
    
    DCD is a file format for storing simulation trajectories.  It is supported by many programs, such
    as CHARMM, NAMD, and X-PLOR.  Note, however, that different programs produce subtly different
    versions of the format.  This class generates the CHARMM version.  Also note that there is no
    standard byte ordering (big-endian or little-endian) for this format.  This class always generates
    files with little-endian ordering.
    
    To use this class, create a DCDFile object, then call writeModel() once for each model in the file."""
    
    def __init__(self, file, topology, dt, firstStep=0, interval=1):
        """Create a DCD file and write out the header.
        
        Parameters:
         - file (file) A file to write to
         - topology (Topology) The Topology defining the molecular system being written
         - dt (time) The time step used in the trajectory
         - firstStep (int=0) The index of the first step in the trajectory
         - interval (int=1) The frequency (measured in time steps) at which states are written to the trajectory
        """
        self._file = file
        self._topology = topology
        self._firstStep = firstStep
        self._interval = interval
        self._modelCount = 0
        if is_quantity(dt):
            dt = dt.value_in_unit(picoseconds)
        dt /= 0.04888821
        boxFlag = 0
        if topology.getUnitCellDimensions() is not None:
            boxFlag = 1
        header = struct.pack('<i4c9if', 84, 'C', 'O', 'R', 'D', 0, firstStep, interval, 0, 0, 0, 0, 0, 0, dt)
        header += struct.pack('<13i', boxFlag, 0, 0, 0, 0, 0, 0, 0, 0, 24, 84, 164, 2)
        header += struct.pack('<80s', 'Created by OpenMM')
        header += struct.pack('<80s', 'Created '+time.asctime(time.localtime(time.time())))
        header += struct.pack('<4i', 164, 4, len(list(topology.atoms())), 4)
        file.write(header)
    
    def writeModel(self, positions, unitCellDimensions=None):
        """Write out a model to the DCD file.
                
        Parameters:
         - positions (list) The list of atomic positions to write
         - unitCellDimensions (Vec3=None) The dimensions of the crystallographic unit cell.  If None, the dimensions specified in
           the Topology will be used.  Regardless of the value specified, no dimensions will be written if the Topology does not
           represent a periodic system.
        """
        if len(list(self._topology.atoms())) != len(positions):
            raise ValueError('The number of positions must match the number of atoms') 
        if is_quantity(positions):
            positions = positions.value_in_unit(nanometers)
        file = self._file
        
        # Update the header.
        
        self._modelCount += 1
        file.seek(8, os.SEEK_SET)
        file.write(struct.pack('<i', self._modelCount))
        file.seek(20, os.SEEK_SET)
        file.write(struct.pack('<i', self._firstStep+self._modelCount*self._interval))
        
        # Write the data.
        
        file.seek(0, os.SEEK_END)
        boxSize = self._topology.getUnitCellDimensions()
        if boxSize is not None:
            if unitCellDimensions is not None:
                boxSize = unitCellDimensions
            size = boxSize.value_in_unit(angstroms)
            file.write(struct.pack('<i6di', 48, size[0], 0, size[1], 0, 0, size[2], 48))
        length = struct.pack('<i', 4*len(positions))
        for i in range(3):
            file.write(length)
            data = array.array('f', (10*x[i] for x in positions))
            data.tofile(file)
            file.write(length)
