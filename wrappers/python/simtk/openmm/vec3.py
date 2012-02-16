"""
vec3.py: Defines the Vec3 class used by OpenMM
"""
__author__ = "Peter Eastman"
__version__ = "1.0"

import simtk.unit as unit

class Vec3(tuple):
    """Vec3 is a 3-element tuple that supports many math operations."""
    
    def __new__(cls, x, y, z):
        """Create a new Vec3."""
        return tuple.__new__(cls, (x, y, z))
    
    def __add__(self, other):
        """Add two Vec3s."""
        return Vec3(self[0]+other[0], self[1]+other[1], self[2]+other[2])
    
    def __radd__(self, other):
        """Add two Vec3s."""
        return Vec3(self[0]+other[0], self[1]+other[1], self[2]+other[2])
    
    def __sub__(self, other):
        """Add two Vec3s."""
        return Vec3(self[0]-other[0], self[1]-other[1], self[2]-other[2])
    
    def __rsub__(self, other):
        """Add two Vec3s."""
        return Vec3(other[0]-self[0], other[1]-self[1], other[2]-self[2])
    
    def __mul__(self, other):
        """Multiply a Vec3 by a constant."""
        if unit.is_unit(other):
            return unit.Quantity(self, other)
        return Vec3(other*self[0], other*self[1], other*self[2])
    
    def __rmul__(self, other):
        """Multiply a Vec3 by a constant."""
        if unit.is_unit(other):
            return unit.Quantity(self, other)
        return Vec3(other*self[0], other*self[1], other*self[2])
    
    def __div__(self, other):
        """Divide a Vec3 by a constant."""
        return Vec3(self[0]/other, self[1]/other, self[2]/other)
    __truediv__ = __div__
    
    def __deepcopy__(self, memo):
        return Vec3(self[0], self[1], self[2])
