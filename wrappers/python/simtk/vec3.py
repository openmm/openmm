#!/bin/env python


"""
vec3.py: Used for managing vectors of lenght three.
"""
__author__ = "Christopher M. Bruns"
__version__ = "1.0"


import sys
import simtk.unit

class Vec3(object):
    """
    Vec3 represents a three-dimensional point or displacement in space.
    
    Examples
    
    >>> v = Vec3((1, 2, 3))
    >>> print v.dot(v)
    14
    >>> print abs(v)
    3.74165738677
    >>> v += Vec3((3,2,1))
    >>> print v
    [4, 4, 4]
    >>> print v**2
    48.0
    >>> print v**1
    6.92820323028
    >>> print v
    [4, 4, 4]
    """
    
    def __init__(self, xyz):
        """
        Create a new Vec3 object, specifying its three elements.
        
        
        Examples:
        
        >>> v = Vec3((1,2,3))
        >>> print v
        [1, 2, 3]
        
        """
        x = xyz[0]
        y = xyz[1]
        z = xyz[2]
        self.data = [x, y, z]

    def __repr__(self):
        """
        Create a string containing a python code representation of this Vec3.
        
        Examples
        
        >>> v = Vec3([1,2,3])
        >>> v
        Vec3([1, 2, 3])
        """
        return self.__class__.__name__ + \
                '([' + repr(self[0]) + ", " + repr(self[1]) + ", " + repr(self[2]) + '])'

    def __len__(self):
        """
        Returns the number of elements in this Vec3.  Always 3.
        
        Examples
        
        >>> v = Vec3((1, 2, 3))
        >>> len(v)
        3
        """
        return 3
        
    def __iter__(self):
        """
        Examples
        
        >>> v = Vec3((1,2,3))
        >>> for c in v:
        ...   print c
        ... 
        1
        2
        3
        """
        for coord in self.data:
            yield coord
        
    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __abs__(self):
        return pow(self.dot(self), 0.5)
        
    def __neg__(self):
        return self.__class__([-self[0], -self[1], -self[2]])
        
    def __pos__(self):
        return self
        
    def __add__(self, other):
        return self.__class__([self[0] + other[0], self[1] + other[1], self[2] + other[2]])

    def __radd__(self, other):
        return self.__class__([other[0] + self[0], other[1] + self[1], other[2] + self[2]])

    def __sub__(self, other):
        return self.__class__([self[0] - other[0], self[1] - other[1], self[2] - other[2]])

    def __rsub__(self, other):
        return self.__class__([other[0] - self[0], other[1] - self[1], other[2] - self[2]])

    def __mul__(self, other):
        """
        Returns NotImplemented if the right hand side is a simtk.unit.Unit,
        so the Unit class can return a proper Quantity unit.
        
        Example
        
        >>> import simtk.unit
        >>> Vec3([1,2,3]) * simtk.unit.meter
        Quantity(Vec3([1, 2, 3]), meter)
        """
        # If other is a Unit, delegate to Unit.__rmul__ multiplication
        if simtk.unit.is_unit(other):
            return NotImplemented
        return self.__class__([self[0] * other, self[1] * other, self[2] * other])

    def __rmul__(self, other):
        return self.__class__([other * self[0], other * self[1], other * self[2]])

    def __div__(self, other):
        # If other is a Unit, delegate to Unit.__rmul__ multiplication
        if simtk.unit.is_unit(other):
            return NotImplemented
        return self.__class__([self[0] / other, self[1] / other, self[2] / other])
            
    def __rdiv__(self, other):
        return self.__class__([other / self[0], other / self[1], other / self[2]])
            
    def __str__(self):
        return "["+str(self[0])+", "+str(self[1])+", "+str(self[2])+"]"

    def __pow__(self, exponent):
        return pow(self.dot(self), exponent * 0.5)

    def dot(self, other):
        assert len(other) == 3
        return self[0] * other[0] + self[1] * other[1] + self[2] * other[2]


# run module directly for testing
if __name__=='__main__':
    
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
