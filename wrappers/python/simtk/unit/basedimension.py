#!/bin/env python


"""
Module simtk.unit.basedimension

BaseDimension class for use by units and quantities.
BaseDimensions are things like "length" and "mass".
"""

__author__ = "Christopher M. Bruns"
__version__ = "0.6"


class BaseDimension(object):
    '''
    A physical dimension such as length, mass, or temperature.
    
    It is unlikely the user will need to create new ones.
    '''
    # Keep deterministic order of dimensions
    _index_by_name = {
        'mass': 1,
        'length': 2,
        'time': 3,
        'temperature': 4,
        'amount': 5,
        'charge': 6,
        'luminous intensity': 7,
        'angle': 8,
    }
    _next_unused_index = 9
    
    def __init__(self, name):
        """Create a new BaseDimension.

        Each new BaseDimension is assumed to be independent of all other BaseDimensions.
        Use the existing BaseDimensions in simtk.dimension instead of creating
        new ones.
        """
        self.name = name
        if not self.name in BaseDimension._index_by_name.keys():
            BaseDimension._index_by_name[name] = BaseDimension._next_unused_index
            BaseDimension._next_unused_index += 1
        self._index = BaseDimension._index_by_name[name]
        
    def __lt__(self, other):
        """
        The implicit order of BaseDimensions is the order in which they were created.
        This method is used for using BaseDimensions as hash keys, and also affects
        the order in which units appear in multi-dimensional Quantities.
        
        Returns True if self < other, False otherwise.
        """
        return self._index < other._index
        
    def __hash__(self):
        """
        Needed for using BaseDimensions as hash keys.
        """
        return self._index

    def __repr__(self):
        return 'BaseDimension("%s")' % self.name


# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])

