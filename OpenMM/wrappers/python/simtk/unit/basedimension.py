#!/bin/env python
"""
Module simtk.unit.basedimension

BaseDimension class for use by units and quantities.
BaseDimensions are things like "length" and "mass".

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Christopher M. Bruns
Contributors: Peter Eastman

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
from __future__ import print_function, division

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
    
    def __eq__(self, other):
        if isinstance(other, BaseDimension):
            return self._index == other._index
        return False
    
    def __ne__(self, other):
        if isinstance(other, BaseDimension):
            return self._index != other._index
        return False


# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])

