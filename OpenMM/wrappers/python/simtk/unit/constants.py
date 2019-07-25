#!/bin/env python
"""
Module simtk.unit.constants

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
from __future__ import print_function, division, absolute_import

__author__ = "Christopher M. Bruns"
__version__ = "0.5"

from .unit_definitions import *

#################
### CONSTANTS ###
#################

# codata 2006
AVOGADRO_CONSTANT_NA = 6.02214179e23 / mole
BOLTZMANN_CONSTANT_kB = 1.3806504e-23 * joule / kelvin
MOLAR_GAS_CONSTANT_R = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB

# From simtkcommon
SPEED_OF_LIGHT_C = 2.99792458e8 * meter / second
GRAVITATIONAL_CONSTANT_G = 6.6742e-11 * newton * meter**2 / kilogram**2

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
