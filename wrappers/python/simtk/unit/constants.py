#!/bin/env python
"""
Module simtk.unit.constants
"""

from __future__ import division

__author__ = "Christopher M. Bruns"
__version__ = "0.5"

from unit_definitions import *

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
