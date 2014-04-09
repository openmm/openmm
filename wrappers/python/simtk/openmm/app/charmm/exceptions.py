"""
This package contains exceptions that may be raised by the CHARMM components of
the OpenMM Application layer

Author: Jason M. Swails
Contributors:
Date: April 9, 2014
"""

class CharmmError(Exception):
    """ Base class for all exceptions raised in this package """

class CharmmWarning(Warning):
    """ Base class for all warnings emitted in this package """

class CharmmPSFError(CharmmError):
    """ If there is a problem parsing CHARMM PSF files """

class SplitResidueWarning(CharmmWarning):
    """ For if a residue with the same number but different names is split """

class ResidueError(CharmmError):
    """ For when there are problems defining a residue """

class CharmmPSFWarning(CharmmWarning):
    """ For non-fatal PSF parsing issues """

class CharmmFileError(CharmmError):
    """ If there is a problem parsing CHARMM files """

class MissingParameter(CharmmError):
    """ If a parameter is missing from a database """

class CmapError(CharmmError):
    """ For an error arising from CMAP grid incompatibilities """

class BondError(CharmmError):
    """ Prevent an atom from bonding to itself """
