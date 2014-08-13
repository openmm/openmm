"""
This package contains exceptions that may be raised by the CHARMM components of
the OpenMM Application layer

This file is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of Biological
Structures at Stanford, funded under the NIH Roadmap for Medical Research,
grant U54 GM072970. See https://simtk.org.  This code was originally part of
the ParmEd program and was ported for use with OpenMM.

Copyright (c) 2014 the Authors

Author: Jason M. Swails
Contributors:
Date: April 18, 2014

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

class CharmmError(Exception):
    """ Base class for all exceptions raised in this package """

class CharmmWarning(Warning):
    """ Base class for all warnings emitted in this package """

class CharmmPSFError(CharmmError):
    """ If there is a problem parsing CHARMM PSF files """

class CharmmPsfEOF(CharmmError):
    """ Raised when the end-of-file is reached on a CHARMM PSF file """

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

class MoleculeError(CharmmError):
    """ For (impossibly) messed up connectivity """
