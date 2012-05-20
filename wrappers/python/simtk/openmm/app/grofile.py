"""
grofile.py: Used for loading GRO files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Lee-Ping Wang
Contributors:

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
__author__ = "Lee-Ping Wang"
__version__ = "1.0"

import os
import sys
#import xml.etree.ElementTree as etree
#from copy import copy
from simtk.openmm import Vec3
#from simtk.openmm.app.internal.pdbstructure import PdbStructure
#from simtk.openmm.app import Topology
from re import sub, match
from simtk.unit import nanometers, angstroms#, is_quantity
#import element as elem
try:
    import numpy
except:
    pass

def isint(word):
    """ONLY matches integers! If you have a decimal point? None shall pass!

    @param[in] word String (for instance, '123', '153.0', '2.', '-354')
    @return answer Boolean which specifies whether the string is an integer (only +/- sign followed by digits)
    
    """
    return match('^[-+]?[0-9]+$',word)

def isfloat(word):
    """Matches ANY number; it can be a decimal, scientific notation, what have you
    CAUTION - this will also match an integer.

    @param[in] word String (for instance, '123', '153.0', '2.', '-354')
    @return answer Boolean which specifies whether the string is any number
    
    """
    return match('^[-+]?[0-9]*\.?[0-9]*([eEdD][-+]?[0-9]+)?$',word)

def is_gro_coord(line):
    """ Determines whether a line contains GROMACS data or not

    @param[in] line The line to be tested
    
    """
    sline = line.split()
    if len(sline) == 6 or len(sline) == 9:
        return all([isint(sline[2]),isfloat(sline[3]),isfloat(sline[4]),isfloat(sline[5])])
    elif len(sline) == 5 or len(sline) == 8:
        return all([isint(line[15:20]),isfloat(sline[2]),isfloat(sline[3]),isfloat(sline[4])])
    else:
        return 0

def is_gro_box(line):
    """ Determines whether a line contains a GROMACS box vector or not

    @param[in] line The line to be tested
    
    """
    sline = line.split()
    if len(sline) == 9 and all([isfloat(i) for i in sline]):
        return 1
    elif len(sline) == 3 and all([isfloat(i) for i in sline]):
        return 1
    else:
        return 0

class GroFile(object):
    """GroFile parses a GROMACS (.gro) file and constructs a set of atom positions from it."""
    
    def __init__(self, file):
        """Load a GRO file.
        
        The atom positions can be retrieved by calling getPositions().
        
        Parameters:
         - file (string) the name of the file to load
        """

        xyzs     = []
        elem     = [] # The element, most useful for quantum chemistry calculations
        atomname = [] # The atom name, for instance 'HW1'
        comms    = []
        resid    = []
        resname  = []
        boxes    = []
        xyz      = []
        ln       = 0
        lna      = 0
        frame    = 0
        for line in open(file):
            sline = line.split()
            if ln == 0:
                comms.append(line.strip())
            elif ln == 1:
                na = int(line.strip())
            elif is_gro_coord(line):
                if frame == 0: # Create the list of residues, atom names etc. only if it's the first frame.
                    # Name of the residue, for instance '153SOL1 -> SOL1' ; strips leading numbers
                    thisresname = sub('^[0-9]*','',sline[0])
                    resname.append(thisresname)
                    resid.append(int(sline[0].replace(thisresname,'')))
                    atomname.append(sline[1])
                    thiselem = sline[1]
                    if len(thiselem) > 1:
                        thiselem = thiselem[0] + sub('[A-Z0-9]','',thiselem[1:])
                    elem.append(thiselem)
                pos = [float(i) for i in sline[-3:]]
                xyz.append(Vec3(pos[0], pos[1], pos[2]))
            elif is_gro_box(line) and ln == na + 2:
                boxes.append([float(i) for i in sline]*nanometers)
                xyzs.append(xyz*nanometers)
                xyz = []
                ln = -1
                frame += 1
            else:
                print "Har, grofile.py encountered a line it didn't expect (line %i)!" % lna
                print line
            ln += 1
            lna += 1
            
        self.positions = xyzs
        self.elem = elem
        self.atomname = atomname
        self.resid = resid
        self.resname = resname
        self.boxes = boxes*nanometers
        self._numpyPositions = None
        
    def getPositions(self, asNumpy=False, index=0):
        """Get the atomic positions.
        
        Parameters:
         - asNumpy (boolean=False) if true, the values are returned as a numpy array instead of a list of Vec3s
         """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = numpy.array(self.positions[index].value_in_unit(nanometers))*nanometers
            return self._numpyPositions
        return self.positions[index]
    
