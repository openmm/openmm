"""
xmlfile.py: Used for loading XML files.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2015 Stanford University and the Authors.
Authors: Carlos Xavier Hernandez
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
__author__ = "Carlos Xavier Hernandez"
__version__ = "1.0"
from simtk.openmm import XmlSerializer
from simtk.openmm.openmm import State, System, Integrator


class XMLFile(object):
    """XMLFile deserializes an OpenMM XML file and returns its contents.

    This base class also provides a method for subclasses to test whether
    the XML returns an expected OpenMM class."""

    def __new__(cls, filename):
        obj = cls._load_xml_(filename)
        assert isinstance(obj, cls.xmlcls), "File does not contain %s." % cls.xmlcls
        return obj

    @staticmethod
    def _load_xml_(filename):
        with open(filename, 'rb') as f:
            xmlfile = XmlSerializer.deserialize(f.read())
        return xmlfile


class XMLStateFile(XMLFile):
    """XMLStateFile deserializes an XML file containing an OpenMM State object."""
    xmlcls = State


class XMLSystemFile(XMLFile):
    """XMLSystemFile deserializes an XML file containing an OpenMM System object."""
    xmlcls = System


class XMLIntegratorFile(XMLFile):
    """XMLIntegratorFile deserializes an XML file containing an OpenMM Integrator object."""
    xmlcls = Integrator
