#!/bin/env python
"""
Module simtk.unit.prefix

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
__version__ = "0.6"

from .baseunit import BaseUnit
from .unit import Unit, ScaledUnit
import sys

###################
### SI PREFIXES ###
###################

class SiPrefix(object):
    """
    Unit prefix that can be multiplied by a unit to yield a new unit.

    e.g. millimeter = milli*meter
    """
    def __init__(self, prefix, factor, symbol):
        self.prefix = prefix
        self.factor = factor
        self.symbol = symbol

    def __mul__(self, unit):
        """
        SiPrefix * BaseUnit yields new BaseUnit
        SiPrefix * ScaledUnit yields new ScaledUnit
        SiPrefix * Unit with exactly one BaseUnit or ScaledUnit yields new Unit
        """
        if isinstance(unit, BaseUnit):
            # BaseUnit version
            symbol = self.symbol + unit.symbol
            name = self.prefix + unit.name
            factor = self.factor
            # TODO - check for existing BaseUnit with same name, symbol, and factor
            new_base_unit = BaseUnit(unit.dimension, name, symbol)
            new_base_unit.define_conversion_factor_to(unit, factor)
            return new_base_unit
        elif isinstance(unit, ScaledUnit):
            # ScaledUnit version
            symbol = self.symbol + unit.symbol
            name = self.prefix + unit.name
            factor = self.factor * unit.factor
            # TODO - check for existing BaseUnit with same name, symbol, and factor
            return ScaledUnit(factor, unit.master, name, symbol)
        elif isinstance(unit, Unit):
            base_units = list(unit.iter_base_or_scaled_units())
            if 1 != len(base_units):
                raise TypeError('Unit prefix "%s" can only be with simple Units containing one component.' % self.prefix)
            if 1 != base_units[0][1]:
                raise TypeError('Unit prefix "%s" can only be with simple Units with an exponent of 1.' % self.prefix)
            base_unit = base_units[0][0]
            # Delegate to Base or Scaled Unit multiply
            new_base_unit = self * base_unit
            new_unit = Unit({new_base_unit: 1.0})
            return new_unit
        else:
            raise TypeError('Unit prefix "%s" can only be applied to a Unit, BaseUnit, or ScaledUnit.' % self.prefix)

yotto = SiPrefix("yotto", 1e-24, 'y')
zepto = SiPrefix("zepto", 1e-21, 'z')
atto  = SiPrefix("atto",  1e-18, 'a')
femto = SiPrefix("femto", 1e-15, 'f')
pico  = SiPrefix("pico",  1e-12, 'p')
nano  = SiPrefix("nano",  1e-9,  'n')
micro = SiPrefix("micro", 1e-6,  'u')
milli = SiPrefix("milli", 1e-3,  'm')
centi = SiPrefix("centi", 1e-2,  'c')
deci  = SiPrefix("deci",  1e-1,  'd')
# two spellings for deka prefix
deka  = SiPrefix("deka",  1e1,   'da')
deca  = SiPrefix("deca",  1e1,   'da')
hecto = SiPrefix("hecto", 1e2,   'h')
kilo  = SiPrefix("kilo",  1e3,   'k')
mega  = SiPrefix("mega",  1e6,   'M')
giga  = SiPrefix("giga",  1e9,   'G')
tera  = SiPrefix("tera",  1e12,  'T')
peta  = SiPrefix("peta",  1e15,  'P')
exa   = SiPrefix("exa",   1e18,  'E')
zetta = SiPrefix("zetta", 1e21,  'Z')
yotta = SiPrefix("yotta", 1e24,  'Y')

si_prefixes = (   yotto
                , zepto
                , atto
                , femto
                , pico
                , nano
                , micro
                , milli
                , centi
                , deci
                , deka
                , deca
                , hecto
                , kilo
                , mega
                , giga
                , tera
                , peta
                , exa
                , zetta
                , yotta)

def define_prefixed_units(base_unit, module = sys.modules[__name__]):
    """
    Create attributes for prefixed units derived from a particular BaseUnit, e.g. "kilometer" from "meter_base_unit"

    Parameters
     - base_unit (BaseUnit) existing base unit to use as a basis for prefixed units
     - module (Module) module which will contain the new attributes.  Defaults to simtk.unit module.
    """
    for prefix in si_prefixes:
        new_base_unit = prefix * base_unit
        name = new_base_unit.name
        new_unit = Unit({new_base_unit: 1.0})
        # Create base_unit attribute, needed for creating UnitSystems
        module.__dict__[name + '_base_unit'] = new_base_unit # e.g. "kilometer_base_unit"
        # Create attribue in this module
        module.__dict__[name] = new_unit # e.g. "kilometer"
        # And plural version
        module.__dict__[name + 's'] = new_unit # e.g. "kilometers"


# Binary prefixes

kibi = SiPrefix("kibi", 2.0**10,  'Ki')
mebi = SiPrefix("mebi", 2.0**20,  'Mi')
gibi = SiPrefix("gibi", 2.0**30,  'Gi')
tebi = SiPrefix("tebi", 2.0**40,  'Ti')
pebi = SiPrefix("pebi", 2.0**50,  'Pi')
exbi = SiPrefix("exbi", 2.0**60,  'Ei')
zebi = SiPrefix("zebi", 2.0**70,  'Zi')
yobi = SiPrefix("yobi", 2.0**80,  'Yi')

binary_prefixes = ( kibi
                  , mebi
                  , gibi
                  , tebi
                  , pebi
                  , exbi
                  , zebi
                  , yobi)


# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
