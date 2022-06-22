#!/bin/env python
"""
Module openmm.unit.unit_definitions

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2020 Stanford University and the Authors.
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

from __future__ import division, print_function, absolute_import

__author__ = "Christopher M. Bruns"
__version__ = "0.6"

from .baseunit import BaseUnit
from .standard_dimensions import *
from .unit import Unit, ScaledUnit, UnitSystem, dimensionless
from .unit_operators import * ; # needed for manipulation of units
from .prefix import *
import math
import sys

# Physical constants in this file are CODATA 2018 values from https://pml.nist.gov/cuu/Constants

#####################
### DIMENSIONLESS ###
#####################

# dimensionless = Unit({}); # defined in unit.py

##############
### LENGTH ###
##############

meter_base_unit = BaseUnit(length_dimension, "meter", "m")
meters = meter = Unit({meter_base_unit: 1.0})
define_prefixed_units(meter_base_unit, module = sys.modules[__name__])

angstrom_base_unit = BaseUnit(length_dimension, "angstrom", "A")
angstrom_base_unit.define_conversion_factor_to(meter_base_unit, 1e-10)
angstroms = angstrom = Unit({angstrom_base_unit: 1.0})

planck_length_base_unit = BaseUnit(length_dimension, "Planck length", "l_P")
planck_length_base_unit.define_conversion_factor_to(meter_base_unit, 1.616255e-35)

inch_base_unit = BaseUnit(length_dimension, "inch", "in")
inch_base_unit.define_conversion_factor_to(centimeter_base_unit, 2.5400)
inch = inches = Unit({inch_base_unit: 1.0})

foot_base_unit = BaseUnit(length_dimension, "foot", "ft")
foot_base_unit.define_conversion_factor_to(inch_base_unit, 12.0)
foot = feet = Unit({foot_base_unit: 1.0})

yard_base_unit = BaseUnit(length_dimension, "yard", "yd")
yard_base_unit.define_conversion_factor_to(foot_base_unit, 3.0)
yard = yards = Unit({yard_base_unit: 1.0})

furlongs = furlong = yard.create_unit(scale=220.0, name="furlong", symbol="furlong")
miles = mile = furlong.create_unit(scale=8.0, name="mile", symbol="mi")

bohrs = bohr = angstrom.create_unit(scale=0.529177210903, name='bohr', symbol='r_0')

############
### MASS ###
############

gram_base_unit = BaseUnit(mass_dimension, "gram", "g")
grams = gram = Unit({gram_base_unit: 1.0})

define_prefixed_units(gram_base_unit, module = sys.modules[__name__])

planck_mass_base_unit = BaseUnit(mass_dimension, "Planck mass", "m_P")
planck_mass_base_unit.define_conversion_factor_to(kilogram_base_unit, 2.176434e-8)

# pound can be mass, force, or currency
pound_mass_base_unit = BaseUnit(mass_dimension, "pound", "lb")
pound_mass_base_unit.define_conversion_factor_to(kilogram_base_unit, 0.3732)
pound_mass = pounds_mass = Unit({pound_mass_base_unit: 1.0})

stone_base_unit = BaseUnit(mass_dimension, "stone", "stone")
stone_base_unit.define_conversion_factor_to(pound_mass_base_unit, 14.0)
stone = stones = Unit({stone_base_unit: 1.0})

############
### TIME ###
############

second_base_unit = BaseUnit(time_dimension, "second", "s")
seconds = second = Unit({second_base_unit: 1.0})

define_prefixed_units(second_base_unit, module = sys.modules[__name__])

planck_time_base_unit = BaseUnit(time_dimension, "Planck time", "t_P")
planck_time_base_unit.define_conversion_factor_to(second_base_unit, 5.391247e-44)

minutes = minute = second.create_unit(scale=60.0, name="minute", symbol="min")
hours = hour = minute.create_unit(scale=60.0, name="hour", symbol="hr")
days = day = hour.create_unit(scale=24.0, name="day", symbol="day")
weeks = week = day.create_unit(scale=7.0, name="week", symbol="week")
years = year = day.create_unit(scale=365.25, name="julian year", symbol="a")
centuries = centurys = century = year.create_unit(scale=100.0, name="century", symbol="century")
millenia = milleniums = millenium = century.create_unit(scale=10.0, name="millenium", symbol="ka")

fortnights = fortnight = day.create_unit(scale=14.0, name="fortnight", symbol="fortnight")

###################
### TEMPERATURE ###
###################

kelvin_base_unit = BaseUnit(temperature_dimension, "kelvin", "K")
kelvins = kelvin = Unit({kelvin_base_unit: 1.0})

planck_temperature_base_unit = BaseUnit(temperature_dimension, "Planck temperature", "T_p")
planck_temperature_base_unit.define_conversion_factor_to(kelvin_base_unit, 1.416784e32)

##############
### CHARGE ###
##############

elementary_charge_base_unit = BaseUnit(charge_dimension, "elementary charge", "e")
elementary_charges = elementary_charge = Unit({elementary_charge_base_unit: 1.0})

coulomb_base_unit = BaseUnit(charge_dimension, "coulomb", "C")
# Exact conversion factor
coulomb_base_unit.define_conversion_factor_to(elementary_charge_base_unit, 6.241509074460763e18)
coulombs = coulomb = Unit({coulomb_base_unit: 1.0})

planck_charge_base_unit = BaseUnit(charge_dimension, "Planck charge", "q_P")
planck_charge_base_unit.define_conversion_factor_to(elementary_charge_base_unit, 1/math.sqrt(7.2973525693e-3)) # Calculated from fine structure constant

##############
### AMOUNT ###
##############

mole_base_unit = BaseUnit(amount_dimension, "mole", "mol")
moles = mole = Unit({mole_base_unit: 1.0})

single_item_amount_base_unit = BaseUnit(amount_dimension, "item", "")
mole_base_unit.define_conversion_factor_to(single_item_amount_base_unit, 6.02214076e23)
items = item = Unit({single_item_amount_base_unit: 1.0})

##########################
### Luminous Intensity ###
##########################

candela_base_unit = BaseUnit(luminous_intensity_dimension, "candela", "cd")
candelas = candela = Unit({candela_base_unit: 1.0})

#############
### ANGLE ###
#############

radian_base_unit = BaseUnit(angle_dimension, "radian", "rad")
radians = radian = Unit({radian_base_unit: 1.0})

degree_base_unit = BaseUnit(angle_dimension, "degree", "deg")
degree_base_unit.define_conversion_factor_to(radian_base_unit, math.pi/180.0)
degrees = degree = Unit({degree_base_unit: 1.0})

arcminutes = arcminute = degree.create_unit(scale=1.0/60.0, name="arcminute", symbol="'")
arcseconds = arcsecond = arcminute.create_unit(scale=1.0/60.0, name="arcsecond", symbol='"')

###################
### INFORMATION ###
###################

bit_base_unit = BaseUnit(information_dimension, "bit", "bit")
bits = bit = Unit({bit_base_unit: 1.0})

byte_base_unit = BaseUnit(information_dimension, "byte", "B")
byte_base_unit.define_conversion_factor_to(bit_base_unit, 8.0)
bytes = byte = Unit({byte_base_unit: 1.0})

nat_base_unit = BaseUnit(information_dimension, "nat", "nat")
nat_base_unit.define_conversion_factor_to(bit_base_unit, 1.0/math.log(2.0))
nats = nat = nits = nit = nepits = nepit = Unit({nat_base_unit: 1.0})

ban_base_unit = BaseUnit(information_dimension, "ban", "ban")
ban_base_unit.define_conversion_factor_to(bit_base_unit, math.log(10.0, 2.0))
bans = ban = hartleys = hartley = dits = dit = Unit({ban_base_unit: 1.0})


###############
### DERIVED ###
###############

# Molar mass
# daltons = dalton = grams / mole
daltons = dalton = Unit({ScaledUnit(1.0, gram/mole, "dalton", "Da"): 1.0})
amus = amu = dalton
atom_mass_units = atomic_mass_unit = dalton

# Volume
liter_base_unit = ScaledUnit(1.0, decimeter**3, "liter", "L")
liter = liters = litre = litres = Unit({liter_base_unit: 1.0})
define_prefixed_units(liter_base_unit, module = sys.modules[__name__])

# Concentration
molar_base_unit = ScaledUnit(1.0, mole/liter, "molar", "M")
molar = molal = Unit({molar_base_unit: 1.0})
define_prefixed_units(molar_base_unit, module = sys.modules[__name__])

# Force
newton_base_unit = ScaledUnit(1.0, kilogram * meter / second / second, "newton", "N")
newtons = newton = Unit({newton_base_unit: 1.0})
define_prefixed_units(newton_base_unit, module = sys.modules[__name__])
# pound can be mass, force, or currency
pound_force_base_unit = ScaledUnit(4.448, newton, "pound", "lb")
pound_force = pounds_force = Unit({pound_force_base_unit: 1.0})
dyne_base_unit = ScaledUnit(1.0, gram * centimeter / second**2, "dyne", "dyn")
dyne = dynes = Unit({dyne_base_unit: 1.0})

# Energy
joule_base_unit = ScaledUnit(1.0, newton * meter, "joule", "J")
joules = joule = Unit({joule_base_unit: 1.0})
define_prefixed_units(joule_base_unit, module = sys.modules[__name__])
erg_base_unit = ScaledUnit(1.0, dyne * centimeter, "erg", "erg")
erg = ergs = Unit({erg_base_unit: 1.0})
hartree_base_unit = ScaledUnit(4.3597447222071e-18, joule, "hartree", "Ha")
hartree = hartrees = Unit({hartree_base_unit: 1.0})

# In molecular simulations, "kilojoules" are in microscopic units
# And you really only want to use kilojoules/mole.
md_kilojoule_raw = gram * nanometer**2 / picosecond**2
md_kilojoules = md_kilojoule = Unit({ScaledUnit(1.0, md_kilojoule_raw, "kilojoule", "kJ"): 1.0})
kilojoules_per_mole = kilojoule_per_mole = md_kilojoule / mole

calorie_base_unit = ScaledUnit(4.184, joule, "calorie", "cal")
calories = calorie = Unit({calorie_base_unit: 1.0})
define_prefixed_units(calorie_base_unit, module = sys.modules[__name__])
md_kilocalories = md_kilocalorie = Unit({ScaledUnit(4.184, md_kilojoule, "kilocalorie", "kcal"): 1.0})
kilocalories_per_mole = kilocalorie_per_mole = md_kilocalorie / mole

# Power
watt_base_unit = ScaledUnit(1.0, joule / second, "watt", "W")
watt = watts = Unit({watt_base_unit: 1.0})

# Current
ampere_base_unit = ScaledUnit(1.0, coulomb / second, "ampere", "A")
ampere = amperes = amp = amps = Unit({ampere_base_unit: 1.0})

# Electrical potential
volt_base_unit = ScaledUnit(1.0, watt / ampere, "volt", "V")
volt = volts = Unit({volt_base_unit: 1.0})

# Magnetic field
tesla_base_unit = ScaledUnit(1.0, newton / (ampere * meter), "tesla", "T")
tesla = teslas = Unit({tesla_base_unit: 1.0})
gauss_base_unit = ScaledUnit(10.0**-4, tesla, "gauss", "G")
gauss = Unit({gauss_base_unit: 1.0})

# Electrical resistance
ohm_base_unit = ScaledUnit(1.0, volt / ampere, "ohm", "O")
ohm = ohms = Unit({ohm_base_unit: 1.0})

# Capacitance
farad_base_unit = ScaledUnit(1.0, coulomb / volt, "farad", "F")
farad = farads = Unit({farad_base_unit: 1.0})

# Inductance
henry_base_unit = ScaledUnit(1.0, volt * second / ampere, "henry", "H")
henry = henries = henrys = Unit({henry_base_unit: 1.0})

# Polarizability
debye_base_unit = ScaledUnit(0.20822678, elementary_charge * angstrom, "debye", "D")
debyes = debye = Unit({debye_base_unit: 1.0})

# Pressure
pascal_base_unit = ScaledUnit(1.0, newton / (meter**2), "pascal", "Pa")
pascals = pascal = Unit({pascal_base_unit: 1.0})
define_prefixed_units(pascal_base_unit, module = sys.modules[__name__])
psi_base_unit = ScaledUnit(1.0, pound_force / (inch**2), "psi", "psi")
psi = Unit({psi_base_unit: 1.0})
bar_base_unit = ScaledUnit(10.0**5, pascal, "bar", "bar")
bar = bars = Unit({bar_base_unit: 1.0})
atmosphere_base_unit = ScaledUnit(101325.0, pascal, "atmosphere", "atm")
atmosphere = atmospheres = Unit({atmosphere_base_unit: 1.0})
torr_base_unit = ScaledUnit(1.0/760.0, atmosphere, "torr", "Torr")
torr = Unit({torr_base_unit: 1.0})
mmHg_base_unit = ScaledUnit(133.322, pascal, "mmHg", "mmHg")
mmHg = Unit({mmHg_base_unit: 1.0})

####################
### Unit Systems ###
####################

ampere_base_unit = ScaledUnit(1.0, coulomb/second, "ampere", "A")

si_unit_system = UnitSystem([
        meter_base_unit,
        kilogram_base_unit,
        second_base_unit,
        ampere_base_unit,
        kelvin_base_unit,
        mole_base_unit,
        candela_base_unit,
        radian_base_unit])

cgs_unit_system = UnitSystem([
        centimeter_base_unit,
        gram_base_unit,
        second_base_unit,
        ampere_base_unit,
        kelvin_base_unit,
        mole_base_unit,
        radian_base_unit])

dalton_base_unit = ScaledUnit(1.0, gram/mole, "dalton", "Da")

md_unit_system = UnitSystem([
        nanometer_base_unit,
        dalton_base_unit,
        picosecond_base_unit,
        elementary_charge_base_unit,
        kelvin_base_unit,
        mole_base_unit,
        radian_base_unit])

planck_unit_system = UnitSystem([\
        planck_length_base_unit,
        planck_mass_base_unit,
        planck_time_base_unit,
        planck_charge_base_unit,
        planck_temperature_base_unit,
        single_item_amount_base_unit,
        radian_base_unit])

########################
### TESTING/EXAMPLES ###
########################

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
