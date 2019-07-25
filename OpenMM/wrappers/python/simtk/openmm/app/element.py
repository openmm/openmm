"""
element.py: Used for managing elements.

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Christopher M. Bruns
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
from __future__ import absolute_import
__author__ = "Christopher M. Bruns"
__version__ = "1.0"

import sys
from collections import OrderedDict
from simtk.unit import daltons, is_quantity
if sys.version_info >= (3, 0):
    import copyreg
else:
    import copy_reg as copyreg


class Element(object):
    """An Element represents a chemical element.

    The simtk.openmm.app.element module contains objects for all the standard chemical elements,
    such as element.hydrogen or element.carbon.  You can also call the static method Element.getBySymbol() to
    look up the Element with a particular chemical symbol.

    Element objects should be considered immutable
    """

    _elements_by_symbol = {}
    _elements_by_atomic_number = {}
    _elements_by_mass = None

    def __init__(self, number, name, symbol, mass):
        """Create a new element

        Parameters
        ----------
        number : int
            The atomic number of the element
        name : string
            The name of the element
        symbol : string
            The chemical symbol of the element
        mass : float
            The atomic mass of the element
        """
        ## The atomic number of the element
        self._atomic_number = number
        ## The name of the element
        self._name = name
        ## The chemical symbol of the element
        self._symbol = symbol
        ## The atomic mass of the element
        self._mass = mass
        # Index this element in a global table
        s = symbol.strip().upper()
        ## If we add a new element, we need to re-hash elements by mass
        Element._elements_by_mass = None

        if s in Element._elements_by_symbol:
            raise ValueError('Duplicate element symbol %s' % s)
        Element._elements_by_symbol[s] = self
        if number in Element._elements_by_atomic_number:
            other_element = Element._elements_by_atomic_number[number]
            if mass < other_element.mass:
                # If two "elements" share the same atomic number, they're
                # probably hydrogen and deuterium, and we want to choose
                # the lighter one to put in the table by atomic_number,
                # since it's the "canonical" element.
                Element._elements_by_atomic_number[number] = self
        else:
            Element._elements_by_atomic_number[number] = self

    @staticmethod
    def getBySymbol(symbol):
        """Get the Element with a particular chemical symbol."""
        s = symbol.strip().upper()
        return Element._elements_by_symbol[s]

    @staticmethod
    def getByAtomicNumber(atomic_number):
        return Element._elements_by_atomic_number[atomic_number]

    @staticmethod
    def getByMass(mass):
        """
        Get the element whose mass is CLOSEST to the requested mass. This method
        should not be used for repartitioned masses

        Parameters
        ----------
        mass : float or Quantity
            Mass of the atom to find the element for. Units assumed to be
            daltons if not specified

        Returns
        -------
        Element
            The element whose atomic mass is closest to the input mass
        """
        # Assume masses are in daltons if they are not units
        if is_quantity(mass):
            mass = mass.value_in_unit(daltons)
        if mass < 0:
            raise ValueError('Invalid Higgs field')
        # If this is our first time calling getByMass (or we added an element
        # since the last call), re-generate the ordered by-mass dict cache
        if Element._elements_by_mass is None:
            Element._elements_by_mass = OrderedDict()
            for elem in sorted(Element._elements_by_symbol.values(),
                               key=lambda x: x.mass):
                Element._elements_by_mass[elem.mass.value_in_unit(daltons)] = elem

        diff = mass
        best_guess = None

        for elemmass, element in _iteritems(Element._elements_by_mass):
            massdiff = abs(elemmass - mass)
            if massdiff < diff:
                best_guess = element
                diff = massdiff
            if elemmass > mass:
                # Elements are only getting heavier, so bail out early
                return best_guess

        # This really should only happen if we wanted ununoctium or something
        # bigger... won't really happen but still make sure we return an Element
        return best_guess

    @property
    def atomic_number(self):
        return self._atomic_number

    @property
    def name(self):
        return self._name

    @property
    def symbol(self):
        return self._symbol

    @property
    def mass(self):
        return self._mass

    def __str__(self):
        return '<Element %s>' % self.name

    def __repr__(self):
        return '<Element %s>' % self.name

# This is for backward compatibility.
def get_by_symbol(symbol):
    """ Get the element with a particular chemical symbol. """
    s = symbol.strip().upper()
    return Element._elements_by_symbol[s]

def _pickle_element(element):
    return (get_by_symbol, (element.symbol,))

copyreg.pickle(Element, _pickle_element)

# NOTE: getElementByMass assumes all masses are Quantity instances with unit
# "daltons". All elements need to obey this assumption, or that method will
# fail. No checking is done in getElementByMass for performance reasons
hydrogen =       Element(  1, "hydrogen", "H", 1.007947*daltons)
deuterium =      Element(  1, "deuterium", "D", 2.01355321270*daltons)
helium =         Element(  2, "helium", "He", 4.003*daltons)
lithium =        Element(  3, "lithium", "Li", 6.9412*daltons)
beryllium =      Element(  4, "beryllium", "Be", 9.0121823*daltons)
boron =          Element(  5, "boron", "B", 10.8117*daltons)
carbon =         Element(  6, "carbon", "C", 12.01078*daltons)
nitrogen =       Element(  7, "nitrogen", "N", 14.00672*daltons)
oxygen =         Element(  8, "oxygen", "O", 15.99943*daltons)
fluorine =       Element(  9, "fluorine", "F", 18.99840325*daltons)
neon =           Element( 10, "neon", "Ne", 20.17976*daltons)
sodium =         Element( 11, "sodium", "Na", 22.989769282*daltons)
magnesium =      Element( 12, "magnesium", "Mg", 24.30506*daltons)
aluminum =       Element( 13, "aluminum", "Al", 26.98153868*daltons)
silicon =        Element( 14, "silicon", "Si", 28.08553*daltons)
phosphorus =     Element( 15, "phosphorus", "P", 30.9737622*daltons)
sulfur =         Element( 16, "sulfur", "S", 32.0655*daltons)
chlorine =       Element( 17, "chlorine", "Cl", 35.4532*daltons)
argon =          Element( 18, "argon", "Ar", 39.9481*daltons)
potassium =      Element( 19, "potassium", "K", 39.09831*daltons)
calcium =        Element( 20, "calcium", "Ca", 40.0784*daltons)
scandium =       Element( 21, "scandium", "Sc", 44.9559126*daltons)
titanium =       Element( 22, "titanium", "Ti", 47.8671*daltons)
vanadium =       Element( 23, "vanadium", "V", 50.94151*daltons)
chromium =       Element( 24, "chromium", "Cr", 51.99616*daltons)
manganese =      Element( 25, "manganese", "Mn", 54.9380455*daltons)
iron =           Element( 26, "iron", "Fe", 55.8452*daltons)
cobalt =         Element( 27, "cobalt", "Co", 58.9331955*daltons)
nickel =         Element( 28, "nickel", "Ni", 58.69342*daltons)
copper =         Element( 29, "copper", "Cu", 63.5463*daltons)
zinc =           Element( 30, "zinc", "Zn", 65.4094*daltons)
gallium =        Element( 31, "gallium", "Ga", 69.7231*daltons)
germanium =      Element( 32, "germanium", "Ge", 72.641*daltons)
arsenic =        Element( 33, "arsenic", "As", 74.921602*daltons)
selenium =       Element( 34, "selenium", "Se", 78.963*daltons)
bromine =        Element( 35, "bromine", "Br", 79.9041*daltons)
krypton =        Element( 36, "krypton", "Kr", 83.7982*daltons)
rubidium =       Element( 37, "rubidium", "Rb", 85.46783*daltons)
strontium =      Element( 38, "strontium", "Sr", 87.621*daltons)
yttrium =        Element( 39, "yttrium", "Y", 88.905852*daltons)
zirconium =      Element( 40, "zirconium", "Zr", 91.2242*daltons)
niobium =        Element( 41, "niobium", "Nb", 92.906382*daltons)
molybdenum =     Element( 42, "molybdenum", "Mo", 95.942*daltons)
technetium =     Element( 43, "technetium", "Tc", 98*daltons)
ruthenium =      Element( 44, "ruthenium", "Ru", 101.072*daltons)
rhodium =        Element( 45, "rhodium", "Rh", 102.905502*daltons)
palladium =      Element( 46, "palladium", "Pd", 106.421*daltons)
silver =         Element( 47, "silver", "Ag", 107.86822*daltons)
cadmium =        Element( 48, "cadmium", "Cd", 112.4118*daltons)
indium =         Element( 49, "indium", "In", 114.8183*daltons)
tin =            Element( 50, "tin", "Sn", 118.7107*daltons)
antimony =       Element( 51, "antimony", "Sb", 121.7601*daltons)
tellurium =      Element( 52, "tellurium", "Te", 127.603*daltons)
iodine =         Element( 53, "iodine", "I", 126.904473*daltons)
xenon =          Element( 54, "xenon", "Xe", 131.2936*daltons)
cesium =         Element( 55, "cesium", "Cs", 132.90545192*daltons)
barium =         Element( 56, "barium", "Ba", 137.3277*daltons)
lanthanum =      Element( 57, "lanthanum", "La", 138.905477*daltons)
cerium =         Element( 58, "cerium", "Ce", 140.1161*daltons)
praseodymium =   Element( 59, "praseodymium", "Pr", 140.907652*daltons)
neodymium =      Element( 60, "neodymium", "Nd", 144.2423*daltons)
promethium =     Element( 61, "promethium", "Pm", 145*daltons)
samarium =       Element( 62, "samarium", "Sm", 150.362*daltons)
europium =       Element( 63, "europium", "Eu", 151.9641*daltons)
gadolinium =     Element( 64, "gadolinium", "Gd", 157.253*daltons)
terbium =        Element( 65, "terbium", "Tb", 158.925352*daltons)
dysprosium =     Element( 66, "dysprosium", "Dy", 162.5001*daltons)
holmium =        Element( 67, "holmium", "Ho", 164.930322*daltons)
erbium =         Element( 68, "erbium", "Er", 167.2593*daltons)
thulium =        Element( 69, "thulium", "Tm", 168.934212*daltons)
ytterbium =      Element( 70, "ytterbium", "Yb", 173.043*daltons)
lutetium =       Element( 71, "lutetium", "Lu", 174.9671*daltons)
hafnium =        Element( 72, "hafnium", "Hf", 178.492*daltons)
tantalum =       Element( 73, "tantalum", "Ta", 180.947882*daltons)
tungsten =       Element( 74, "tungsten", "W", 183.841*daltons)
rhenium =        Element( 75, "rhenium", "Re", 186.2071*daltons)
osmium =         Element( 76, "osmium", "Os", 190.233*daltons)
iridium =        Element( 77, "iridium", "Ir", 192.2173*daltons)
platinum =       Element( 78, "platinum", "Pt", 195.0849*daltons)
gold =           Element( 79, "gold", "Au", 196.9665694*daltons)
mercury =        Element( 80, "mercury", "Hg", 200.592*daltons)
thallium =       Element( 81, "thallium", "Tl", 204.38332*daltons)
lead =           Element( 82, "lead", "Pb", 207.21*daltons)
bismuth =        Element( 83, "bismuth", "Bi", 208.980401*daltons)
polonium =       Element( 84, "polonium", "Po", 209*daltons)
astatine =       Element( 85, "astatine", "At", 210*daltons)
radon =          Element( 86, "radon", "Rn", 222.018*daltons)
francium =       Element( 87, "francium", "Fr", 223*daltons)
radium =         Element( 88, "radium", "Ra", 226*daltons)
actinium =       Element( 89, "actinium", "Ac", 227*daltons)
thorium =        Element( 90, "thorium", "Th", 232.038062*daltons)
protactinium =   Element( 91, "protactinium", "Pa", 231.035882*daltons)
uranium =        Element( 92, "uranium", "U", 238.028913*daltons)
neptunium =      Element( 93, "neptunium", "Np", 237*daltons)
plutonium =      Element( 94, "plutonium", "Pu", 244*daltons)
americium =      Element( 95, "americium", "Am", 243*daltons)
curium =         Element( 96, "curium", "Cm", 247*daltons)
berkelium =      Element( 97, "berkelium", "Bk", 247*daltons)
californium =    Element( 98, "californium", "Cf", 251*daltons)
einsteinium =    Element( 99, "einsteinium", "Es", 252*daltons)
fermium =        Element(100, "fermium", "Fm", 257*daltons)
mendelevium =    Element(101, "mendelevium", "Md", 258*daltons)
nobelium =       Element(102, "nobelium", "No", 259*daltons)
lawrencium =     Element(103, "lawrencium",     "Lr", 262*daltons)
rutherfordium =  Element(104, "rutherfordium",  "Rf", 261*daltons)
dubnium =        Element(105, "dubnium",        "Db", 262*daltons)
seaborgium =     Element(106, "seaborgium",     "Sg", 266*daltons)
bohrium =        Element(107, "bohrium",        "Bh", 264*daltons)
hassium =        Element(108, "hassium",        "Hs", 269*daltons)
meitnerium =     Element(109, "meitnerium",     "Mt", 268*daltons)
darmstadtium =   Element(110, "darmstadtium",   "Ds", 281*daltons)
roentgenium =    Element(111, "roentgenium",    "Rg", 272*daltons)
ununbium =       Element(112, "ununbium",       "Uub", 285*daltons)
ununtrium =      Element(113, "ununtrium",      "Uut", 284*daltons)
ununquadium =    Element(114, "ununquadium",    "Uuq", 289*daltons)
ununpentium =    Element(115, "ununpentium",    "Uup", 288*daltons)
ununhexium =     Element(116, "ununhexium",     "Uuh", 292*daltons)

# Aliases to recognize common alternative spellings. Both the '==' and 'is'
# relational operators will work with any chosen name
sulphur = sulfur
aluminium = aluminum

if sys.version_info >= (3, 0):
    def _iteritems(dict):
        return dict.items()
else:
    def _iteritems(dict):
        return dict.iteritems()
