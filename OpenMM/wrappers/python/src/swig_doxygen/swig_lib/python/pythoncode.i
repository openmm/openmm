%pythoncode %{

try:
    import numpy
except ImportError:
    numpy = None

import copy
import sys
import math
import functools
import operator
RMIN_PER_SIGMA=math.pow(2, 1/6.0)
RVDW_PER_SIGMA=math.pow(2, 1/6.0)/2.0
if sys.version_info[0] == 2:
    _string_types = (basestring,)
else:
    _string_types = (bytes, str)

import simtk.unit as unit
from simtk.openmm.vec3 import Vec3

class State(_object):
    """
     A State object records a snapshot of the
     current state of a simulation at a point
     in time.  You create it by calling
     getState() on a Context.

     When a State is created, you specify what
     information should be stored in it.  This
     saves time and memory by only copying in
     the information that you actually want.
     This is especially important for forces
     and energies, since they may need to be
     calculated.  If you query a State object
     for a piece of information which is not
     available (because it was not requested
     when the State was created), it will
     return None.

     In general return values are Python Units
     (https://simtk.org/home/python_units).
     Among other things Python Units provides a
     container class, Quantity, which holds a
     value and a representation of the value's
     unit.  Values can be integers, floats,
     lists, numarrays, etc.  Quantity objects
     can be used in arithmetic operation just
     like number, except they also keep track
     of units.   To extract the value from a
     quantity, us the value_in_unit() method.
     For example, to extract the value from a
     length quantity, in units of nanometers,
     do the following:
     myLengthQuantity.value_in_unit(unit.nanometer)

"""
    def __init__(self,
                 simTime=None,
                 energy=None,
                 coordList=None,
                 velList=None,
                 forceList=None,
                 periodicBoxVectorsList=None,
                 paramMap=None,
                 paramDerivMap=None):
        self._simTime=simTime
        self._periodicBoxVectorsList=periodicBoxVectorsList
        self._periodicBoxVectorsListNumpy=None
        if energy:
            self._eK0=energy[0]
            self._eP0=energy[1]
        else:
            self._eK0=None
            self._eP0=None
        self._coordList=coordList
        self._coordListNumpy=None
        self._velList=velList
        self._velListNumpy=None
        self._forceList=forceList
        self._forceListNumpy=None
        self._paramMap=paramMap
        self._paramDerivMap=paramDerivMap

    def __getstate__(self):
        serializationString = XmlSerializer.serialize(self)
        return serializationString

    def __setstate__(self, serializationString):
        dState = XmlSerializer.deserialize(serializationString)
        # Safe provided no __slots__ or other weird things are used
        self.__dict__.update(dState.__dict__)

    def getTime(self):
        """Get the time for which this State was created."""
        return self._simTime * unit.picosecond

    def getPeriodicBoxVectors(self, asNumpy=False):
        """Get the three periodic box vectors if this state is from a
           simulation using PBC ."""
        if self._periodicBoxVectorsList is None:
            raise TypeError('periodic box vectors were not available.')

        if asNumpy:
            if self._periodicBoxVectorsListNumpy is None:
                self._periodicBoxVectorsListNumpy = \
                     numpy.array(self._periodicBoxVectorsList)
            returnValue=self._periodicBoxVectorsListNumpy
        else:
            returnValue=self._periodicBoxVectorsList

        returnValue = unit.Quantity(returnValue, unit.nanometers)
        return returnValue

    def getPeriodicBoxVolume(self):
        """Get the volume of the periodic box."""
        a = self._periodicBoxVectorsList[0]
        b = self._periodicBoxVectorsList[1]
        c = self._periodicBoxVectorsList[2]
        bcrossc = Vec3(b[1]*c[2]-b[2]*c[1], b[2]*c[0]-b[0]*c[2], b[0]*c[1]-b[1]*c[0])
        return unit.Quantity(unit.dot(a, bcrossc), unit.nanometers*unit.nanometers*unit.nanometers)

    def getPositions(self, asNumpy=False):
        """Get the position of each particle with units.
           Raises an exception if postions where not requested in
           the context.getState() call.
           Returns a list of tuples, unless asNumpy is True, in
           which  case a Numpy array of arrays will be returned.
           To remove the units, divide return value by unit.angstrom
           or unit.nanometer.  See the following for details:
           https://simtk.org/home/python_units
           """
        if self._coordList is None:
            raise TypeError('Positions were not requested in getState() call, so are not available.')

        if asNumpy:
            if self._coordListNumpy is None:
                self._coordListNumpy=numpy.array(self._coordList)
            returnValue=self._coordListNumpy
        else:
            returnValue=self._coordList

        returnValue = unit.Quantity(returnValue, unit.nanometers)
        return returnValue

    def getVelocities(self, asNumpy=False):
        """Get the velocity of each particle with units.
           Raises an exception if velocities where not requested in
           the context.getState() call.
           Returns a list of tuples, unless asNumpy is True, in
           which  case a Numpy array of arrays will be returned.
           To remove the units, you can divide the return value by
           unit.angstrom/unit.picosecond or unit.meter/unit.second,
           etc.  See the following for details:
           https://simtk.org/home/python_units
           """
        if self._velList is None:
            raise TypeError('Velocities were not requested in getState() call, so are not available.')

        if asNumpy:
            if self._velListNumpy is None:
                self._velListNumpy=numpy.array(self._velList)
            returnValue=self._velListNumpy
        else:
            returnValue=self._velList

        returnValue = unit.Quantity(returnValue, unit.nanometers/unit.picosecond)
        return returnValue

    def getForces(self, asNumpy=False):
        """Get the force acting on each particle with units.
           Raises an exception if forces where not requested in
           the context.getState() call.
           Returns a list of tuples, unless asNumpy is True, in
           which  case a Numpy array of arrays will be returned.
           To remove the units, you can divide the return value by
           unit.kilojoule_per_mole/unit.angstrom or
           unit.calorie_per_mole/unit.nanometer, etc.
           See the following for details:
           https://simtk.org/home/python_units
           """
        if self._forceList is None:
            raise TypeError('Forces were not requested in getState() call, so are not available.')

        if asNumpy:
            if self._forceListNumpy is None:
                self._forceListNumpy=numpy.array(self._forceList)
            returnValue=self._forceListNumpy
        else:
            returnValue=self._forceList

        returnValue = unit.Quantity(returnValue,
                                    unit.kilojoule_per_mole/unit.nanometer)
        return returnValue

    def getKineticEnergy(self):
        """Get the total kinetic energy of the system with units.
           To remove the units, you can divide the return value by
           unit.kilojoule_per_mole or unit.calorie_per_mole, etc.
           See the following for details:
           https://simtk.org/home/python_units
        """
        if self._eK0 is None:
            raise TypeError('Energy was not requested in getState() call, so it is not available.')
        return self._eK0 * unit.kilojoule_per_mole

    def getPotentialEnergy(self):
        """Get the total potential energy of the system with units.
           To remove the units, you can divide the return value by
           unit.kilojoule_per_mole or unit.kilocalorie_per_mole, etc.
           See the following for details:
           https://simtk.org/home/python_units
        """
        if self._eP0 is None:
            raise TypeError('Energy was not requested in getState() call, so it is not available.')
        return self._eP0 * unit.kilojoule_per_mole

    def getParameters(self):
        """Get a map containing the values of all parameters.
        """
        if self._paramMap is None:
            raise TypeError('Parameters were not requested in getState() call, so are not available.')
        return self._paramMap

    def getEnergyParameterDerivatives(self):
        """Get a map containing derivatives of the potential energy with respect to context parameters.

        In most cases derivatives are only calculated if the corresponding Force objects have been
        specifically told to compute them.  Otherwise, the values in the map will be zero.  Likewise,
        if multiple Forces depend on the same parameter but only some have been told to compute
        derivatives with respect to it, the returned value will include only the contributions from
        the Forces that were told to compute it."""
        if self._paramDerivMap is None:
            raise TypeError('Parameter derivatives were not requested in getState() call, so are not available.')
        return self._paramDerivMap

%}

%pythonappend OpenMM::Context::Context %{
    self._system = args[0]
    self._integrator = args[1]
%}

%pythonprepend OpenMM::AmoebaAngleForce::addAngle %{
    try:
        length = args[3]
        if isinstance(args, tuple):
            args = list(args)
    except (NameError, UnboundLocalError):
        if unit.is_quantity(length):
            length = length.value_in_unit(unit.degree)
    else:
        if unit.is_quantity(length):
            args[3] = length.value_in_unit(unit.degree)
%}

%pythonprepend OpenMM::AmoebaAngleForce::setAngleParameters %{
    try:
        length = args[4]
        if isinstance(args, tuple):
            args = list(args)
    except (NameError, UnboundLocalError):
        if unit.is_quantity(length):
            length = length.value_in_unit(unit.degree)
    else:
        if unit.is_quantity(length):
            args[4] = length.value_in_unit(unit.degree)
%}

%pythonprepend OpenMM::AmoebaTorsionTorsionForce::setTorsionTorsionGrid %{
    def deunitize_grid(grid):
        if isinstance(grid, tuple):
            grid = list(grid)
        for i, row in enumerate(grid):
            if isinstance(row, tuple):
                row = list(row)
                grid[i] = row
            for i, column in enumerate(row):
                if isinstance(column, tuple):
                    column = list(column)
                    row[i] = column
                # Data is angle, angle, energy, de/dang1, de/dang2, d^2e/dang1dang2
                if unit.is_quantity(column[0]):
                    column[0] = column[0].value_in_unit(unit.degree)
                if unit.is_quantity(column[1]):
                    column[1] = column[1].value_in_unit(unit.degree)
                if unit.is_quantity(column[2]):
                    column[2] = column[2].value_in_unit(unit.kilojoule_per_mole)
                if len(column) > 3 and unit.is_quantity(column[3]):
                    column[3] = column[3].value_in_unit(unit.kilojoule_per_mole/unit.radians)
                if len(column) > 4 and unit.is_quantity(column[4]):
                    column[4] = column[4].value_in_unit(unit.kilojoule_per_mole/unit.radians)
                if len(column) > 5 and unit.is_quantity(column[5]):
                    column[5] = column[5].value_in_unit(unit.kilojoule_per_mole/unit.radians**2)
        return grid
    try:
        grid = copy.deepcopy(args[1])
        if isinstance(args, tuple):
            args = list(args)
    except (NameError, UnboundLocalError):
        try:
            # Support numpy arrays
            grid = grid.tolist()
        except AttributeError:
            grid = copy.deepcopy(grid)
        grid = deunitize_grid(grid)
    else:
        args[1] = deunitize_grid(grid)
%}
