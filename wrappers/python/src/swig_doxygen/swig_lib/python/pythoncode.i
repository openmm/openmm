%pythoncode %{

try:
    import numpy
except:
    pass

import math
RMIN_PER_SIGMA=math.pow(2, 1/6.0)
RVDW_PER_SIGMA=math.pow(2, 1/6.0)/2.0

import simtk.unit as unit


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
                 paramMap=None):
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
        return self._eK0 * unit.kilojoule_per_mole

    def getPotentialEnergy(self):
        """Get the total potential energy of the system with units.
           To remove the units, you can divide the return value by
           unit.kilojoule_per_mole or unit.kilocalorie_per_mole, etc.
           See the following for details:
           https://simtk.org/home/python_units
        """
        return self._eP0 * unit.kilojoule_per_mole

    def getParameters(self):
        """Get a map containing the values of all parameters.
        """
        return self._paramMap


# Strings can cause trouble
# as can any container that has infinite levels of containment
def _is_string(x):
     # step 1) String is always a container
     # and its contents are themselves containers.
     try:
         first_item = iter(x).next()
         inner_item = iter(first_item).next()
         if first_item == inner_item:
             return True
         else:
             return False
     except TypeError:
         return False
     except StopIteration:
         return False

def stripUnits(args):
    """
    getState(self, quantity) 
          -> value with *no* units

    Examples
    >>> import simtk

    >>> x = 5
    >>> print x
    5

    >>> x = stripUnits((5*simtk.unit.nanometer,))
    >>> x
    (5,)

    >>> arg1 = 5*simtk.unit.angstrom
    >>> x = stripUnits((arg1,))
    >>> x
    (0.5,)

    >>> arg1 = 5
    >>> x = stripUnits((arg1,))
    >>> x
    (5,)

    >>> arg1 = (1*simtk.unit.angstrom, 5*simtk.unit.angstrom)
    >>> x = stripUnits((arg1,))
    >>> x
    ((0.10000000000000001, 0.5),)

    >>> arg1 = (1*simtk.unit.angstrom,
    ...         5*simtk.unit.kilojoule_per_mole,
    ...         1*simtk.unit.kilocalorie_per_mole)
    >>> y = stripUnits((arg1,))
    >>> y
    ((0.10000000000000001, 5, 4.1840000000000002),)

    """
    newArgList=[]
    for arg in args:
        if unit.is_quantity(arg):
            # JDC: Ugly workaround for OpenMM using 'bar' for fundamental pressure unit.
            if arg.unit.is_compatible(unit.bar):
                arg = arg / unit.bar
            else:
                arg=arg.value_in_unit_system(unit.md_unit_system)                
            # JDC: End workaround.
            #arg=arg.value_in_unit_system(unit.md_unit_system)
        elif isinstance(arg, dict):
            newArg = {}
            for key in arg:
                newKey = key
                newValue = arg[key]
                if not _is_string(newKey):
                    newKey = stripUnits(newKey)
                if not _is_string(newValue):
                    newValue = stripUnits(newValue)
                newArg[newKey] = newValue
            arg = newArg
        elif not _is_string(arg):
            try:
                iter(arg)
                # Reclusively strip units from all quantities
                arg=stripUnits(arg)
            except TypeError:
                pass
        newArgList.append(arg)
    return tuple(newArgList)
%}

%pythonappend OpenMM::Context::Context %{
    self._system = args[0]
    self._integrator = args[1]
%}
