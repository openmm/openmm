%inline %{
    #include <cstring>
    #include <numpy/arrayobject.h>
%}

%extend OpenMM::Context {

  %pythoncode %{
    def getIntegrator(self):
        return self._integrator

    def getState(self, getPositions=False, getVelocities=False,
                 getForces=False, getEnergy=False, getParameters=False,
                 getParameterDerivatives=False, enforcePeriodicBox=False, groups=-1):
        """Get a State object recording the current state information stored in this context.

        Parameters
        ----------
        getPositions : bool=False
            whether to store particle positions in the State
        getVelocities : bool=False
            whether to store particle velocities in the State
        getForces : bool=False
            whether to store the forces acting on particles in the State
        getEnergy : bool=False
            whether to store potential and kinetic energy in the State
        getParameters : bool=False
            whether to store context parameters in the State
        getParameterDerivatives : bool=False
            whether to store parameter derivatives in the State
        enforcePeriodicBox : bool=False
            if false, the position of each particle will be whatever position
            is stored in the Context, regardless of periodic boundary conditions.
            If true, particle positions will be translated so the center of
            every molecule lies in the same periodic box.
        groups : set={0,1,2,...,31}
            a set of indices for which force groups to include when computing
            forces and energies. The default value includes all groups. groups
            can also be passed as an unsigned integer interpreted as a bitmask,
            in which case group i will be included if (groups&(1<<i)) != 0.
        """
        try:
            # is the input integer-like?
            groups_mask = int(groups)
        except TypeError:
            if isinstance(groups, set):
                # nope, okay, then it should be an set
                groups_mask = functools.reduce(operator.or_,
                        ((1<<x) & 0xffffffff for x in groups))
            else:
                raise TypeError('%s is neither an int nor set' % groups)
        if groups_mask >= 0x80000000:
            groups_mask -= 0x100000000
        types = 0
        if getPositions:
            types += State.Positions
        if getVelocities:
            types += State.Velocities
        if getForces:
            types += State.Forces
        if getEnergy:
            types += State.Energy
        if getParameters:
            types += State.Parameters
        if getParameterDerivatives:
            types += State.ParameterDerivatives
        state = _openmm.Context_getState(self, types, enforcePeriodicBox, groups_mask)
        return state

  %}

  %feature("docstring") createCheckpoint "Create a checkpoint recording the current state of the Context.
This should be treated as an opaque block of binary data.  See loadCheckpoint() for more details.

Returns: a string containing the checkpoint data
"
  std::string createCheckpoint() {
    std::stringstream stream(std::ios_base::in | std::ios_base::out | std::ios_base::binary);
    self->createCheckpoint(stream);
    return stream.str();
  }

  %feature ("docstring") loadCheckpoint "Load a checkpoint that was written by createCheckpoint().

A checkpoint contains not only publicly visible data such as the particle positions and
velocities, but also internal data such as the states of random number generators.  Ideally,
loading a checkpoint should restore the Context to an identical state to when it was written,
such that continuing the simulation will produce an identical trajectory.  This is not strictly
guaranteed to be true, however, and should not be relied on.  For most purposes, however, the
internal state should be close enough to be reasonably considered equivalent.

A checkpoint contains data that is highly specific to the Context from which it was created.
It depends on the details of the System, the Platform being used, and the hardware and software
of the computer it was created on.  If you try to load it on a computer with different hardware,
or for a System that is different in any way, loading is likely to fail.  Checkpoints created
with different versions of OpenMM are also often incompatible.  If a checkpoint cannot be loaded,
that is signaled by throwing an exception.

Parameters:
 - checkpoint (string) the checkpoint data to load
"
  void loadCheckpoint(std::string checkpoint) {
    std::stringstream stream(std::ios_base::in | std::ios_base::out | std::ios_base::binary);
    stream << checkpoint;
    self->loadCheckpoint(stream);
  }
}

%extend OpenMM::RPMDIntegrator {
  %pythoncode %{
    def getState(self,
                 copy,
                 getPositions=False,
                 getVelocities=False,
                 getForces=False,
                 getEnergy=False,
                 getParameters=False,
                 getParameterDerivatives=False,
                 enforcePeriodicBox=False,
                 groups=-1):
        """Get a State object recording the current state information about one copy of the system.

        Parameters
        ----------
        copy : int
            the index of the copy for which to retrieve state information
        getPositions : bool=False
            whether to store particle positions in the State
        getVelocities : bool=False
            whether to store particle velocities in the State
        getForces : bool=False
            whether to store the forces acting on particles in the State
        getEnergy : bool=False
            whether to store potential and kinetic energy in the State
        getParameters : bool=False
            whether to store context parameters in the State
        getParameterDerivatives : bool=False
            whether to store parameter derivatives in the State
        enforcePeriodicBox : bool=False
            if false, the position of each particle will be whatever position
            is stored in the Context, regardless of periodic boundary conditions.
            If true, particle positions will be translated so the center of
            every molecule lies in the same periodic box.
        groups : set={0,1,2,...,31}
            a set of indices for which force groups to include when computing
            forces and energies. The default value includes all groups. groups
            can also be passed as an unsigned integer interpreted as a bitmask,
            in which case group i will be included if (groups&(1<<i)) != 0.
        """
        getP, getV, getF, getE, getPa, getPd, enforcePeriodic = map(bool,
            (getPositions, getVelocities, getForces, getEnergy, getParameters,
             getParameterDerivatives, enforcePeriodicBox))

        try:
            # is the input integer-like?
            groups_mask = int(groups)
        except TypeError:
            if isinstance(groups, set):
                groups_mask = functools.reduce(operator.or_,
                        ((1<<x) & 0xffffffff for x in groups))
            else:
                raise TypeError('%s is neither an int nor set' % groups)
        if groups_mask >= 0x80000000:
            groups_mask -= 0x100000000
        types = 0
        if getPositions:
            types += State.Positions
        if getVelocities:
            types += State.Velocities
        if getForces:
            types += State.Forces
        if getEnergy:
            types += State.Energy
        if getParameters:
            types += State.Parameters
        if getParameterDerivatives:
            types += State.ParameterDerivatives
        state = _openmm.RPMDIntegrator_getState(self, copy, types, enforcePeriodicBox, groups_mask)
        return state
  %}
}

%extend OpenMM::NonbondedForce {
  %pythoncode %{
    def addParticle_usingRVdw(self, charge, rVDW, epsilon):
        """Add particle using elemetrary charge.  Rvdw and epsilon,
           which is consistent with AMBER parameter file usage.
           Note that the sum of the radii of the two interacting atoms is
           the minimum energy point in the Lennard Jones potential and
           is often called rMin.  The conversion from sigma follows:
           rVDW = 2^1/6 * sigma/2
        """
        return self.addParticle(charge, rVDW/RVDW_PER_SIGMA, epsilon)

    def addException_usingRMin(self, particle1, particle2,
                               chargeProd, rMin, epsilon):
        """Add interaction exception using the product of the two atoms'
           elementary charges, rMin and epsilon, which is standard for AMBER
           force fields.  Note that rMin is the minimum energy point in the
           Lennard Jones potential.  The conversion from sigma is:
           rMin = 2^1/6 * sigma.
        """
        return self.addException(particle1, particle2,
                                 chargeProd, rMin/RMIN_PER_SIGMA, epsilon)
  %}
}

%extend OpenMM::System {
  %pythoncode %{
    def __getstate__(self):
        serializationString = XmlSerializer.serializeSystem(self)
        return serializationString

    def __setstate__(self, serializationString):
        system = XmlSerializer.deserializeSystem(serializationString)
        self.this = system.this
    def __deepcopy__(self, memo):
        return self.__copy__()
    def getForces(self):
        """Get the list of Forces in this System"""
        return [self.getForce(i) for i in range(self.getNumForces())]
  %}
  %newobject __copy__;
  OpenMM::System* __copy__() {
      return OpenMM::XmlSerializer::clone<OpenMM::System>(*self);
  }
}

%extend OpenMM::XmlSerializer {
  %feature(docstring, "This method exists only for backward compatibility.\n@deprecated Use serialize() instead.") serializeSystem;
  static std::string serializeSystem(const OpenMM::System* object) {
      std::stringstream ss;
      OpenMM::XmlSerializer::serialize<OpenMM::System>(object, "System", ss);
      return ss.str();
  }

  %feature(docstring, "This method exists only for backward compatibility.\n@deprecated Use deserialize() instead.") deserializeSystem;
  %newobject deserializeSystem;
  static OpenMM::System* deserializeSystem(const char* inputString) {
      std::stringstream ss;
      ss << inputString;
      return OpenMM::XmlSerializer::deserialize<OpenMM::System>(ss);
  }

  static std::string _serializeForce(const OpenMM::Force* object) {
      std::stringstream ss;
      OpenMM::XmlSerializer::serialize<OpenMM::Force>(object, "Force", ss);
      return ss.str();
  }

  %newobject _deserializeForce;
  static OpenMM::Force* _deserializeForce(const char* inputString) {
      std::stringstream ss;
      ss << inputString;
      return OpenMM::XmlSerializer::deserialize<OpenMM::Force>(ss);
  }

  static std::string _serializeIntegrator(const OpenMM::Integrator* object) {
      std::stringstream ss;
      OpenMM::XmlSerializer::serialize<OpenMM::Integrator>(object, "Integrator", ss);
      return ss.str();
  }

  %newobject _deserializeIntegrator;
  static OpenMM::Integrator* _deserializeIntegrator(const char* inputString) {
      std::stringstream ss;
      ss << inputString;
      return OpenMM::XmlSerializer::deserialize<OpenMM::Integrator>(ss);
  }

  static std::string _serializeTabulatedFunction(const OpenMM::TabulatedFunction* object) {
      std::stringstream ss;
      OpenMM::XmlSerializer::serialize<OpenMM::TabulatedFunction>(object, "TabulatedFunction", ss);
      return ss.str();
  }

  %newobject _deserializeTabulatedFunction;
  static OpenMM::TabulatedFunction* _deserializeTabulatedFunction(const char* inputString) {
      std::stringstream ss;
      ss << inputString;
      return OpenMM::XmlSerializer::deserialize<OpenMM::TabulatedFunction>(ss);
  }

  static std::string _serializeState(const OpenMM::State* object) {
      std::stringstream ss;
      OpenMM::XmlSerializer::serialize<OpenMM::State>(object, "State", ss);
      return ss.str();
  }

  %newobject _deserializeState;
  static OpenMM::State* _deserializeState(const char* inputString) {
      std::stringstream ss;
      ss << inputString;
      return OpenMM::XmlSerializer::deserialize<OpenMM::State>(ss);
  }

  %pythoncode %{
    @staticmethod
    def serialize(object):
      """Serialize an object as XML."""
      if isinstance(object, System):
        return XmlSerializer.serializeSystem(object)
      elif isinstance(object, Force):
        return XmlSerializer._serializeForce(object)
      elif isinstance(object, Integrator):
        return XmlSerializer._serializeIntegrator(object)
      elif isinstance(object, State):
        return XmlSerializer._serializeState(object)
      elif isinstance(object, TabulatedFunction):
        return XmlSerializer._serializeTabulatedFunction(object)
      raise ValueError("Unsupported object type")

    @staticmethod
    def deserialize(inputString):
      """Reconstruct an object that has been serialized as XML."""
      import re
      match = re.search("<([^?]\S*)", inputString)
      if match is None:
        raise ValueError("Invalid input string")
      type = match.groups()[0]
      if type == "System":
        return XmlSerializer.deserializeSystem(inputString)
      if type == "Force":
        return XmlSerializer._deserializeForce(inputString)
      if type == "Integrator":
        return XmlSerializer._deserializeIntegrator(inputString)
      if type == "State":
        return XmlSerializer._deserializeState(inputString)
      if type == "TabulatedFunction":
        return XmlSerializer._deserializeTabulatedFunction(inputString)
      raise ValueError("Unsupported object type")
  %}
}

%extend OpenMM::CustomIntegrator {
    PyObject* getPerDofVariable(int index) const {
        std::vector<Vec3> values;
        self->getPerDofVariable(index, values);
        return copyVVec3ToList(values);
    }
}

%extend OpenMM::Force {
  %pythoncode %{
    def __getstate__(self):
        serializationString = XmlSerializer.serialize(self)
        return serializationString

    def __setstate__(self, serializationString):
        system = XmlSerializer.deserialize(serializationString)
        self.this = system.this

    def __deepcopy__(self, memo):
        return self.__copy__()
  %}
  %newobject __copy__;
  OpenMM::Force* __copy__() {
      return OpenMM::XmlSerializer::clone<OpenMM::Force>(*self);
  }
}

%extend OpenMM::Integrator {
  %pythoncode %{
    def __getstate__(self):
        serializationString = XmlSerializer.serialize(self)
        return serializationString

    def __setstate__(self, serializationString):
        system = XmlSerializer.deserialize(serializationString)
        self.this = system.this

    def __deepcopy__(self, memo):
        return self.__copy__()
  %}
  %newobject __copy__;
  OpenMM::Integrator* __copy__() {
      return OpenMM::XmlSerializer::clone<OpenMM::Integrator>(*self);
  }
}

%extend OpenMM::TabulatedFunction {
  %pythoncode %{
    def __getstate__(self):
        serializationString = XmlSerializer.serialize(self)
        return serializationString

    def __setstate__(self, serializationString):
        system = XmlSerializer.deserialize(serializationString)
        self.this = system.this

    def __deepcopy__(self, memo):
        return self.__copy__()
  %}
  %newobject __copy__;
  OpenMM::TabulatedFunction* __copy__() {
      return OpenMM::XmlSerializer::clone<OpenMM::TabulatedFunction>(*self);
  }
}

%extend OpenMM::State {
  %pythoncode %{
    def __getstate__(self):
        serializationString = XmlSerializer.serialize(self)
        return serializationString

    def __setstate__(self, serializationString):
        system = XmlSerializer.deserialize(serializationString)
        self.this = system.this

    def __deepcopy__(self, memo):
        return self.__copy__()

    def getPeriodicBoxVectors(self, asNumpy=False):
        """Get the vectors defining the axes of the periodic box."""
        vectors = _openmm.State_getPeriodicBoxVectors(self)
        if asNumpy:
            vectors = numpy.array(vectors)
        return vectors*unit.nanometers

    def getPositions(self, asNumpy=False):
        """Get the position of each particle with units.
           Raises an exception if positions where not requested in
           the context.getState() call.
           Returns a list of Vec3s, unless asNumpy is True, in
           which  case a Numpy array of arrays will be returned.
           """
        if asNumpy:
            if '_positionsNumpy' not in dir(self):
                self._positionsNumpy = numpy.empty([self._getNumParticles(), 3], numpy.float64)
                self._getVectorAsNumpy(State.Positions, self._positionsNumpy)
                self._positionsNumpy = self._positionsNumpy*unit.nanometers
            return self._positionsNumpy
        if '_positions' not in dir(self):
            self._positions = self._getVectorAsVec3(State.Positions)*unit.nanometers
        return self._positions

    def getVelocities(self, asNumpy=False):
        """Get the velocity of each particle with units.
           Raises an exception if velocities where not requested in
           the context.getState() call.
           Returns a list of Vec3s if asNumpy is False, or a Numpy
           array if asNumpy is True.
           """
        if asNumpy:
            if '_velocitiesNumpy' not in dir(self):
                self._velocitiesNumpy = numpy.empty([self._getNumParticles(), 3], numpy.float64)
                self._getVectorAsNumpy(State.Velocities, self._velocitiesNumpy)
                self._velocitiesNumpy = self._velocitiesNumpy*unit.nanometers/unit.picosecond
            return self._velocitiesNumpy
        if '_velocities' not in dir(self):
            self._velocities = self._getVectorAsVec3(State.Velocities)*unit.nanometers/unit.picosecond
        return self._velocities

    def getForces(self, asNumpy=False):
        """Get the force acting on each particle with units.
           Raises an exception if forces where not requested in
           the context.getState() call.
           Returns a list of Vec3s if asNumpy is False, or a Numpy
           array if asNumpy is True.
           """
        if asNumpy:
            if '_forcesNumpy' not in dir(self):
                self._forcesNumpy = numpy.empty([self._getNumParticles(), 3], numpy.float64)
                self._getVectorAsNumpy(State.Forces, self._forcesNumpy)
                self._forcesNumpy = self._forcesNumpy*unit.kilojoules_per_mole/unit.nanometer
            return self._forcesNumpy
        if '_forces' not in dir(self):
            self._forces = self._getVectorAsVec3(State.Forces)*unit.kilojoules_per_mole/unit.nanometer
        return self._forces
  %}
  
  int _getNumParticles() {
      if ((self->getDataTypes() & State::Positions) != 0)
          return self->getPositions().size();
      if ((self->getDataTypes() & State::Velocities) != 0)
          return self->getVelocities().size();
      if ((self->getDataTypes() & State::Forces) != 0)
          return self->getForces().size();
      return 0;
  }
  
  PyObject* _getVectorAsVec3(State::DataType type) {
      if (type == State::Positions)
          return copyVVec3ToList(self->getPositions());
      if (type == State::Velocities)
          return copyVVec3ToList(self->getVelocities());
      if (type == State::Forces)
          return copyVVec3ToList(self->getForces());
      PyErr_SetString(PyExc_ValueError, "Illegal type specified in _getVectorAsVec3");
      return NULL;
  }
  
  void _getVectorAsNumpy(State::DataType type, PyObject* output) {
      const std::vector<Vec3>* array;
      if (type == State::Positions)
          array = &self->getPositions();
      else if (type == State::Velocities)
          array = &self->getVelocities();
      else if (type == State::Forces)
          array = &self->getForces();
      else {
        PyErr_SetString(PyExc_ValueError, "Illegal type specified in _getVectorAsNumpy");
        return;
      }
      void* data = PyArray_DATA((PyArrayObject*) output);
      memcpy(data, &array[0][0], 3*sizeof(double)*array->size());
  }

  %newobject __copy__;
  OpenMM::State* __copy__() {
      return OpenMM::XmlSerializer::clone<OpenMM::State>(*self);
  }
}
