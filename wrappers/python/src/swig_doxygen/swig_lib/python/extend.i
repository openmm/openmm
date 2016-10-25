%inline %{
    typedef int bitmask32t;
%}

%typemap(in) bitmask32t %{
    $1 = 0;
#if PY_VERSION_HEX >= 0x03000000
    if (PyLong_Check($input)) {
        unsigned long u = PyLong_AsUnsignedLongMask($input);
#else
    if (PyInt_Check($input)) {
        unsigned long u = PyInt_AsUnsignedLongMask($input);
#endif
        // 64-bit Windows has 32-bit longs, but other platforms have
        // 64-bit longs
        $1 = u & 0xffffffff;
    } else {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
%}



%extend OpenMM::Context {
  PyObject *_getStateAsLists(int getPositions,
                            int getVelocities,
                            int getForces,
                            int getEnergy,
                            int getParameters,
                            int getParameterDerivatives,
                            int enforcePeriodic,
                            bitmask32t groups) {
    State state;
    PyThreadState* _savePythonThreadState = PyEval_SaveThread();
    int types = 0;
    if (getPositions) types |= State::Positions;
    if (getVelocities) types |= State::Velocities;
    if (getForces) types |= State::Forces;
    if (getEnergy) types |= State::Energy;
    if (getParameters) types |= State::Parameters;
    if (getParameterDerivatives) types |= State::ParameterDerivatives;
    try {
        state = self->getState(types, enforcePeriodic, groups);
    }
    catch (...) {
        PyEval_RestoreThread(_savePythonThreadState);
        throw;
    }
    PyEval_RestoreThread(_savePythonThreadState);
    return _convertStateToLists(state);
  }


  %pythoncode %{
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
        getP, getV, getF, getE, getPa, getPd, enforcePeriodic = map(bool,
            (getPositions, getVelocities, getForces, getEnergy, getParameters,
             getParameterDerivatives, enforcePeriodicBox))

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

        (simTime, periodicBoxVectorsList, energy, coordList, velList,
         forceList, paramMap, paramDerivMap) = \
            self._getStateAsLists(getP, getV, getF, getE, getPa, getPd, enforcePeriodic, groups_mask)

        state = State(simTime=simTime,
                      energy=energy,
                      coordList=coordList,
                      velList=velList,
                      forceList=forceList,
                      periodicBoxVectorsList=periodicBoxVectorsList,
                      paramMap=paramMap,
                      paramDerivMap=paramDerivMap)
        return state

    def setState(self, state):
        """
        setState(Context self, State state)

        Copy information from a State object into this Context.  This restores the Context to
        approximately the same state it was in when the State was created.  If the State does not include
        a piece of information (e.g. positions or velocities), that aspect of the Context is
        left unchanged.

        Even when all possible information is included in the State, the effect of calling this method
        is still less complete than loadCheckpoint().  For example, it does not restore the internal
        states of random number generators.  On the other hand, it has the advantage of not being hardware
        specific.
        """
        self.setTime(state._simTime)
        self.setPeriodicBoxVectors(state._periodicBoxVectorsList[0], state._periodicBoxVectorsList[1], state._periodicBoxVectorsList[2])
        if state._coordList is not None:
             self.setPositions(state._coordList)
        if state._velList is not None:
             self.setVelocities(state._velList)
        if state._paramMap is not None:
             for param in state._paramMap:
                 self.setParameter(param, state._paramMap[param])
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
  PyObject *_getStateAsLists(int copy,
                            int getPositions,
                            int getVelocities,
                            int getForces,
                            int getEnergy,
                            int getParameters,
                            int getParameterDerivatives, 
                            int enforcePeriodic,
                            int groups) {
    State state;
    PyThreadState* _savePythonThreadState = PyEval_SaveThread();
    int types = 0;
    if (getPositions) types |= State::Positions;
    if (getVelocities) types |= State::Velocities;
    if (getForces) types |= State::Forces;
    if (getEnergy) types |= State::Energy;
    if (getParameters) types |= State::Parameters;
    if (getParameterDerivatives) types |= State::ParameterDerivatives;
    try {
        state = self->getState(copy, types, enforcePeriodic, groups);
    }
    catch (...) {
        PyEval_RestoreThread(_savePythonThreadState);
        throw;
    }
    PyEval_RestoreThread(_savePythonThreadState);
    return _convertStateToLists(state);
  }


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

        (simTime, periodicBoxVectorsList, energy, coordList, velList,
          forceList, paramMap, paramDerivMap) = \
            self._getStateAsLists(getP, getV, getF, getE, getPa, getPd, enforcePeriodic, groups_mask)

        state = State(simTime=simTime,
                      energy=energy,
                      coordList=coordList,
                      velList=velList,
                      forceList=forceList,
                      periodicBoxVectorsList=periodicBoxVectorsList,
                      paramMap=paramMap,
                      paramDerivMap=paramDerivMap)
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
  %feature(docstring, "This method exists only for backward compatibility. @deprecated Use serialize() instead.") serializeSystem;
  static std::string serializeSystem(const OpenMM::System* object) {
      std::stringstream ss;
      OpenMM::XmlSerializer::serialize<OpenMM::System>(object, "System", ss);
      return ss.str();
  }

  %feature(docstring, "This method exists only for backward compatibility. @deprecated Use deserialize() instead.") deserializeSystem;
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

  static std::string _serializeStateAsLists(
                                const std::vector<Vec3>& pos,
                                const std::vector<Vec3>& vel,
                                const std::vector<Vec3>& forces,
                                double kineticEnergy,
                                double potentialEnergy,
                                double time,
                                const std::vector<Vec3>& boxVectors,
                                const std::map<string, double>& params,
                                const std::map<string, double>& paramDerivs,
                                int types) {
    OpenMM::State myState =  _convertListsToState(pos,vel,forces,kineticEnergy,potentialEnergy,time,boxVectors,params,paramDerivs,types);
    std::stringstream buffer;
    OpenMM::XmlSerializer::serialize<OpenMM::State>(&myState, "State", buffer);
    return buffer.str();
  }

  static PyObject* _deserializeStringIntoLists(const std::string &stateAsString) {
    std::stringstream ss;
    ss << stateAsString;
    OpenMM::State* deserializedState = OpenMM::XmlSerializer::deserialize<OpenMM::State>(ss);
    PyObject* obj = _convertStateToLists(*deserializedState);
    delete deserializedState;
    return obj;
  }

  %pythoncode %{
    @staticmethod
    def _serializeState(pythonState):
      positions = []
      velocities = []
      forces = []
      kineticEnergy = 0.0
      potentialEnergy = 0.0
      params = {}
      paramDerivs = {}
      types = 0
      try:
        positions = pythonState.getPositions().value_in_unit(unit.nanometers)
        types |= 1
      except:
        pass
      try:
        velocities = pythonState.getVelocities().value_in_unit(unit.nanometers/unit.picoseconds)
        types |= 2
      except:
        pass
      try:
        forces = pythonState.getForces().value_in_unit(unit.kilojoules_per_mole/unit.nanometers)
        types |= 4
      except:
        pass
      try:
        kineticEnergy = pythonState.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
        potentialEnergy = pythonState.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        types |= 8
      except:
        pass
      try:
        params = pythonState.getParameters()
        types |= 16
      except:
        pass
      try:
        params = pythonState.getEnergyParameterDerivatives()
        types |= 32
      except:
        pass
      time = pythonState.getTime().value_in_unit(unit.picoseconds)
      boxVectors = pythonState.getPeriodicBoxVectors().value_in_unit(unit.nanometers)
      string = XmlSerializer._serializeStateAsLists(positions, velocities, forces, kineticEnergy, potentialEnergy, time, boxVectors, params, paramDerivs, types)
      return string

    @staticmethod
    def _deserializeState(pythonString):

      (simTime, periodicBoxVectorsList, energy, coordList, velList,
       forceList, paramMap, paramDerivMap) = XmlSerializer._deserializeStringIntoLists(pythonString)

      state = State(simTime=simTime,
                    energy=energy,
                    coordList=coordList,
                    velList=velList,
                    forceList=forceList,
                    periodicBoxVectorsList=periodicBoxVectorsList,
                    paramMap=paramMap,
                    paramDerivMap=paramDerivMap)
      return state

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
