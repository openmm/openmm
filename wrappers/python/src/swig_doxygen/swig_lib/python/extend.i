%extend OpenMM::Context {
  PyObject *_getStateAsLists(int getPositions,
                            int getVelocities,
                            int getForces,
                            int getEnergy,
                            int getParameters,
                            int enforcePeriodic,
                            int groups) {
    int types = 0;
    if (getPositions) types |= State::Positions;
    if (getVelocities) types |= State::Velocities;
    if (getForces) types |= State::Forces;
    if (getEnergy) types |= State::Energy;
    if (getParameters) types |= State::Parameters;
    State state = self->getState(types, enforcePeriodic, groups);
    return _convertStateToLists(state);
  }


  %pythoncode {
    def getState(self,
                 getPositions=False,
                 getVelocities=False,
                 getForces=False,
                 getEnergy=False,
                 getParameters=False,
                 enforcePeriodicBox=False,
                 groups=-1):
        """
        getState(self,
                 getPositions = False,
                 getVelocities = False,
                 getForces = False,
                 getEnergy = False,
                 getParameters = False,
                 enforcePeriodicBox = False,
                 groups = -1)
              -> State
        Get a State object recording the current state information stored in this context.
           getPositions -- whether to store particle positions in the State
           getVelocities -- whether to store particle velocities in the State
           getForces -- whether to store the forces acting on particles in the State
           getEnergy -- whether to store potential and kinetic energy in the State
           getParameter -- whether to store context parameters in the State
           enforcePeriodicBox -- if false, the position of each particle will be whatever position is stored in the Context, regardless of periodic boundary conditions.  If true, particle positions will be translated so the center of every molecule lies in the same periodic box.
           groups -- a set of bit flags for which force groups to include when computing forces and energies.  Group i will be included if (groups&(1<<i)) != 0.  The default value includes all groups.
        """
        
        if getPositions: getP=1
        else: getP=0
        if getVelocities: getV=1
        else: getV=0
        if getForces: getF=1
        else: getF=0
        if getEnergy: getE=1
        else: getE=0
        if getParameters: getPa=1
        else: getPa=0
        if enforcePeriodicBox: enforcePeriodic=1
        else: enforcePeriodic=0

        (simTime, periodicBoxVectorsList, energy, coordList, velList,
         forceList, paramMap) = \
            self._getStateAsLists(getP, getV, getF, getE, getPa, enforcePeriodic, groups)
        
        state = State(simTime=simTime,
                      energy=energy,
                      coordList=coordList,
                      velList=velList,
                      forceList=forceList,
                      periodicBoxVectorsList=periodicBoxVectorsList,
                      paramMap=paramMap)
        return state
  }

}

%extend OpenMM::RPMDIntegrator {
  PyObject *_getStateAsLists(int copy,
                            int getPositions,
                            int getVelocities,
                            int getForces,
                            int getEnergy,
                            int getParameters,
                            int enforcePeriodic,
                            int groups) {
    int types = 0;
    if (getPositions) types |= State::Positions;
    if (getVelocities) types |= State::Velocities;
    if (getForces) types |= State::Forces;
    if (getEnergy) types |= State::Energy;
    if (getParameters) types |= State::Parameters;
    State state = self->getState(copy, types, enforcePeriodic, groups);
    return _convertStateToLists(state);
  }


  %pythoncode {
    def getState(self,
                 copy,
                 getPositions=False,
                 getVelocities=False,
                 getForces=False,
                 getEnergy=False,
                 getParameters=False,
                 enforcePeriodicBox=False,
                 groups=-1):
        """
        getState(self,
                 copy,
                 getPositions = False,
                 getVelocities = False,
                 getForces = False,
                 getEnergy = False,
                 getParameters = False,
                 enforcePeriodicBox = False,
                 groups = -1)
              -> State
        Get a State object recording the current state information about one copy of the system.
           copy -- the index of the copy for which to retrieve state information
           getPositions -- whether to store particle positions in the State
           getVelocities -- whether to store particle velocities in the State
           getForces -- whether to store the forces acting on particles in the State
           getEnergy -- whether to store potential and kinetic energy in the State
           getParameter -- whether to store context parameters in the State
           enforcePeriodicBox -- if false, the position of each particle will be whatever position is stored in the Context, regardless of periodic boundary conditions.  If true, particle positions will be translated so the center of every molecule lies in the same periodic box.
           groups -- a set of bit flags for which force groups to include when computing forces and energies.  Group i will be included if (groups&(1<<i)) != 0.  The default value includes all groups.
        """
        
        if getPositions: getP=1
        else: getP=0
        if getVelocities: getV=1
        else: getV=0
        if getForces: getF=1
        else: getF=0
        if getEnergy: getE=1
        else: getE=0
        if getParameters: getPa=1
        else: getPa=0
        if enforcePeriodicBox: enforcePeriodic=1
        else: enforcePeriodic=0

        (simTime, periodicBoxVectorsList, energy, coordList, velList,
         forceList, paramMap) = \
            self._getStateAsLists(copy, getP, getV, getF, getE, getPa, enforcePeriodic, groups)
        
        state = State(simTime=simTime,
                      energy=energy,
                      coordList=coordList,
                      velList=velList,
                      forceList=forceList,
                      periodicBoxVectorsList=periodicBoxVectorsList,
                      paramMap=paramMap)
        return state
  }

}


%extend OpenMM::NonbondedForce {
  %pythoncode {
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
  }
}


%extend OpenMM::System {
  %pythoncode {
    def __getstate__(self):
        serializationString = XmlSerializer.serializeSystem(self)
        return serializationString

    def __setstate__(self, serializationString):
        system = XmlSerializer.deserializeSystem(serializationString)
        self.this = system.this
  }
}

%extend OpenMM::CustomIntegrator {
    PyObject* getPerDofVariable(int index) const {
        std::vector<Vec3> values;
        self->getPerDofVariable(index, values);
        return copyVVec3ToList(values);
    }
}