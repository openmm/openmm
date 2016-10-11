
%header %{
namespace OpenMM {

PyObject *copyVVec3ToList(std::vector<Vec3> vVec3) {
  int i, n;
  PyObject *pyList;

  n=vVec3.size(); 
  pyList=PyList_New(n);
  PyObject* mm = PyImport_AddModule("simtk.openmm");
  PyObject* vec3 = PyObject_GetAttrString(mm, "Vec3");
  for (i=0; i<n; i++) {
    OpenMM::Vec3& v = vVec3.at(i);
    PyObject* args = Py_BuildValue("(d,d,d)", v[0], v[1], v[2]);
    PyObject* pyVec = PyObject_CallObject(vec3, args);
    Py_DECREF(args);
    PyList_SET_ITEM(pyList, i, pyVec);
  }
  return pyList;
}

State _convertListsToState( const std::vector<Vec3> &pos, 
                            const std::vector<Vec3> &vel, 
                            const std::vector<Vec3> &forces,
                            double kineticEnergy,
                            double potentialEnergy,
                            double time,
                            const std::vector<Vec3> &boxVectors,
                            const std::map<std::string, double> &params,
                            const std::map<std::string, double> &paramDerivs,
                            int types ) {  
    State::StateBuilder sb(time); 
    if(types & State::Positions)
      sb.setPositions(pos);
    if(types & State::Velocities)
      sb.setVelocities(vel);
    if(types & State::Forces)
      sb.setForces(forces);
    if(types & State::Energy)
      sb.setEnergy(kineticEnergy, potentialEnergy);
    if(types & State::Parameters)
      sb.setParameters(params);
    if(types & State::ParameterDerivatives)
      sb.setEnergyParameterDerivatives(paramDerivs);
    sb.setPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    return sb.getState();
}

PyObject *_convertStateToLists(const State& state) {
    double simTime;
    PyObject *pPeriodicBoxVectorsList;
    PyObject *pEnergy;
    PyObject *pPositions;
    PyObject *pVelocities;
    PyObject *pForces;
    PyObject *pyTuple;
    PyObject *pParameters;
    PyObject *pParameterDerivs;
    simTime=state.getTime();

    OpenMM::Vec3 myVecA;
    OpenMM::Vec3 myVecB;
    OpenMM::Vec3 myVecC;
    state.getPeriodicBoxVectors(myVecA, myVecB, myVecC);
    PyObject* mm = PyImport_AddModule("simtk.openmm");
    PyObject* vec3 = PyObject_GetAttrString(mm, "Vec3");
    PyObject* args1 = Py_BuildValue("(d,d,d)", myVecA[0], myVecA[1], myVecA[2]);
    PyObject* args2 = Py_BuildValue("(d,d,d)", myVecB[0], myVecB[1], myVecB[2]);
    PyObject* args3 = Py_BuildValue("(d,d,d)", myVecC[0], myVecC[1], myVecC[2]);
    PyObject* pyVec1 = PyObject_CallObject(vec3, args1);
    PyObject* pyVec2 = PyObject_CallObject(vec3, args2);
    PyObject* pyVec3 = PyObject_CallObject(vec3, args3);
    Py_DECREF(args1);
    Py_DECREF(args2);
    Py_DECREF(args3);
    pPeriodicBoxVectorsList = Py_BuildValue("N,N,N", pyVec1, pyVec2, pyVec3);

    try {
      pPositions = copyVVec3ToList(state.getPositions());
    }
    catch (std::exception& ex) {
      pPositions = Py_None;
      Py_INCREF(Py_None);
    }
    try {
      pVelocities = copyVVec3ToList(state.getVelocities());
    }
    catch (std::exception& ex) {
      pVelocities = Py_None;
      Py_INCREF(Py_None);
    }
    try {
      pForces = copyVVec3ToList(state.getForces());
    }
    catch (std::exception& ex) {
      pForces = Py_None;
      Py_INCREF(Py_None);
    }
    try {
      pEnergy = Py_BuildValue("(d,d)",
                             state.getKineticEnergy(),
                             state.getPotentialEnergy());
    }
    catch (std::exception& ex) {
      pEnergy = Py_None;
      Py_INCREF(Py_None);
    }
    try {
      pParameters = PyDict_New();
      const std::map<std::string, double>& params = state.getParameters();
      for (std::map<std::string, double>::const_iterator iter = params.begin(); iter != params.end(); ++iter)
          PyDict_SetItemString(pParameters, iter->first.c_str(), Py_BuildValue("d", iter->second));
    }
    catch (std::exception& ex) {
      pParameters = Py_None;
      Py_INCREF(Py_None);
    }
    try {
      pParameterDerivs = PyDict_New();
      const std::map<std::string, double>& params = state.getEnergyParameterDerivatives();
      for (std::map<std::string, double>::const_iterator iter = params.begin(); iter != params.end(); ++iter)
          PyDict_SetItemString(pParameterDerivs, iter->first.c_str(), Py_BuildValue("d", iter->second));
    }
    catch (std::exception& ex) {
      pParameterDerivs = Py_None;
      Py_INCREF(Py_None);
    }
  
    pyTuple=Py_BuildValue("(d,N,N,N,N,N,N,N)",
                          simTime, pPeriodicBoxVectorsList, pEnergy,
                          pPositions, pVelocities,
                          pForces, pParameters, pParameterDerivs);
  
    return pyTuple;
}

} // namespace OpenMM
%}


