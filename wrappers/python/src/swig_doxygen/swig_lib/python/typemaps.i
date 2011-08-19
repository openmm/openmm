
/* Convert python list of tuples to C++ std::vector of Vec3 objects */
%typemap(in) std::vector<Vec3>& (std::vector<OpenMM::Vec3> vVec) {
  // typemap -- %typemap(in) std::vector<Vec3>& (std::vector<OpenMM::Vec3> vVec)
  int i, pLength;
  double x, y, z;
  PyObject *o;
  PyObject *o1;

  pLength=(int)PySequence_Length($input);
  for (i=0; i<pLength; i++) {
    o=PySequence_GetItem($input, i);

    o1=PySequence_GetItem(o, 0);
    x=PyFloat_AsDouble(o1);
    Py_DECREF(o1);

    o1=PySequence_GetItem(o, 1);
    y=PyFloat_AsDouble(o1);
    Py_DECREF(o1);

    o1=PySequence_GetItem(o, 2);
    z=PyFloat_AsDouble(o1);
    Py_DECREF(o1);

    Py_DECREF(o);
    vVec.push_back( OpenMM::Vec3(x, y, z) );
  }
  $1 = &vVec;
}


/* Convert python tuple to C++ Vec3 object*/
%typemap(in) Vec3 {
  // typemap -- %typemap(in) Vec3
  double x, y, z;
  PyObject *o;

  o=PySequence_GetItem($input, 0);
  x=PyFloat_AsDouble(o);
  Py_DECREF(o);

  o=PySequence_GetItem($input, 1);
  y=PyFloat_AsDouble(o);
  Py_DECREF(o);

  o=PySequence_GetItem($input, 2);
  z=PyFloat_AsDouble(o);
  Py_DECREF(o);

  $1 = OpenMM::Vec3(x, y, z);
}
%typemap(in) const Vec3& (OpenMM::Vec3 myVec) {
  // typemap -- %typemap(in) Vec3
  double x, y, z;
  PyObject *o;

  o=PySequence_GetItem($input, 0);
  x=PyFloat_AsDouble(o);
  Py_DECREF(o);

  o=PySequence_GetItem($input, 1);
  y=PyFloat_AsDouble(o);
  Py_DECREF(o);

  o=PySequence_GetItem($input, 2);
  z=PyFloat_AsDouble(o);
  Py_DECREF(o);

  myVec = OpenMM::Vec3(x, y, z);
  $1 = &myVec;
}

/* Convert C++ (Vec3&, Vec3&, Vec3&) object to python tuple or tuples */
%typemap(argout) (Vec3& a, Vec3& b, Vec3& c) {
  // %typemap(argout) (Vec3& a, Vec3& b, Vec3& c)
  PyObject* mm = PyImport_AddModule("simtk.openmm");
  PyObject* vec3 = PyObject_GetAttrString(mm, "Vec3");
  PyObject* args1 = Py_BuildValue("(d,d,d)", (*$1)[0], (*$1)[1], (*$1)[2]);
  PyObject* args2 = Py_BuildValue("(d,d,d)", (*$2)[0], (*$2)[1], (*$2)[2]);
  PyObject* args3 = Py_BuildValue("(d,d,d)", (*$3)[0], (*$3)[1], (*$3)[2]);
  PyObject* pyVec1 = PyObject_CallObject(vec3, args1);
  PyObject* pyVec2 = PyObject_CallObject(vec3, args2);
  PyObject* pyVec3 = PyObject_CallObject(vec3, args3);
  Py_DECREF(args1);
  Py_DECREF(args2);
  Py_DECREF(args3);
  PyObject *o, *o2, *o3;
  o = Py_BuildValue("[N, N, N]", pyVec1, pyVec2, pyVec3);
  if ((!$result) || ($result == Py_None)) {
    $result = o;
  } else {
    if (!PyTuple_Check($result)) {
      PyObject *o2 = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($result, 0, o2);
    }
  o3 = PyTuple_New(1);
  PyTuple_SetItem(o3, 0, o);
  o2 = $result;
  $result = PySequence_Concat(o2, o3);
  Py_DECREF(o2);
  Py_DECREF(o3);
  }
}

%typemap(in, numinputs=0) (Vec3& a, Vec3& b, Vec3& c) (Vec3 tempA, Vec3 tempB, Vec3 tempC) {
  $1 = &tempA;
  $2 = &tempB;
  $3 = &tempC;
}

