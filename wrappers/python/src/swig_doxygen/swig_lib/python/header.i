
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

int isNumpyAvailable() {
    static bool initialized = false;
    static bool available = false;
    if (!initialized) {
        initialized = true;
        available = (_import_array() >= 0);
    }
    return available;
}

} // namespace OpenMM
%}


