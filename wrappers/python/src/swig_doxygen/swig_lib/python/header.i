
%header %{
namespace OpenMM {

PyObject *copyVVec3ToList(std::vector<Vec3> vVec3) {
  int i, n;
  PyObject *pyTuple;
  PyObject *pyList;
  OpenMM::Vec3 vec3;

  n=vVec3.size(); 
  pyList=PyList_New(n);
  for (i=0; i<n; i++) {
    vec3=vVec3.at(i);
    pyTuple=Py_BuildValue("(d,d,d)",
                          vec3[0], vec3[1], vec3[2]);
    PyList_SET_ITEM(pyList, i, pyTuple);
  }
  return pyList;
}


} // namespace OpenMM
%}


