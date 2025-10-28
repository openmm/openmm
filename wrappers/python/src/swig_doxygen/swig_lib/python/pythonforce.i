%newobject new_PythonForce;
%newobject OpenMM::PythonForce::PythonForce;

%inline %{

namespace OpenMM {
    /**
     * This is the PythonForceComputation that performs the computation for a PythonForce.  It invokes the function
     * provided by the user, validates the outputs, and converts them to the required format.
     */
    class ComputationWrapper : public PythonForceComputation {
    public:
        ComputationWrapper(PyObject* computation) : computation(computation) {
        }
        void compute(const State& state, double& energy, std::vector<Vec3>& forces) const {
            PyGILState_STATE gstate;
            gstate = PyGILState_Ensure();

            // Invoke the function.

            swig_type_info* info = SWIGTYPE_p_OpenMM__State;
            PyObject* wrappedState = SWIG_NewPointerObj((void*) &state, info, 0);
            PyObject* result = PyObject_CallFunctionObjArgs(computation, wrappedState, NULL);

            // Extract the return values.

            if (!PyTuple_Check(result) || PyTuple_Size(result) != 2) {
                PyGILState_Release(gstate);
                throw OpenMMException("PythonForce: Expected two return values");
            }
            PyObject* pyenergy = Py_StripOpenMMUnits(PyTuple_GetItem(result, 0));
            PyObject* pyforces = Py_StripOpenMMUnits(PyTuple_GetItem(result, 1));
            energy = PyFloat_AsDouble(pyenergy);

            // Copy the forces to the output vector.

            if (!PyArray_Check(pyforces) || PyArray_NDIM((PyArrayObject*) pyforces) != 2) {
                PyGILState_Release(gstate);
                throw OpenMMException("PythonForce: The forces must be returned in a 2-dimensional NumPy array");
            }
            double** forceData;
            npy_intp dims[2];
            PyArray_Descr* descr = PyArray_DescrFromType(NPY_DOUBLE);
            if (PyArray_AsCArray(&pyforces, (void **) &forceData, dims, 2, descr) < 0) {
                PyGILState_Release(gstate);
                throw OpenMMException("PythonForce: Error accessing force array data");
            }
            if (dims[0] != forces.size() || dims[1] != 3) {
                PyGILState_Release(gstate);
                throw OpenMMException("PythonForce: The forces must be returned in a NumPy array of shape (# particles, 3)");
            }
            for (int i = 0; i < dims[0]; i++)
                for (int j = 0; j < 3; j++)
                    forces[i][j] = forceData[i][j];
            PyGILState_Release(gstate);
        }
    private:
        PyObject* computation;
    };
}

%}

%extend OpenMM::PythonForce {
    PythonForce(PyObject* computation, const std::map<std::string, double>& globalParameters={}) {
        return new PythonForce(new ComputationWrapper(computation), globalParameters);
    }
}
