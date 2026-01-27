%newobject OpenMM::PythonForce::PythonForce;

%inline %{

#include <iomanip>

namespace OpenMM {
    /**
     * This is the PythonForceComputation that performs the computation for a PythonForce.  It invokes the function
     * provided by the user, validates the outputs, and converts them to the required format.
     */
    class ComputationWrapper : public PythonForceComputation {
    public:
        ComputationWrapper(PyObject* computation) : computation(computation), batchComputation(nullptr) {
            Py_INCREF(computation);
        }
        ComputationWrapper(PyObject* computation, PyObject* batchComputation) : computation(computation), batchComputation(batchComputation) {
            Py_INCREF(computation);
            if (batchComputation != nullptr && batchComputation != Py_None)
                Py_INCREF(batchComputation);
            else
                this->batchComputation = nullptr;
        }
        ~ComputationWrapper() {
            Py_XDECREF(computation);
            if (batchComputation != nullptr)
                Py_XDECREF(batchComputation);
        }
        bool supportsBatchedEvaluation() const override {
            return batchComputation != nullptr;
        }
        void compute(const State& state, double& energy, void* forces, bool forcesAreDouble) const {
            PyGILState_STATE gstate;
            gstate = PyGILState_Ensure();

            // Invoke the function.

            swig_type_info* info = SWIGTYPE_p_OpenMM__State;
            PyObject* wrappedState = SWIG_NewPointerObj((void*) &state, info, 0);
            PyObject* result = PyObject_CallFunctionObjArgs(computation, wrappedState, NULL);
            if (result == NULL) {
                // The function raised an exception.  Convert it to an OpenMMException.

#if PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 12
                PyObject *type;
                PyObject *exception;
                PyObject *traceback;
                PyErr_Fetch(&type, &exception, &traceback);
#else
                PyObject *exception = PyErr_GetRaisedException();
#endif
                PyObject *message = PyObject_Str(exception);
                std::string *ptr;
                SWIG_AsPtr_std_string(message, &ptr);
                Py_XDECREF(message);
                PyGILState_Release(gstate);
                throw OpenMMException(*ptr);
            }

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
            npy_intp* dims = PyArray_DIMS((PyArrayObject*) pyforces);
            int numParticles = state.getPositions().size();
            if (dims[0] != numParticles || dims[1] != 3) {
                PyGILState_Release(gstate);
                throw OpenMMException("PythonForce: The forces must be returned in a NumPy array of shape (# particles, 3)");
            }
            PyObject* array;
            int targetType = (forcesAreDouble ? NPY_DOUBLE : NPY_FLOAT);
            if (PyArray_CHKFLAGS((PyArrayObject*) pyforces, NPY_ARRAY_CARRAY_RO) && PyArray_DESCR((PyArrayObject*) pyforces)->type_num == targetType)
                array = pyforces;
            else
                array = PyArray_FromAny(pyforces, PyArray_DescrFromType(targetType), 2, 2, NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_FORCECAST, NULL);
            int elementSize = (forcesAreDouble ? sizeof(double) : sizeof(float));
            void* data = PyArray_DATA((PyArrayObject*) array);
            memcpy(forces, data, 3*elementSize*numParticles);
            Py_XDECREF(wrappedState);
            Py_XDECREF(result);
            Py_XDECREF(pyenergy);
            Py_XDECREF(pyforces);
            if (array != pyforces)
                Py_XDECREF(array);
            PyGILState_Release(gstate);
        }
        void computeBatch(const std::vector<State>& states, double& energy, void* forces, bool forcesAreDouble) const override {
            if (batchComputation == nullptr)
                throw OpenMMException("PythonForce: Batch computation not configured");
            
            PyGILState_STATE gstate;
            gstate = PyGILState_Ensure();
            
            // Convert states vector to Python list
            PyObject* statesList = PyList_New(states.size());
            swig_type_info* info = SWIGTYPE_p_OpenMM__State;
            for (size_t i = 0; i < states.size(); i++) {
                PyObject* wrappedState = SWIG_NewPointerObj((void*) &states[i], info, 0);
                PyList_SetItem(statesList, i, wrappedState);
            }
            
            // Call batch computation function
            PyObject* result = PyObject_CallFunctionObjArgs(batchComputation, statesList, NULL);
            if (result == NULL) {
                // Handle exception
#if PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 12
                PyObject *type;
                PyObject *exception;
                PyObject *traceback;
                PyErr_Fetch(&type, &exception, &traceback);
#else
                PyObject *exception = PyErr_GetRaisedException();
#endif
                PyObject *message = PyObject_Str(exception);
                std::string *ptr;
                SWIG_AsPtr_std_string(message, &ptr);
                Py_XDECREF(message);
                Py_XDECREF(statesList);
                PyGILState_Release(gstate);
                throw OpenMMException(*ptr);
            }
            
            // Extract return values: (energy, forces_array)
            if (!PyTuple_Check(result) || PyTuple_Size(result) != 2) {
                Py_XDECREF(statesList);
                Py_XDECREF(result);
                PyGILState_Release(gstate);
                throw OpenMMException("PythonForce batch: Expected two return values (energy, forces)");
            }
            
            PyObject* pyenergy = Py_StripOpenMMUnits(PyTuple_GetItem(result, 0));
            PyObject* pyforces = Py_StripOpenMMUnits(PyTuple_GetItem(result, 1));
            energy = PyFloat_AsDouble(pyenergy);
            
            // Copy forces - should be (numCopies, numParticles, 3)
            if (!PyArray_Check(pyforces) || PyArray_NDIM((PyArrayObject*) pyforces) != 3) {
                Py_XDECREF(statesList);
                Py_XDECREF(result);
                PyGILState_Release(gstate);
                throw OpenMMException("PythonForce batch: Forces must be (numCopies, numParticles, 3) array");
            }
            
            npy_intp* dims = PyArray_DIMS((PyArrayObject*) pyforces);
            int numParticles = states[0].getPositions().size();
            if (dims[0] != (npy_intp)states.size() || dims[1] != numParticles || dims[2] != 3) {
                Py_XDECREF(statesList);
                Py_XDECREF(result);
                PyGILState_Release(gstate);
                throw OpenMMException("PythonForce batch: Forces array has wrong shape");
            }
            
            PyObject* array;
            int targetType = (forcesAreDouble ? NPY_DOUBLE : NPY_FLOAT);
            if (PyArray_CHKFLAGS((PyArrayObject*) pyforces, NPY_ARRAY_CARRAY_RO) && PyArray_DESCR((PyArrayObject*) pyforces)->type_num == targetType)
                array = pyforces;
            else
                array = PyArray_FromAny(pyforces, PyArray_DescrFromType(targetType), 3, 3, NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_FORCECAST, NULL);
            
            int elementSize = (forcesAreDouble ? sizeof(double) : sizeof(float));
            void* data = PyArray_DATA((PyArrayObject*) array);
            memcpy(forces, data, 3 * elementSize * numParticles * states.size());
            
            Py_XDECREF(statesList);
            Py_XDECREF(result);
            Py_XDECREF(pyenergy);
            Py_XDECREF(pyforces);
            if (array != pyforces)
                Py_XDECREF(array);
            PyGILState_Release(gstate);
        }
    private:
        PyObject* computation;
        PyObject* batchComputation;
    };

    /**
     * Construct a new PythonForce.
     */
    PythonForce* _createPythonForce(PyObject* computation, const std::map<std::string, double>& globalParameters={}, PyObject* batchComputation=nullptr) {
        PythonForce* force = new PythonForce(new ComputationWrapper(computation, batchComputation), globalParameters);
        PyObject* pickle = PyImport_ImportModule("pickle");
        PyObject* dumps = PyUnicode_FromString("dumps");
        PyObject* result = PyObject_CallMethodOneArg(pickle, dumps, computation);
        if (result == NULL) {
            // It couldn't be pickled.  It will still work, but can't be serialized.  Clear the error flag.
            PyErr_Clear();
        }
        else {
            char* buffer;
            Py_ssize_t len;
            if (PyBytes_AsStringAndSize(result, &buffer, &len) == 0)
                force->setPickledFunction(buffer, len);
        }
        return force;
    }

    /**
     * This is the serialization proxy used to serialize PythonForce objects.
     */
    class PythonForceProxy : public SerializationProxy {
    public:
        PythonForceProxy() : SerializationProxy("PythonForce") {
        }

        static std::string hexEncode(const std::vector<char>& input) {
            std::stringstream ss;
            ss << std::hex << std::setfill('0');
            for (unsigned char i : input)
                ss << std::setw(2) << static_cast<uint64_t>(i);
            return ss.str();
        }

        static std::vector<char> hexDecode(const std::string& input) {
            std::vector<char> res;
            res.reserve(input.size() / 2);
            for (size_t i = 0; i < input.length(); i += 2) {
                std::istringstream iss(input.substr(i, 2));
                uint64_t temp;
                iss >> std::hex >> temp;
                res.push_back(static_cast<unsigned char>(temp));
            }
            return res;
        }

        void serialize(const void* object, SerializationNode& node) const {
            node.setIntProperty("version", 1);
            const PythonForce& force = *reinterpret_cast<const PythonForce*>(object);
            if (force.getPickledFunction().size() == 0)
                throw OpenMMException("PythonForceProxy: Could not serialize PythonForce because its function could not be pickled.");
            node.setStringProperty("function", hexEncode(force.getPickledFunction()));
            node.setIntProperty("forceGroup", force.getForceGroup());
            node.setBoolProperty("usesPeriodic", force.usesPeriodicBoundaryConditions());
            SerializationNode& globalParams = node.createChildNode("GlobalParameters");
            for (auto param : force.getGlobalParameters())
                globalParams.createChildNode("Parameter").setStringProperty("name", param.first).setDoubleProperty("default", param.second);
        }

        void* deserialize(const SerializationNode& node) const {
            int version = node.getIntProperty("version");
            if (version != 1)
                throw OpenMMException("Unsupported version number");
            std::vector<char> pickledFunction = hexDecode(node.getStringProperty("function"));
            PyObject* pickle = PyImport_ImportModule("pickle");
            PyObject* loads = PyUnicode_FromString("loads");
            PyObject *pythonBytes = PyBytes_FromStringAndSize(pickledFunction.data(), pickledFunction.size());
            PyObject *function = PyObject_CallMethodOneArg(pickle, loads, pythonBytes);
            Py_XDECREF(pythonBytes);
            const SerializationNode& paramsNode = node.getChildNode("GlobalParameters");
            std::map<std::string, double> params;
            for (auto& parameter : paramsNode.getChildren())
                params[parameter.getStringProperty("name")] = parameter.getDoubleProperty("default");
            PythonForce* force = _createPythonForce(function, params);
            if (node.hasProperty("forceGroup"))
                force->setForceGroup(node.getIntProperty("forceGroup", 0));
            if (node.hasProperty("usesPeriodic"))
                force->setUsesPeriodicBoundaryConditions(node.getBoolProperty("usesPeriodic"));
            return force;
        }
    };

    /**
     * Register the serialization proxy.  This function is invoked automatically when the openmm module is imported.
     */
    void registerPythonForceProxy() {
        SerializationProxy::registerProxy(typeid(PythonForce), new PythonForceProxy());
    }
}

%}

%extend OpenMM::PythonForce {
    %feature("docstring") PythonForce "Create a PythonForce.

Parameters
----------
computation : function
    A function that performs the computation.  It should take a State as its argument
    and return two values: the potential energy (a scalar) and the forces (a NumPy array).
globalParameters : dict
    Any global parameters the function depends on.  Keys are the parameter names, and the
    corresponding values are their default values.
batchComputation : function, optional
    Optional batched computation function for RPMD. Takes a list of States and returns
    (total_energy, forces_array) where forces_array has shape (numCopies, numParticles, 3).
"
    PythonForce(PyObject* computation, const std::map<std::string, double>& globalParameters={}, PyObject* batchComputation=nullptr) {
        return _createPythonForce(computation, globalParameters, batchComputation);
    }
}
