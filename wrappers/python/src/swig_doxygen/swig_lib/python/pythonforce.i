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
        ComputationWrapper(PyObject* computation) : computation(computation) {
        }
        void compute(const State& state, double& energy, std::vector<Vec3>& forces) const {
            PyGILState_STATE gstate;
            gstate = PyGILState_Ensure();

            // Invoke the function.

            swig_type_info* info = SWIGTYPE_p_OpenMM__State;
            PyObject* wrappedState = SWIG_NewPointerObj((void*) &state, info, 0);
            PyObject* result = PyObject_CallFunctionObjArgs(computation, wrappedState, NULL);
            if (result == NULL) {
                // The function raised an exception.  Convert it to an OpenMMException.

                PyObject *exception = PyErr_GetRaisedException();
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

    PythonForce* _createPythonForce(PyObject* computation, const std::map<std::string, double>& globalParameters={}) {
        PythonForce* force = new PythonForce(new ComputationWrapper(computation), globalParameters);
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

    void registerPythonForceProxy() {
        SerializationProxy::registerProxy(typeid(PythonForce), new PythonForceProxy());
    }
}

%}

%extend OpenMM::PythonForce {
    PythonForce(PyObject* computation, const std::map<std::string, double>& globalParameters={}) {
        return _createPythonForce(computation, globalParameters);
    }
}
