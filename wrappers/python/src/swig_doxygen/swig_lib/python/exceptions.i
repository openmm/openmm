
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyObject* mm = PyImport_AddModule("openmm");
        PyObject* openmm_exception = PyObject_GetAttrString(mm, "OpenMMException");
        PyErr_SetString(openmm_exception, const_cast<char*>(e.what()));
        return NULL;
    }
}

%exception *::getState {
    PyThreadState* _savePythonThreadState = PyEval_SaveThread();
    try {
        $action
    } catch (std::exception &e) {
        PyEval_RestoreThread(_savePythonThreadState);
        PyObject* mm = PyImport_AddModule("openmm");
        PyObject* openmm_exception = PyObject_GetAttrString(mm, "OpenMMException");
        PyErr_SetString(openmm_exception, const_cast<char*>(e.what()));
        return NULL;
    }
    PyEval_RestoreThread(_savePythonThreadState);
}

%exception *::step {
    PyThreadState* _savePythonThreadState = PyEval_SaveThread();
    try {
        $action
    } catch (std::exception &e) {
        PyEval_RestoreThread(_savePythonThreadState);
        PyObject* mm = PyImport_AddModule("openmm");
        PyObject* openmm_exception = PyObject_GetAttrString(mm, "OpenMMException");
        PyErr_SetString(openmm_exception, const_cast<char*>(e.what()));
        return NULL;
    }
    PyEval_RestoreThread(_savePythonThreadState);
}

%exception OpenMM::LocalEnergyMinimizer::minimize {
    bool releaseGIL = (nobjs < 4 || swig_obj[3] == Py_None);
    PyThreadState* _savePythonThreadState = (releaseGIL ? PyEval_SaveThread() : nullptr);
    try {
        $action
    }
    catch (std::exception &e) {
        if (releaseGIL)
            PyEval_RestoreThread(_savePythonThreadState);
        if (dynamic_cast<Swig::DirectorException*>(&e) != NULL) {
            SWIG_fail;
        }
        else {
            PyObject* mm = PyImport_AddModule("openmm");
            PyObject* openmm_exception = PyObject_GetAttrString(mm, "OpenMMException");
            PyErr_SetString(openmm_exception, const_cast<char*>(e.what()));
            return NULL;
        }
    }
    if (releaseGIL)
        PyEval_RestoreThread(_savePythonThreadState);
}

%exception OpenMM::Context::setVelocitiesToTemperature {
    PyThreadState* _savePythonThreadState = PyEval_SaveThread();
    try {
        $action
    } catch (std::exception &e) {
        PyEval_RestoreThread(_savePythonThreadState);
        PyObject* mm = PyImport_AddModule("openmm");
        PyObject* openmm_exception = PyObject_GetAttrString(mm, "OpenMMException");
        PyErr_SetString(openmm_exception, const_cast<char*>(e.what()));
        return NULL;
    }
    PyEval_RestoreThread(_savePythonThreadState);
}

%exception OpenMM::CustomCVForce::getCollectiveVariableValues {
    PyThreadState* _savePythonThreadState = PyEval_SaveThread();
    try {
        $action
    } catch (std::exception &e) {
        PyEval_RestoreThread(_savePythonThreadState);
        PyObject* mm = PyImport_AddModule("openmm");
        PyObject* openmm_exception = PyObject_GetAttrString(mm, "OpenMMException");
        PyErr_SetString(openmm_exception, const_cast<char*>(e.what()));
        return NULL;
    }
    PyEval_RestoreThread(_savePythonThreadState);
}