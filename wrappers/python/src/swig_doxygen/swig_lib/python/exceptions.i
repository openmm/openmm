
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}

%exception *::step {
    PyThreadState* _savePythonThreadState = PyEval_SaveThread();
    try {
        $action
    } catch (std::exception &e) {
        PyEval_RestoreThread(_savePythonThreadState);
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
    PyEval_RestoreThread(_savePythonThreadState);
}

%exception OpenMM::LocalEnergyMinimizer::minimize {
    PyThreadState* _savePythonThreadState = PyEval_SaveThread();
    try {
        $action
    } catch (std::exception &e) {
        PyEval_RestoreThread(_savePythonThreadState);
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
    PyEval_RestoreThread(_savePythonThreadState);
}
