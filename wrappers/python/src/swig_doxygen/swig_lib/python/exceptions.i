
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}

%exception *::step {
    try {
        Py_BEGIN_ALLOW_THREADS
        $action
        Py_END_ALLOW_THREADS
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}

%exception OpenMM::LocalEnergyMinimizer::minimize {
    try {
        Py_BEGIN_ALLOW_THREADS
        $action
        Py_END_ALLOW_THREADS
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}
