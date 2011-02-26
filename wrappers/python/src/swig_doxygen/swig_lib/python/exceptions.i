
%exception {
    try {
        $action
    } catch (OpenMM::OpenMMException &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}


