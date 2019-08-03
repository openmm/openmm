%fragment("Vec3_to_PyVec3", "header") {
/**
 * Convert an OpenMM::Vec3 into a Python simtk.openmm.Vec3 object
 *
 * Returns a new reference.
 */
PyObject* Vec3_to_PyVec3(const OpenMM::Vec3& v) {
    static PyObject *__s_mm = NULL;
    static PyObject *__s_Vec3 = NULL;
    if (__s_mm == NULL) {
        __s_mm = PyImport_AddModule("simtk.openmm");
        __s_Vec3 = PyObject_GetAttrString(__s_mm, "Vec3");
    }
    PyObject* tuple = Py_BuildValue("(d,d,d)", v[0], v[1], v[2]);
    PyObject* PyVec3 = PyObject_CallObject(__s_Vec3, tuple);
    Py_DECREF(tuple);
    return PyVec3;
}
}

%fragment("Py_StripOpenMMUnits", "header") {

/**
 * Strip any OpenMM units of an input PyObject.
 *
 * This is equivalent to the following Python code
 *
 * >>> from simtk import unit
 * >>> if isinstance(input, unit.Quantity)
 * ...     if input.is_compatible(unit.bar)
 * ...         return input.value_in_unit(unit.bar)
 * ...     return input.value_in_unit_system(unit.md_input_system)
 * ... return input
 *
 * Returns a new reference.
 */
PyObject* Py_StripOpenMMUnits(PyObject *input) {
    static PyObject *__s_Quantity = NULL;
    static PyObject *__s_md_unit_system_tuple = NULL;
    static PyObject *__s_bar_tuple = NULL;

    if (__s_Quantity == NULL) {
        PyObject* module = NULL;
        module = PyImport_ImportModule("simtk.unit");
        if (!module) {
            PyErr_SetString(PyExc_ImportError, "simtk.unit"); Py_CLEAR(module); return NULL;
        }

        __s_Quantity = PyObject_GetAttrString(module, "Quantity");
        if (!__s_Quantity) {
            PyErr_SetString(PyExc_AttributeError, "'module' object has no attribute 'Quantity'");
            Py_CLEAR(module);
            Py_CLEAR(__s_Quantity);
            return NULL;
        }

        PyObject* bar = NULL;
        bar = PyObject_GetAttrString(module, "bar");
        if (!bar) {
            PyErr_SetString(PyExc_AttributeError, "'module' object has no attribute 'bar'");
            Py_CLEAR(module);
            Py_CLEAR(__s_Quantity);
            Py_CLEAR(bar);
            return NULL;
        }

        PyObject* md_unit_system = NULL;
        md_unit_system = PyObject_GetAttrString(module, "md_unit_system");
        if (!md_unit_system) {
            PyErr_SetString(PyExc_AttributeError, "'module' object has no attribute 'md_unit_system'");
            Py_CLEAR(module);
            Py_CLEAR(__s_Quantity);
            Py_CLEAR(bar);
           Py_CLEAR(md_unit_system);
        }
        __s_md_unit_system_tuple = PyTuple_Pack(1, md_unit_system);
        __s_bar_tuple = PyTuple_Pack(1, bar);
        Py_DECREF(md_unit_system);
        Py_DECREF(bar);
        Py_DECREF(module);
    }
    PyObject *val;

    if (PyObject_IsInstance(input, __s_Quantity)) {
        PyObject* input_unit = NULL, *is_compatible = NULL, *compatible_with_bar = NULL;
        input_unit = PyObject_GetAttrString(input, "unit");
        is_compatible = PyObject_GetAttrString(input_unit, "is_compatible");
        compatible_with_bar = PyObject_Call(is_compatible, __s_bar_tuple, NULL);
        if (PyObject_IsTrue(compatible_with_bar)) {
            // input.in_units_of(unit.bar)
            PyObject* value_in_unit = PyObject_GetAttrString(input, "value_in_unit");
            val = PyObject_Call(value_in_unit, __s_bar_tuple, NULL);
            Py_DECREF(value_in_unit);
        } else {
            // input.value_in_unit_system(md_unit_system_tuple)
            PyObject* value_in_unit_system = PyObject_GetAttrString(input, "value_in_unit_system");
            val = PyObject_Call(value_in_unit_system, __s_md_unit_system_tuple, NULL);
            Py_DECREF(value_in_unit_system);
        }
        Py_CLEAR(input_unit);
        Py_CLEAR(is_compatible);
        Py_CLEAR(compatible_with_bar);
        if (PyErr_Occurred() != NULL) {
            return NULL;
        }
    } else {
        val = input;
        Py_INCREF(val);
    }
    return val;
}
}


%fragment("Py_SequenceToVec3", "header", fragment="Py_StripOpenMMUnits") {
OpenMM::Vec3 Py_SequenceToVec3(PyObject* obj, int& status) {
    PyObject* s, *o, *o1;
    double x[3];
    int i, length;
    s = Py_StripOpenMMUnits(obj);
    if (s == NULL) {
        status = SWIG_ERROR;
        return OpenMM::Vec3(0, 0, 0);
    }

    length = (int) PySequence_Length(s);
    if (length != 3) {
        Py_DECREF(s);
        PyErr_SetString(PyExc_TypeError, "Item must have length 3");
        status = SWIG_ERROR;
        return OpenMM::Vec3(0, 0, 0);
    }

    for (i = 0; i < 3; i++ ) {
        o = PySequence_GetItem(s, i);
        o1 = Py_StripOpenMMUnits(o);
        if (o1 != NULL) {
            x[i] = PyFloat_AsDouble(o1);
        }
        if (o1 == NULL || PyErr_Occurred() != NULL) {
            Py_DECREF(s);
            Py_DECREF(o);
            Py_XDECREF(o1);
            status = SWIG_ERROR;
            return OpenMM::Vec3(0, 0, 0);
        }
        Py_DECREF(o);
        Py_DECREF(o1);
    }

    status = SWIG_OK;
    Py_DECREF(s);
    return OpenMM::Vec3(x[0], x[1], x[2]);
  }
}

%fragment("Py_SequenceToVecDouble", "header", fragment="Py_StripOpenMMUnits") {
int Py_SequenceToVecDouble(PyObject* obj, std::vector<double>& out) {
    PyObject* stripped = Py_StripOpenMMUnits(obj);
    PyObject* item = NULL;
    PyObject* item1 = NULL;

    if (isNumpyAvailable()) {
        if (PyArray_Check(stripped) && PyArray_ISCARRAY_RO(stripped) && PyArray_NDIM(stripped) == 1) {
            int type = PyArray_TYPE(stripped);
            int length = PyArray_SIZE(stripped);
            void* data = PyArray_DATA((PyArrayObject*) stripped);
            if (type == NPY_DOUBLE) {
                out.resize(length);
                memcpy(&out[0], data, sizeof(double)*length);
                Py_DECREF(stripped);
                return SWIG_OK;
            }
            if (type == NPY_FLOAT) {
                out.resize(length);
                float* floatData = (float*) data;
                for (int i = 0; i < length; i++)
                    out[i] = floatData[i];
                Py_DECREF(stripped);
                return SWIG_OK;
            }
            if (type == NPY_INT32) {
                out.resize(length);
                int* intData = (int*) data;
                for (int i = 0; i < length; i++)
                    out[i] = intData[i];
                Py_DECREF(stripped);
                return SWIG_OK;
            }
            if (type == NPY_INT64) {
                out.resize(length);
                long long* longData = (long long*) data;
                for (int i = 0; i < length; i++)
                    out[i] = longData[i];
                Py_DECREF(stripped);
                return SWIG_OK;
            }
        }
    }

    PyObject* iterator = PyObject_GetIter(stripped);
    if (iterator == NULL) {
        Py_DECREF(stripped);
        return SWIG_ERROR;
    }

    while ((item = PyIter_Next(iterator))) {
        item1 = Py_StripOpenMMUnits(item);
        if (item1 == NULL) {
            Py_DECREF(stripped);
            Py_DECREF(iterator);
            Py_DECREF(item);
            return SWIG_ERROR;
        }
        double d = PyFloat_AsDouble(item1);
        Py_DECREF(item);
        Py_DECREF(item1);

        if (PyErr_Occurred() != NULL) {
            Py_DECREF(stripped);
            Py_DECREF(iterator);
            return SWIG_ERROR;
        }
        out.push_back(d);
    }
    Py_DECREF(iterator);
    Py_DECREF(stripped);
    return SWIG_OK;
}
}

%fragment("Py_SequenceToVecVec3", "header", fragment="Py_SequenceToVec3") {
int Py_SequenceToVecVec3(PyObject* obj, std::vector<Vec3>& out) {
    PyObject* stripped = Py_StripOpenMMUnits(obj);      // new reference
    if (isNumpyAvailable()) {
        if (PyArray_Check(stripped) && PyArray_ISCARRAY_RO(stripped) && PyArray_NDIM(stripped) == 2 && PyArray_DIM(stripped, 1) == 3) {
            int type = PyArray_TYPE(stripped);
            int length = PyArray_DIM(stripped, 0);
            void* data = PyArray_DATA((PyArrayObject*) stripped);
            if (type == NPY_DOUBLE) {
                out.resize(length);
                memcpy(&out[0][0], data, 3*sizeof(double)*length);
                Py_DECREF(stripped);
                return SWIG_OK;
            }
            if (type == NPY_FLOAT) {
                out.resize(length);
                float* floatData = (float*) data;
                for (int i = 0; i < length; i++)
                    out[i] = Vec3(floatData[3*i], floatData[3*i+1], floatData[3*i+2]);
                Py_DECREF(stripped);
                return SWIG_OK;
            }
            if (type == NPY_INT32) {
                out.resize(length);
                int* intData = (int*) data;
                for (int i = 0; i < length; i++)
                    out[i] = Vec3(intData[3*i], intData[3*i+1], intData[3*i+2]);
                Py_DECREF(stripped);
                return SWIG_OK;
            }
            if (type == NPY_INT64) {
                out.resize(length);
                long long* longData = (long long*) data;
                for (int i = 0; i < length; i++)
                    out[i] = Vec3(longData[3*i], longData[3*i+1], longData[3*i+2]);
                Py_DECREF(stripped);
                return SWIG_OK;
            }
        }
    }
    int ret = 0;
    PyObject* item = NULL;
    PyObject* item1 = NULL;
    PyObject* iterator = NULL;
    iterator = PyObject_GetIter(stripped);    // new reference

    if (iterator == NULL) {
        Py_DECREF(stripped);
        return SWIG_ERROR;
    }

    while ((item = PyIter_Next(iterator))) {  // new reference
        item1 = Py_StripOpenMMUnits(item);    // new reference
        if (item1 == NULL) {
            Py_DECREF(stripped);
            Py_DECREF(iterator);
            Py_DECREF(item);
            return SWIG_ERROR;
        }

        OpenMM::Vec3 v = Py_SequenceToVec3(item1, ret);
        Py_DECREF(item);
        Py_DECREF(item1);

        if (!SWIG_IsOK(ret)) {
            Py_DECREF(stripped);
            Py_DECREF(iterator);
            return SWIG_ERROR;
        }

        out.push_back(v);
    }

    Py_DECREF(iterator);
    Py_DECREF(stripped);
    return SWIG_OK;
}
}

%fragment("Py_SequenceToVecVecDouble", "header", fragment="Py_SequenceToVecDouble") {
int Py_SequenceToVecVecDouble(PyObject* obj, std::vector<std::vector<double> >& out) {
    PyObject* stripped = NULL;
    PyObject* item = NULL;
    PyObject* item1 = NULL;
    PyObject* iterator = NULL;
    stripped = Py_StripOpenMMUnits(obj);
    iterator = PyObject_GetIter(stripped);

    if (iterator == NULL) {
        Py_DECREF(stripped);
        return SWIG_ERROR;
    }

    while ((item = PyIter_Next(iterator))) {
        item1 = Py_StripOpenMMUnits(item);
        if (item1 == NULL) {
            Py_DECREF(stripped);
            Py_DECREF(iterator);
            Py_DECREF(item);
            return SWIG_ERROR;
        }
        std::vector<double> v;
        int r2 = Py_SequenceToVecDouble(item1, v);
        Py_DECREF(item);
        Py_DECREF(item1);

        if (!SWIG_IsOK(r2) || PyErr_Occurred() != NULL) {
            Py_DECREF(stripped);
            Py_DECREF(iterator);
            return SWIG_ERROR;
        }
        out.push_back(v);
    }
    Py_DECREF(iterator);
    Py_DECREF(stripped);
    return SWIG_OK;
}
}

%fragment("Py_SequenceToVecVecVecDouble", "header", fragment="Py_SequenceToVecVecDouble") {
int Py_SequenceToVecVecVecDouble(PyObject* obj, std::vector<std::vector<std::vector<double> > >& out) {
    PyObject* stripped = NULL;
    PyObject* item = NULL;
    PyObject* item1 = NULL;
    PyObject* iterator = NULL;
    stripped = Py_StripOpenMMUnits(obj);
    iterator = PyObject_GetIter(stripped);

    if (iterator == NULL) {
        Py_DECREF(stripped);
        return SWIG_ERROR;
    }

    while ((item = PyIter_Next(iterator))) {
        item1 = Py_StripOpenMMUnits(item);
        if (item1 == NULL) {
            Py_DECREF(stripped);
            Py_DECREF(iterator);
            Py_DECREF(item);
            return SWIG_ERROR;
        }
        std::vector<std::vector<double> >v;
        int r2 = Py_SequenceToVecVecDouble(item1, v);
        Py_DECREF(item);
        Py_DECREF(item1);

        if (!SWIG_IsOK(r2) || PyErr_Occurred() != NULL) {
            Py_DECREF(stripped);
            Py_DECREF(iterator);
            return SWIG_ERROR;
        }
        out.push_back(v);
    }
    Py_DECREF(iterator);
    Py_DECREF(stripped);
    return SWIG_OK;
}
}


// ------ typemap for double ----
%typemap(typecheck, precedence=SWIG_TYPECHECK_DOUBLE, fragment="Py_StripOpenMMUnits") double {
    double argp = 0;
    PyObject* s = NULL;
    s = Py_StripOpenMMUnits($input);
    $1 = (s != NULL) ? SWIG_IsOK(SWIG_AsVal_double(s, &argp)) : 0;
    Py_DECREF(s);
}
%typemap(in, noblock=1, fragment="Py_StripOpenMMUnits") double (double argp = 0, int res = 0,
    PyObject* stripped = NULL) {

    stripped = Py_StripOpenMMUnits($input);
    if (stripped == NULL) { SWIG_fail; }
    res = SWIG_AsVal_double(stripped, &argp);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
    $1 = ($ltype)(argp);
    Py_CLEAR(stripped);
}


// ------ typemap for Vec3
%typemap(in, fragment="Py_SequenceToVec3") Vec3 (int res=0){
    // typemap -- %typemap(in) Vec3
    $1 = Py_SequenceToVec3($input, res);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
}
%typemap(typecheck, fragment="Py_SequenceToVec3") Vec3 {
    int res = 0;
    Py_SequenceToVec3($input, res);
    $1 = SWIG_IsOK(res);
}


// typemap for const Vec3&
%typemap(in, fragment="Py_SequenceToVec3") const Vec3& (OpenMM::Vec3 myVec, int res=0) {
    // typemap -- %typemap(in) Vec3
    myVec = Py_SequenceToVec3($input, res);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
    $1 = &myVec;
}
%typemap(typecheck, fragment="Py_SequenceToVec3") const Vec3& {
    int res = 0;
    Py_SequenceToVec3($input, res);
    $1 = SWIG_IsOK(res);
}


// typemap for const vector<vector<double> >
%typemap(in, fragment="Py_SequenceToVecVecDouble") const std::vector<std::vector<double> >(std::vector<std::vector<double> > v, int res=0) {
    res = Py_SequenceToVecVecDouble($input, v);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
    $1 = &v;
}
%typemap(typecheck, fragment="Py_SequenceToVecVecDouble") const std::vector<std::vector<double> > {
    std::vector<std::vector<double> > v;
    int res = Py_SequenceToVecVecDouble($input, v);
    $1 = SWIG_IsOK(res);
}


// typemap for const vector<vector<vector<double> > >
%typemap(in, fragment="Py_SequenceToVecVecVecDouble") const std::vector<std::vector<std::vector<double> > >(std::vector<std::vector<std::vector<double> > > v, int res=0) {
    res = Py_SequenceToVecVecVecDouble($input, v);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
    $1 = &v;
}
%typemap(typecheck, fragment="Py_SequenceToVecVecVecDouble") const std::vector<std::vector<std::vector<double> > >{
    std::vector<std::vector<std::vector<double> > > v;
    int res = Py_SequenceToVecVecVecDouble($input, v);
    $1 = SWIG_IsOK(res);
}


// typemap for vector<Vec3>
%typemap(in, fragment="Py_SequenceToVecVec3") const std::vector<Vec3>& (std::vector<Vec3> v, int res=0) {
    res = Py_SequenceToVecVec3($input, v);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
    $1 = &v;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_DOUBLE_ARRAY, fragment="Py_SequenceToVecVec3") const std::vector<Vec3>& {
    std::vector<Vec3> v;
    int res=0;
    res = Py_SequenceToVecVec3($input, v);
    $1 = SWIG_IsOK(res);
}
%typemap(out) const std::vector<Vec3>& {
    $result = copyVVec3ToList(*$1);
}


// typemap for const vector<double>
%typemap(in, fragment="Py_SequenceToVecDouble") const std::vector<double> & (std::vector<double> v, int res=0) {
    res = Py_SequenceToVecDouble($input, v);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
    $1 = &v;
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_DOUBLE_ARRAY) const std::vector<double> & {
    std::vector<double> v;
    int res = 0;
    res = Py_SequenceToVecDouble($input, v);
    $1 = SWIG_IsOK(res);
}



/* The following two typemaps cause a non-const vector<Vec3>& to become a return value. */
%typemap(in, numinputs=0) std::vector<Vec3>& (std::vector<Vec3> temp) {
    $1 = &temp;
}


%typemap(argout, fragment="Vec3_to_PyVec3") std::vector<Vec3>& {
    int n = (*$1).size();
    PyObject * pyList = PyList_New(n);
    for (int i=0; i<n; i++) {
        OpenMM::Vec3& v = (*$1).at(i);
        PyObject* pyVec = Vec3_to_PyVec3(v);
        PyList_SET_ITEM(pyList, i, pyVec);
    }
    $result = pyList;
}


/* const vector<Vec3> should NOT become an output. */
%typemap(argout) const std::vector<Vec3>& {
}


%typemap(out, fragment="Vec3_to_PyVec3") Vec3 {
    $result = Vec3_to_PyVec3(*$1);
}


%typemap(out, fragment="Vec3_to_PyVec3") const Vec3& {
    $result = Vec3_to_PyVec3(*$1);
}

/* Convert C++ (Vec3&, Vec3&, Vec3&) object to python tuple or tuples */
%typemap(argout, fragment="Vec3_to_PyVec3") (Vec3& a, Vec3& b, Vec3& c) {
    PyObject* pyVec1 = Vec3_to_PyVec3(*$1);
    PyObject* pyVec2 = Vec3_to_PyVec3(*$2);
    PyObject* pyVec3 = Vec3_to_PyVec3(*$3);
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


%typemap(out) std::string OpenMM::Context::createCheckpoint{
    // createCheckpoint returns a bytes object
    $result = PyBytes_FromStringAndSize($1.c_str(), $1.length());
}


%typemap(in) std::string {
    // if we have a C++ method that takes in a std::string, we're most happy
    // to accept a python bytes object. But if the user passes in a unicode
    // object we'll try to recover by encoding it to UTF-8 bytes
    PyObject* temp = NULL;
    char* c_str = NULL;
    Py_ssize_t len = 0;

    if (PyUnicode_Check($input)) {
        temp = PyUnicode_AsUTF8String($input);
        if (temp == NULL) {
            SWIG_exception_fail(SWIG_TypeError, "'utf-8' codec can't decode byte");
        }
        PyBytes_AsStringAndSize(temp, &c_str, &len);
        Py_XDECREF(temp);
    } else if (PyBytes_Check($input)) {
        PyBytes_AsStringAndSize($input, &c_str, &len);
    } else {
         SWIG_exception_fail(SWIG_TypeError, "argument must be str or bytes");
    }

    if (c_str == NULL) {
        SWIG_exception_fail(SWIG_TypeError, "argument must be str or bytes");
    }

    $1 = std::string(c_str, len);
}
