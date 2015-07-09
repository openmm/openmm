%fragment("Py_StripOpenMMUnits", "header") {

static PyObject *__s_Quantity = NULL;
static PyObject *__s_md_unit_system_tuple = NULL;
static PyObject *__s_bar_tuple = NULL;

PyObject* Py_StripOpenMMUnits(PyObject *input) {
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
            Py_DECREF(item);
            return SWIG_ERROR;
        }
        double d = PyFloat_AsDouble(item1);
        Py_DECREF(item);
        Py_DECREF(item1);

        if (PyErr_Occurred() != NULL) {
            return SWIG_ERROR;
        }
        out.push_back(d);
    }
    Py_DECREF(iterator);
    Py_DECREF(stripped);
    return SWIG_OK;
 }
}


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


%typemap(in, fragment="Py_SequenceToVec3") Vec3 (int res=0){
    // typemap -- %typemap(in) Vec3
    $1 = Py_SequenceToVec3($input, res);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
}


%typemap(in, fragment="Py_SequenceToVec3") const Vec3& (OpenMM::Vec3 myVec, int res=0) {
    // typemap -- %typemap(in) Vec3
    myVec = Py_SequenceToVec3($input, res);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
    $1 = &myVec;
}


/* Convert python list of tuples to C++ std::vector of Vec3 objects */
%typemap(in, fragment="Py_SequenceToVec3") const std::vector<Vec3>& (std::vector<OpenMM::Vec3> vVec, PyObject* s=NULL, PyObject* o=NULL) {
    int i, pLength, ret;
    s = Py_StripOpenMMUnits($input);
    pLength = (int)PySequence_Length(s);
    for (i = 0; i < pLength; i++) {
        o = PySequence_GetItem(s, i);
        OpenMM::Vec3 v = Py_SequenceToVec3(o, ret);
        if (!SWIG_IsOK(ret)) {
          Py_DECREF(s);
          Py_DECREF(o);
          PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
          SWIG_fail;
        }
        vVec.push_back(v);
    }
    $1 = &vVec;
    Py_DECREF(s);
}



%typemap(in, fragment="Py_SequenceToVecDouble") const std::vector<double> & (std::vector<double> v, int res=0) {
    res = Py_SequenceToVecDouble($input, v);
    if (!SWIG_IsOK(res)) {
        PyErr_SetString(PyExc_ValueError, "in method $symname, argument $argnum could not be converted to type $type");
        SWIG_fail;
    }
    $1 = &v;
}

/* The following two typemaps cause a non-const vector<Vec3>& to become a return value. */
%typemap(in, numinputs=0) std::vector<Vec3>& (std::vector<Vec3> temp) {
    $1 = &temp;
}


%typemap(argout) std::vector<Vec3>& {
    int i, n;
    PyObject *pyList;

    n=(*$1).size();
    pyList=PyList_New(n);
    PyObject* mm = PyImport_AddModule("simtk.openmm");
    PyObject* vec3 = PyObject_GetAttrString(mm, "Vec3");
    for (i=0; i<n; i++) {
        OpenMM::Vec3& v = (*$1).at(i);
        PyObject* args = Py_BuildValue("(d,d,d)", v[0], v[1], v[2]);
        PyObject* pyVec = PyObject_CallObject(vec3, args);
        Py_DECREF(args);
        PyList_SET_ITEM(pyList, i, pyVec);
    }
    $result = pyList;
}


/* const vector<Vec3> should NOT become an output. */
%typemap(argout) const std::vector<Vec3>& {
}

/* Convert python tuple to C++ Vec3 object*/
%typemap(typecheck) Vec3 {
    // typemap -- %typemap(typecheck) Vec3
    $1 = (PySequence_Length($input) >= 3 ? 1 : 0);
}


%typemap(typecheck) const Vec3& {
    // typemap -- %typemap(typecheck) Vec3
    $1 = (PySequence_Length($input) >= 3 ? 1 : 0);
}


%typemap(out) Vec3 {
    PyObject* mm = PyImport_AddModule("simtk.openmm");
    PyObject* vec3 = PyObject_GetAttrString(mm, "Vec3");
    PyObject* args = Py_BuildValue("(d,d,d)", ($1)[0], ($1)[1], ($1)[2]);
    $result = PyObject_CallObject(vec3, args);
    Py_DECREF(args);
}


%typemap(out) const Vec3& {
    PyObject* mm = PyImport_AddModule("simtk.openmm");
    PyObject* vec3 = PyObject_GetAttrString(mm, "Vec3");
    PyObject* args = Py_BuildValue("(d,d,d)", (*$1)[0], (*$1)[1], (*$1)[2]);
    $result = PyObject_CallObject(vec3, args);
    Py_DECREF(args);
}

/* Convert C++ (Vec3&, Vec3&, Vec3&) object to python tuple or tuples */
%typemap(argout) (Vec3& a, Vec3& b, Vec3& c) {
    // %typemap(argout) (Vec3& a, Vec3& b, Vec3& c)
    PyObject* mm = PyImport_AddModule("simtk.openmm");
    PyObject* vec3 = PyObject_GetAttrString(mm, "Vec3");
    PyObject* args1 = Py_BuildValue("(d,d,d)", (*$1)[0], (*$1)[1], (*$1)[2]);
    PyObject* args2 = Py_BuildValue("(d,d,d)", (*$2)[0], (*$2)[1], (*$2)[2]);
    PyObject* args3 = Py_BuildValue("(d,d,d)", (*$3)[0], (*$3)[1], (*$3)[2]);
    PyObject* pyVec1 = PyObject_CallObject(vec3, args1);
    PyObject* pyVec2 = PyObject_CallObject(vec3, args2);
    PyObject* pyVec3 = PyObject_CallObject(vec3, args3);
    Py_DECREF(args1);
    Py_DECREF(args2);
    Py_DECREF(args3);
    Py_DECREF(mm);
    Py_DECREF(vec3);
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
