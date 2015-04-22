%include pythoncode.i
%include exceptions.i
%include extend.i
%include header.i
%include pythonprepend_all.i
%include pythonprepend.i
%include pythonappend.i
%include typemaps.i
/* SWIG 3.x resolved a bug in which all wrapped C++ functions took *args as its
 * default argument list. OpenMM then exploited this bug by doing stuff like
 * passing args to stripUnits (and all added code assumed that the arguments
 * were in an "args" list). So in order to restore this arguably buggy behavior
 * from SWIG 2, enable the "compactdefaultargs" feature globally.
 *
 * See https://github.com/swig/swig/issues/387
 */
%feature("compactdefaultargs");
