%module openmm
%include "factory.i"

%include "std_string.i"
%include "std_iostream.i"
%include "std_map.i"
%include "std_pair.i"
%include "std_set.i"
%include "std_vector.i"

namespace std {
  %template(pairii) pair<int,int>;
  %template(vectord) vector<double>;
  %template(vectorddd) vector< vector< vector<double> > >;
  %template(vectori) vector<int>;
  %template(vectorii) vector < vector<int> >;
  %template(vectorpairii) vector< pair<int,int> >;
  %template(vectorstring) vector<string>;
  %template(mapstringstring) map<string,string>;
  %template(mapstringdouble) map<string,double>;
  %template(mapii) map<int,int>;
  %template(seti) set<int>;
};

%include "typemaps.i"
%include "windows.i"

%{
#define SWIG_FILE_WITH_INIT

#include <sstream>

#include <exception>
#include <fstream>
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
#include "OpenMMDrude.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/SerializationProxy.h"
#include "openmm/serialization/XmlSerializer.h"

using namespace OpenMM;

%}

%feature("autodoc", "0");
%nodefaultctor;

%include features.i

%include OpenMM_docstring.i

%include OpenMMSwigHeaders.i

%pythoncode %{
  # when we import * from the python module, we only want to import the
  # actual classes, and not the swigregistration methods, which have already
  # been called, and are now unneeded by the user code, and only pollute the
  # namespace
  __all__ = [k for k in locals().keys() if not (k.endswith('_swigregister') or k.startswith('_'))]
%}
