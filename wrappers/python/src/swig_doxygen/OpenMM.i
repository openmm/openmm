
%define DOCSTRING
"PyOpenMM is a Python application programming interface (API) to be
used for performing molecular dynamics (MD) simulations on various
computer architectures (including GPUs).  It is implemented in Python
and C/C++, and provides a Python interface to the OpenMM libraries
(see https://simtk.org/home/openmm for OpenMM details).  The primary
motivation for creating PyOpenMM is to make it possible to write
GPU-accelerated MD code in pure Python.

See https://simtk.org/home/pyopenmm for details"
%enddef

%module (docstring=DOCSTRING) openmm

%include "factory.i"
%include "std_string.i"
%include "std_iostream.i"
%include "typemaps.i"

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

%include "windows.i"

%{
#define SWIG_FILE_WITH_INIT

#include <sstream>

#include <exception>
#include <fstream>
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "openmm/RPMDIntegrator.h"
#include "OpenMMDrude.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/SerializationProxy.h"
#include "openmm/serialization/XmlSerializer.h"

using namespace OpenMM;

%}

%feature("autodoc", "1");
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

/*
%extend OpenMM::XmlSerializer {
    %template(XmlSerializer_serialize_AndersenThermostat) XmlSerializer::serialize<AndersenThermostat>;
    %template(XmlSerializer_serialize_RBTorsionForce) XmlSerializer::serialize<RBTorsionForce>;
    %template(XmlSerializer_serialize_CMAPTorsionForce) XmlSerializer::serialize<CMAPTorsionForce>;
    %template(XmlSerializer_serialize_CMMotionRemover) XmlSerializer::serialize<CMMotionRemover>;
    %template(XmlSerializer_serialize_CustomAngleForce) XmlSerializer::serialize<CustomAngleForce>;
    %template(XmlSerializer_serialize_CustomBondForce) XmlSerializer::serialize<CustomBondForce>;
    %template(XmlSerializer_serialize_CustomExternalForce) XmlSerializer::serialize<CustomExternalForce>;
    %template(XmlSerializer_serialize_CustomGBForce) XmlSerializer::serialize<CustomGBForce>;
    %template(XmlSerializer_serialize_CustomHbondForce) XmlSerializer::serialize<CustomHbondForce>;
    %template(XmlSerializer_serialize_CustomNonbondedForce) XmlSerializer::serialize<CustomNonbondedForce>;
    %template(XmlSerializer_serialize_CustomTorsionForce) XmlSerializer::serialize<CustomTorsionForce>;
    %template(XmlSerializer_serialize_GBSAOBCForce) XmlSerializer::serialize<GBSAOBCForce>;
    %template(XmlSerializer_serialize_GBVIForce) XmlSerializer::serialize<GBVIForce>;
    %template(XmlSerializer_serialize_HarmonicAngleForce) XmlSerializer::serialize<HarmonicAngleForce>;
    %template(XmlSerializer_serialize_HarmonicBondForce) XmlSerializer::serialize<HarmonicBondForce>;
    %template(XmlSerializer_serialize_MonteCarloBarostat) XmlSerializer::serialize<MonteCarloBarostat>;
    %template(XmlSerializer_serialize_MonteCarloAnisotropicBarostat) XmlSerializer::serialize<MonteCarloAnisotropicBarostat>;
    %template(XmlSerializer_serialize_NonbondedForce) XmlSerializer::serialize<NonbondedForce>;
    %template(XmlSerializer_serialize_RBTorsionForce) XmlSerializer::serialize<RBTorsionForce>;
    %template(XmlSerializer_serialize_System) XmlSerializer::serialize<System>;

    %template(XmlSerializer_deserialize_AndersenThermostat) XmlSerializer::deserialize<AndersenThermostat>;
    %template(XmlSerializer_deserialize_RBTorsionForce) XmlSerializer::deserialize<RBTorsionForce>;
    %template(XmlSerializer_deserialize_CMAPTorsionForce) XmlSerializer::deserialize<CMAPTorsionForce>;
    %template(XmlSerializer_deserialize_CMMotionRemover) XmlSerializer::deserialize<CMMotionRemover>;
    %template(XmlSerializer_deserialize_CustomAngleForce) XmlSerializer::deserialize<CustomAngleForce>;
    %template(XmlSerializer_deserialize_CustomBondForce) XmlSerializer::deserialize<CustomBondForce>;
    %template(XmlSerializer_deserialize_CustomExternalForce) XmlSerializer::deserialize<CustomExternalForce>;
    %template(XmlSerializer_deserialize_CustomGBForce) XmlSerializer::deserialize<CustomGBForce>;
    %template(XmlSerializer_deserialize_CustomHbondForce) XmlSerializer::deserialize<CustomHbondForce>;
    %template(XmlSerializer_deserialize_CustomNonbondedForce) XmlSerializer::deserialize<CustomNonbondedForce>;
    %template(XmlSerializer_deserialize_CustomTorsionForce) XmlSerializer::deserialize<CustomTorsionForce>;
    %template(XmlSerializer_deserialize_GBSAOBCForce) XmlSerializer::deserialize<GBSAOBCForce>;
    %template(XmlSerializer_deserialize_GBVIForce) XmlSerializer::deserialize<GBVIForce>;
    %template(XmlSerializer_deserialize_HarmonicAngleForce) XmlSerializer::deserialize<HarmonicAngleForce>;
    %template(XmlSerializer_deserialize_HarmonicBondForce) XmlSerializer::deserialize<HarmonicBondForce>;
    %template(XmlSerializer_deserialize_MonteCarloBarostat) XmlSerializer::deserialize<MonteCarloBarostat>;
    %template(XmlSerializer_deserialize_MonteCarloAnisotropicBarostat) XmlSerializer::deserialize<MonteCarloAnisotropicBarostat>;
    %template(XmlSerializer_deserialize_NonbondedForce) XmlSerializer::deserialize<NonbondedForce>;
    %template(XmlSerializer_deserialize_RBTorsionForce) XmlSerializer::deserialize<RBTorsionForce>;
    %template(XmlSerializer_deserialize_System) XmlSerializer::deserialize<System>;
};
*/
