
#include "OpenMMCWrapper.h"
#include "OpenMM.h"
#include <cstring>
#include <vector>

using namespace OpenMM;
using namespace std;

#if defined(__cplusplus)
extern "C" {
#endif

/* OpenMM_Vec3 */
void openmm_vec3_scale_(const OpenMM_Vec3& vec, double const& scale, OpenMM_Vec3& result) {
    result = OpenMM_Vec3_scale(vec, scale);
}
void OPENMM_VEC3_SCALE(const OpenMM_Vec3& vec, double const& scale, OpenMM_Vec3& result) {
    result = OpenMM_Vec3_scale(vec, scale);
}

/* OpenMM_Vec3Array */
void openmm_vec3array_create_(OpenMM_Vec3Array*& result, const int& size) {
    result = OpenMM_Vec3Array_create(size);
}
void OPENMM_VEC3ARRAY_CREATE(OpenMM_Vec3Array*& result, const int& size) {
    result = OpenMM_Vec3Array_create(size);
}
void openmm_vec3array_destroy_(OpenMM_Vec3Array*& array) {
    OpenMM_Vec3Array_destroy(array);
    array = 0;
}
void OPENMM_VEC3ARRAY_DESTROY(OpenMM_Vec3Array*& array) {
    OpenMM_Vec3Array_destroy(array);
    array = 0;
}
int openmm_vec3array_getsize_(const OpenMM_Vec3Array* const& array) {
    return OpenMM_Vec3Array_getSize(array);
}
int OPENMM_VEC3ARRAY_GETSIZE(const OpenMM_Vec3Array* const& array) {
    return OpenMM_Vec3Array_getSize(array);
}
void openmm_vec3array_resize_(OpenMM_Vec3Array* const& array, const int& size) {
    OpenMM_Vec3Array_resize(array, size);
}
void OPENMM_VEC3ARRAY_RESIZE(OpenMM_Vec3Array* const& array, const int& size) {
    OpenMM_Vec3Array_resize(array, size);
}
void openmm_vec3array_append_(OpenMM_Vec3Array* const& array, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_append(array, vec);
}
void OPENMM_VEC3ARRAY_APPEND(OpenMM_Vec3Array* const& array, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_append(array, vec);
}
void openmm_vec3array_set_(OpenMM_Vec3Array* const& array, const int& index, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_set(array, index-1, vec);
}
void OPENMM_VEC3ARRAY_SET(OpenMM_Vec3Array* const& array, const int& index, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_set(array, index-1, vec);
}
void openmm_vec3array_get_(const OpenMM_Vec3Array* const& array, const int& index, OpenMM_Vec3& result) {
    result = *OpenMM_Vec3Array_get(array, index-1);
}
void OPENMM_VEC3ARRAY_GET(const OpenMM_Vec3Array* const& array, const int& index, OpenMM_Vec3& result) {
    result = *OpenMM_Vec3Array_get(array, index-1);
}

/* OpenMM_StringArray */
void openmm_stringarray_create_(OpenMM_StringArray*& result, const int& size) {
    result = OpenMM_StringArray_create(size);
}
void OPENMM_STRINGARRAY_CREATE(OpenMM_StringArray*& result, const int& size) {
    result = OpenMM_StringArray_create(size);
}
void openmm_stringarray_destroy_(OpenMM_StringArray*& array) {
    OpenMM_StringArray_destroy(array);
    array = 0;
}
void OPENMM_STRINGARRAY_DESTROY(OpenMM_StringArray*& array) {
    OpenMM_StringArray_destroy(array);
    array = 0;
}
int openmm_stringarray_getsize_(const OpenMM_StringArray* const& array) {
    return OpenMM_StringArray_getSize(array);
}
int OPENMM_STRINGARRAY_GETSIZE(const OpenMM_StringArray* const& array) {
    return OpenMM_StringArray_getSize(array);
}
void openmm_stringarray_resize_(OpenMM_StringArray* const& array, const int& size) {
    OpenMM_StringArray_resize(array, size);
}
void OPENMM_STRINGARRAY_RESIZE(OpenMM_StringArray* const& array, const int& size) {
    OpenMM_StringArray_resize(array, size);
}
void openmm_stringarray_append_(OpenMM_StringArray* const& array, const char* str, int length) {
    OpenMM_StringArray_append(array, string(str, length).c_str());
}
void OPENMM_STRINGARRAY_APPEND(OpenMM_StringArray* const& array, const char* str, int length) {
    OpenMM_StringArray_append(array, string(str, length).c_str());
}
void openmm_stringarray_set_(OpenMM_StringArray* const& array, const int& index, const char* str, int length) {
    OpenMM_StringArray_set(array, index-1, string(str, length).c_str());
}
void OPENMM_STRINGARRAY_SET(OpenMM_StringArray* const& array, const int& index, const char* str, int length) {
    OpenMM_StringArray_set(array, index-1, string(str, length).c_str());
}
void openmm_stringarray_get_(const OpenMM_StringArray* const& array, const int& index, char* result, int length) {
    const char* str = OpenMM_StringArray_get(array, index-1);
    strncpy(result, str, length);
}
void OPENMM_STRINGARRAY_GET(const OpenMM_StringArray* const& array, const int& index, char* result, int length) {
    const char* str = OpenMM_StringArray_get(array, index-1);
    strncpy(result, str, length);
}

/* OpenMM_BondArray */
void openmm_bondarray_create_(OpenMM_BondArray*& result, const int& size) {
    result = OpenMM_BondArray_create(size);
}
void OPENMM_BONDARRAY_CREATE(OpenMM_BondArray*& result, const int& size) {
    result = OpenMM_BondArray_create(size);
}
void openmm_bondarray_destroy_(OpenMM_BondArray*& array) {
    OpenMM_BondArray_destroy(array);
    array = 0;
}
void OPENMM_BONDARRAY_DESTROY(OpenMM_BondArray*& array) {
    OpenMM_BondArray_destroy(array);
    array = 0;
}
int openmm_bondarray_getsize_(const OpenMM_BondArray* const& array) {
    return OpenMM_BondArray_getSize(array);
}
int OPENMM_BONDARRAY_GETSIZE(const OpenMM_BondArray* const& array) {
    return OpenMM_BondArray_getSize(array);
}
void openmm_bondarray_resize_(OpenMM_BondArray* const& array, const int& size) {
    OpenMM_BondArray_resize(array, size);
}
void OPENMM_BONDARRAY_RESIZE(OpenMM_BondArray* const& array, const int& size) {
    OpenMM_BondArray_resize(array, size);
}
void openmm_bondarray_append_(OpenMM_BondArray* const& array, const int& particle1, const int& particle2) {
    OpenMM_BondArray_append(array, particle1, particle2);
}
void OPENMM_BONDARRAY_APPEND(OpenMM_BondArray* const& array, const int& particle1, const int& particle2) {
    OpenMM_BondArray_append(array, particle1, particle2);
}
void openmm_bondarray_set_(OpenMM_BondArray* const& array, const int& index, const int& particle1, const int& particle2) {
    OpenMM_BondArray_set(array, index-1, particle1, particle2);
}
void OPENMM_BONDARRAY_SET(OpenMM_BondArray* const& array, const int& index, const int& particle1, const int& particle2) {
    OpenMM_BondArray_set(array, index-1, particle1, particle2);
}
void openmm_bondarray_get_(const OpenMM_BondArray* const& array, const int& index, int* particle1, int* particle2) {
    OpenMM_BondArray_get(array, index-1, particle1, particle2);
}
void OPENMM_BONDARRAY_GET(const OpenMM_BondArray* const& array, const int& index, int* particle1, int* particle2) {
    OpenMM_BondArray_get(array, index-1, particle1, particle2);
}

/* OpenMM_ParameterArray */
int openmm_parameterarray_getsize_(const OpenMM_ParameterArray* const& array) {
    return OpenMM_ParameterArray_getSize(array);
}
int OPENMM_PARAMETERARRAY_GETSIZE(const OpenMM_ParameterArray* const& array) {
    return OpenMM_ParameterArray_getSize(array);
}
double openmm_parameterarray_get_(const OpenMM_ParameterArray* const& array, const char* name, int length) {
    return OpenMM_ParameterArray_get(array, string(name, length).c_str());
}
double OPENMM_PARAMETERARRAY_GET(const OpenMM_ParameterArray* const& array, const char* name, int length) {
    return OpenMM_ParameterArray_get(array, string(name, length).c_str());
}

/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
void openmm_context_getstate_(const OpenMM_Context*& target, int const& types, OpenMM_State*& result) {
    result = OpenMM_Context_getState(target, types);
};
void OPENMM_CONTEXT_GETSTATE(const OpenMM_Context*& target, int const& types, OpenMM_State*& result) {
    result = OpenMM_Context_getState(target, types);
};
void openmm_platform_loadpluginsfromdirectory_(const char* directory, OpenMM_StringArray*& result, int length) {
    result = OpenMM_Platform_loadPluginsFromDirectory(string(directory, length).c_str());
};
void OPENMM_PLATFORM_LOADPLUGINSFROMDIRECTORY(const char* directory, OpenMM_StringArray*& result, int length) {
    result = OpenMM_Platform_loadPluginsFromDirectory(string(directory, length).c_str());
};

 
 

/* OpenMM::HarmonicBondForce*/
void openmm_harmonicbondforce_create_(OpenMM_HarmonicBondForce*& result) {
    result = OpenMM_HarmonicBondForce_create();
}
void OPENMM_HARMONICBONDFORCE_CREATE(OpenMM_HarmonicBondForce*& result) {
    result = OpenMM_HarmonicBondForce_create();
}
void openmm_harmonicbondforce_destroy_(OpenMM_HarmonicBondForce*& destroy) {
    OpenMM_HarmonicBondForce_destroy(destroy);
    destroy = 0;
}
void OPENMM_HARMONICBONDFORCE_DESTROY(OpenMM_HarmonicBondForce*& destroy) {
    OpenMM_HarmonicBondForce_destroy(destroy);
    destroy = 0;
}
int openmm_harmonicbondforce_getnumbonds_(const OpenMM_HarmonicBondForce*& target) {
    return OpenMM_HarmonicBondForce_getNumBonds(target);
};
int OPENMM_HARMONICBONDFORCE_GETNUMBONDS(const OpenMM_HarmonicBondForce*& target) {
    return OpenMM_HarmonicBondForce_getNumBonds(target);
};
int openmm_harmonicbondforce_addbond_(OpenMM_HarmonicBondForce*& target, int const& particle1, int const& particle2, double const& length, double const& k) {
    return OpenMM_HarmonicBondForce_addBond(target, particle1, particle2, length, k);
};
int OPENMM_HARMONICBONDFORCE_ADDBOND(OpenMM_HarmonicBondForce*& target, int const& particle1, int const& particle2, double const& length, double const& k) {
    return OpenMM_HarmonicBondForce_addBond(target, particle1, particle2, length, k);
};
void openmm_harmonicbondforce_getbondparameters_(const OpenMM_HarmonicBondForce*& target, int const& index, int* particle1, int* particle2, double* length, double* k) {
    OpenMM_HarmonicBondForce_getBondParameters(target, index, particle1, particle2, length, k);
};
void OPENMM_HARMONICBONDFORCE_GETBONDPARAMETERS(const OpenMM_HarmonicBondForce*& target, int const& index, int* particle1, int* particle2, double* length, double* k) {
    OpenMM_HarmonicBondForce_getBondParameters(target, index, particle1, particle2, length, k);
};
void openmm_harmonicbondforce_setbondparameters_(OpenMM_HarmonicBondForce*& target, int const& index, int const& particle1, int const& particle2, double const& length, double const& k) {
    OpenMM_HarmonicBondForce_setBondParameters(target, index, particle1, particle2, length, k);
};
void OPENMM_HARMONICBONDFORCE_SETBONDPARAMETERS(OpenMM_HarmonicBondForce*& target, int const& index, int const& particle1, int const& particle2, double const& length, double const& k) {
    OpenMM_HarmonicBondForce_setBondParameters(target, index, particle1, particle2, length, k);
};


/* OpenMM::BrownianIntegrator*/
void openmm_brownianintegrator_create_(OpenMM_BrownianIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_BrownianIntegrator_create(temperature, frictionCoeff, stepSize);
}
void OPENMM_BROWNIANINTEGRATOR_CREATE(OpenMM_BrownianIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_BrownianIntegrator_create(temperature, frictionCoeff, stepSize);
}
void openmm_brownianintegrator_destroy_(OpenMM_BrownianIntegrator*& destroy) {
    OpenMM_BrownianIntegrator_destroy(destroy);
    destroy = 0;
}
void OPENMM_BROWNIANINTEGRATOR_DESTROY(OpenMM_BrownianIntegrator*& destroy) {
    OpenMM_BrownianIntegrator_destroy(destroy);
    destroy = 0;
}
double openmm_brownianintegrator_gettemperature_(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getTemperature(target);
};
double OPENMM_BROWNIANINTEGRATOR_GETTEMPERATURE(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getTemperature(target);
};
void openmm_brownianintegrator_settemperature_(OpenMM_BrownianIntegrator*& target, double const& temp) {
    OpenMM_BrownianIntegrator_setTemperature(target, temp);
};
void OPENMM_BROWNIANINTEGRATOR_SETTEMPERATURE(OpenMM_BrownianIntegrator*& target, double const& temp) {
    OpenMM_BrownianIntegrator_setTemperature(target, temp);
};
double openmm_brownianintegrator_getfriction_(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getFriction(target);
};
double OPENMM_BROWNIANINTEGRATOR_GETFRICTION(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getFriction(target);
};
void openmm_brownianintegrator_setfriction_(OpenMM_BrownianIntegrator*& target, double const& coeff) {
    OpenMM_BrownianIntegrator_setFriction(target, coeff);
};
void OPENMM_BROWNIANINTEGRATOR_SETFRICTION(OpenMM_BrownianIntegrator*& target, double const& coeff) {
    OpenMM_BrownianIntegrator_setFriction(target, coeff);
};
int openmm_brownianintegrator_getrandomnumberseed_(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getRandomNumberSeed(target);
};
int OPENMM_BROWNIANINTEGRATOR_GETRANDOMNUMBERSEED(const OpenMM_BrownianIntegrator*& target) {
    return OpenMM_BrownianIntegrator_getRandomNumberSeed(target);
};
void openmm_brownianintegrator_setrandomnumberseed_(OpenMM_BrownianIntegrator*& target, int const& seed) {
    OpenMM_BrownianIntegrator_setRandomNumberSeed(target, seed);
};
void OPENMM_BROWNIANINTEGRATOR_SETRANDOMNUMBERSEED(OpenMM_BrownianIntegrator*& target, int const& seed) {
    OpenMM_BrownianIntegrator_setRandomNumberSeed(target, seed);
};
void openmm_brownianintegrator_step_(OpenMM_BrownianIntegrator*& target, int const& steps) {
    OpenMM_BrownianIntegrator_step(target, steps);
};
void OPENMM_BROWNIANINTEGRATOR_STEP(OpenMM_BrownianIntegrator*& target, int const& steps) {
    OpenMM_BrownianIntegrator_step(target, steps);
};


/* OpenMM::OpenMMException*/
void openmm_openmmexception_create_(OpenMM_OpenMMException*& result, const char* message, int message_length) {
    result = OpenMM_OpenMMException_create(string(message, message_length).c_str());
}
void OPENMM_OPENMMEXCEPTION_CREATE(OpenMM_OpenMMException*& result, const char* message, int message_length) {
    result = OpenMM_OpenMMException_create(string(message, message_length).c_str());
}
void openmm_openmmexception_destroy_(OpenMM_OpenMMException*& destroy) {
    OpenMM_OpenMMException_destroy(destroy);
    destroy = 0;
}
void OPENMM_OPENMMEXCEPTION_DESTROY(OpenMM_OpenMMException*& destroy) {
    OpenMM_OpenMMException_destroy(destroy);
    destroy = 0;
}
void openmm_openmmexception_what_(const OpenMM_OpenMMException*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_OpenMMException_what(target);
    strncpy(result, result_chars, result_length);
 
};
void OPENMM_OPENMMEXCEPTION_WHAT(const OpenMM_OpenMMException*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_OpenMMException_what(target);
    strncpy(result, result_chars, result_length);
 
};


/* OpenMM::NonbondedForce*/
void openmm_nonbondedforce_create_(OpenMM_NonbondedForce*& result) {
    result = OpenMM_NonbondedForce_create();
}
void OPENMM_NONBONDEDFORCE_CREATE(OpenMM_NonbondedForce*& result) {
    result = OpenMM_NonbondedForce_create();
}
void openmm_nonbondedforce_destroy_(OpenMM_NonbondedForce*& destroy) {
    OpenMM_NonbondedForce_destroy(destroy);
    destroy = 0;
}
void OPENMM_NONBONDEDFORCE_DESTROY(OpenMM_NonbondedForce*& destroy) {
    OpenMM_NonbondedForce_destroy(destroy);
    destroy = 0;
}
int openmm_nonbondedforce_getnumparticles_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumParticles(target);
};
int OPENMM_NONBONDEDFORCE_GETNUMPARTICLES(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumParticles(target);
};
int openmm_nonbondedforce_getnumexceptions_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumExceptions(target);
};
int OPENMM_NONBONDEDFORCE_GETNUMEXCEPTIONS(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getNumExceptions(target);
};
void openmm_nonbondedforce_getnonbondedmethod_(const OpenMM_NonbondedForce*& target, int& result) {
    result = OpenMM_NonbondedForce_getNonbondedMethod(target);
};
void OPENMM_NONBONDEDFORCE_GETNONBONDEDMETHOD(const OpenMM_NonbondedForce*& target, int& result) {
    result = OpenMM_NonbondedForce_getNonbondedMethod(target);
};
void openmm_nonbondedforce_setnonbondedmethod_(OpenMM_NonbondedForce*& target, int const& method) {
    OpenMM_NonbondedForce_setNonbondedMethod(target, (OpenMM_NonbondedForce_NonbondedMethod) method);
};
void OPENMM_NONBONDEDFORCE_SETNONBONDEDMETHOD(OpenMM_NonbondedForce*& target, int const& method) {
    OpenMM_NonbondedForce_setNonbondedMethod(target, (OpenMM_NonbondedForce_NonbondedMethod) method);
};
double openmm_nonbondedforce_getcutoffdistance_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getCutoffDistance(target);
};
double OPENMM_NONBONDEDFORCE_GETCUTOFFDISTANCE(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getCutoffDistance(target);
};
void openmm_nonbondedforce_setcutoffdistance_(OpenMM_NonbondedForce*& target, double const& distance) {
    OpenMM_NonbondedForce_setCutoffDistance(target, distance);
};
void OPENMM_NONBONDEDFORCE_SETCUTOFFDISTANCE(OpenMM_NonbondedForce*& target, double const& distance) {
    OpenMM_NonbondedForce_setCutoffDistance(target, distance);
};
double openmm_nonbondedforce_getreactionfielddielectric_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getReactionFieldDielectric(target);
};
double OPENMM_NONBONDEDFORCE_GETREACTIONFIELDDIELECTRIC(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getReactionFieldDielectric(target);
};
void openmm_nonbondedforce_setreactionfielddielectric_(OpenMM_NonbondedForce*& target, double const& dielectric) {
    OpenMM_NonbondedForce_setReactionFieldDielectric(target, dielectric);
};
void OPENMM_NONBONDEDFORCE_SETREACTIONFIELDDIELECTRIC(OpenMM_NonbondedForce*& target, double const& dielectric) {
    OpenMM_NonbondedForce_setReactionFieldDielectric(target, dielectric);
};
double openmm_nonbondedforce_getewalderrortolerance_(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getEwaldErrorTolerance(target);
};
double OPENMM_NONBONDEDFORCE_GETEWALDERRORTOLERANCE(const OpenMM_NonbondedForce*& target) {
    return OpenMM_NonbondedForce_getEwaldErrorTolerance(target);
};
void openmm_nonbondedforce_setewalderrortolerance_(OpenMM_NonbondedForce*& target, double const& tol) {
    OpenMM_NonbondedForce_setEwaldErrorTolerance(target, tol);
};
void OPENMM_NONBONDEDFORCE_SETEWALDERRORTOLERANCE(OpenMM_NonbondedForce*& target, double const& tol) {
    OpenMM_NonbondedForce_setEwaldErrorTolerance(target, tol);
};
void openmm_nonbondedforce_getperiodicboxvectors_(const OpenMM_NonbondedForce*& target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c) {
    OpenMM_NonbondedForce_getPeriodicBoxVectors(target, a, b, c);
};
void OPENMM_NONBONDEDFORCE_GETPERIODICBOXVECTORS(const OpenMM_NonbondedForce*& target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c) {
    OpenMM_NonbondedForce_getPeriodicBoxVectors(target, a, b, c);
};
void openmm_nonbondedforce_setperiodicboxvectors_(OpenMM_NonbondedForce*& target, OpenMM_Vec3 const& a, OpenMM_Vec3 const& b, OpenMM_Vec3 const& c) {
    OpenMM_NonbondedForce_setPeriodicBoxVectors(target, a, b, c);
};
void OPENMM_NONBONDEDFORCE_SETPERIODICBOXVECTORS(OpenMM_NonbondedForce*& target, OpenMM_Vec3 const& a, OpenMM_Vec3 const& b, OpenMM_Vec3 const& c) {
    OpenMM_NonbondedForce_setPeriodicBoxVectors(target, a, b, c);
};
int openmm_nonbondedforce_addparticle_(OpenMM_NonbondedForce*& target, double const& charge, double const& sigma, double const& epsilon) {
    return OpenMM_NonbondedForce_addParticle(target, charge, sigma, epsilon);
};
int OPENMM_NONBONDEDFORCE_ADDPARTICLE(OpenMM_NonbondedForce*& target, double const& charge, double const& sigma, double const& epsilon) {
    return OpenMM_NonbondedForce_addParticle(target, charge, sigma, epsilon);
};
void openmm_nonbondedforce_getparticleparameters_(const OpenMM_NonbondedForce*& target, int const& index, double* charge, double* sigma, double* epsilon) {
    OpenMM_NonbondedForce_getParticleParameters(target, index, charge, sigma, epsilon);
};
void OPENMM_NONBONDEDFORCE_GETPARTICLEPARAMETERS(const OpenMM_NonbondedForce*& target, int const& index, double* charge, double* sigma, double* epsilon) {
    OpenMM_NonbondedForce_getParticleParameters(target, index, charge, sigma, epsilon);
};
void openmm_nonbondedforce_setparticleparameters_(OpenMM_NonbondedForce*& target, int const& index, double const& charge, double const& sigma, double const& epsilon) {
    OpenMM_NonbondedForce_setParticleParameters(target, index, charge, sigma, epsilon);
};
void OPENMM_NONBONDEDFORCE_SETPARTICLEPARAMETERS(OpenMM_NonbondedForce*& target, int const& index, double const& charge, double const& sigma, double const& epsilon) {
    OpenMM_NonbondedForce_setParticleParameters(target, index, charge, sigma, epsilon);
};
int openmm_nonbondedforce_addexception_(OpenMM_NonbondedForce*& target, int const& particle1, int const& particle2, double const& chargeProd, double const& sigma, double const& epsilon, OpenMM_Boolean const& replace) {
    return OpenMM_NonbondedForce_addException(target, particle1, particle2, chargeProd, sigma, epsilon, replace);
};
int OPENMM_NONBONDEDFORCE_ADDEXCEPTION(OpenMM_NonbondedForce*& target, int const& particle1, int const& particle2, double const& chargeProd, double const& sigma, double const& epsilon, OpenMM_Boolean const& replace) {
    return OpenMM_NonbondedForce_addException(target, particle1, particle2, chargeProd, sigma, epsilon, replace);
};
void openmm_nonbondedforce_getexceptionparameters_(const OpenMM_NonbondedForce*& target, int const& index, int* particle1, int* particle2, double* chargeProd, double* sigma, double* epsilon) {
    OpenMM_NonbondedForce_getExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon);
};
void OPENMM_NONBONDEDFORCE_GETEXCEPTIONPARAMETERS(const OpenMM_NonbondedForce*& target, int const& index, int* particle1, int* particle2, double* chargeProd, double* sigma, double* epsilon) {
    OpenMM_NonbondedForce_getExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon);
};
void openmm_nonbondedforce_setexceptionparameters_(OpenMM_NonbondedForce*& target, int const& index, int const& particle1, int const& particle2, double const& chargeProd, double const& sigma, double const& epsilon) {
    OpenMM_NonbondedForce_setExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon);
};
void OPENMM_NONBONDEDFORCE_SETEXCEPTIONPARAMETERS(OpenMM_NonbondedForce*& target, int const& index, int const& particle1, int const& particle2, double const& chargeProd, double const& sigma, double const& epsilon) {
    OpenMM_NonbondedForce_setExceptionParameters(target, index, particle1, particle2, chargeProd, sigma, epsilon);
};
void openmm_nonbondedforce_createexceptionsfrombonds_(OpenMM_NonbondedForce*& target, const OpenMM_BondArray*& bonds, double const& coulomb14Scale, double const& lj14Scale) {
    OpenMM_NonbondedForce_createExceptionsFromBonds(target, bonds, coulomb14Scale, lj14Scale);
};
void OPENMM_NONBONDEDFORCE_CREATEEXCEPTIONSFROMBONDS(OpenMM_NonbondedForce*& target, const OpenMM_BondArray*& bonds, double const& coulomb14Scale, double const& lj14Scale) {
    OpenMM_NonbondedForce_createExceptionsFromBonds(target, bonds, coulomb14Scale, lj14Scale);
};


/* OpenMM::VariableLangevinIntegrator*/
void openmm_variablelangevinintegrator_create_(OpenMM_VariableLangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& errorTol) {
    result = OpenMM_VariableLangevinIntegrator_create(temperature, frictionCoeff, errorTol);
}
void OPENMM_VARIABLELANGEVININTEGRATOR_CREATE(OpenMM_VariableLangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& errorTol) {
    result = OpenMM_VariableLangevinIntegrator_create(temperature, frictionCoeff, errorTol);
}
void openmm_variablelangevinintegrator_destroy_(OpenMM_VariableLangevinIntegrator*& destroy) {
    OpenMM_VariableLangevinIntegrator_destroy(destroy);
    destroy = 0;
}
void OPENMM_VARIABLELANGEVININTEGRATOR_DESTROY(OpenMM_VariableLangevinIntegrator*& destroy) {
    OpenMM_VariableLangevinIntegrator_destroy(destroy);
    destroy = 0;
}
double openmm_variablelangevinintegrator_gettemperature_(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getTemperature(target);
};
double OPENMM_VARIABLELANGEVININTEGRATOR_GETTEMPERATURE(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getTemperature(target);
};
void openmm_variablelangevinintegrator_settemperature_(OpenMM_VariableLangevinIntegrator*& target, double const& temp) {
    OpenMM_VariableLangevinIntegrator_setTemperature(target, temp);
};
void OPENMM_VARIABLELANGEVININTEGRATOR_SETTEMPERATURE(OpenMM_VariableLangevinIntegrator*& target, double const& temp) {
    OpenMM_VariableLangevinIntegrator_setTemperature(target, temp);
};
double openmm_variablelangevinintegrator_getfriction_(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getFriction(target);
};
double OPENMM_VARIABLELANGEVININTEGRATOR_GETFRICTION(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getFriction(target);
};
void openmm_variablelangevinintegrator_setfriction_(OpenMM_VariableLangevinIntegrator*& target, double const& coeff) {
    OpenMM_VariableLangevinIntegrator_setFriction(target, coeff);
};
void OPENMM_VARIABLELANGEVININTEGRATOR_SETFRICTION(OpenMM_VariableLangevinIntegrator*& target, double const& coeff) {
    OpenMM_VariableLangevinIntegrator_setFriction(target, coeff);
};
double openmm_variablelangevinintegrator_geterrortolerance_(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getErrorTolerance(target);
};
double OPENMM_VARIABLELANGEVININTEGRATOR_GETERRORTOLERANCE(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getErrorTolerance(target);
};
void openmm_variablelangevinintegrator_seterrortolerance_(OpenMM_VariableLangevinIntegrator*& target, double const& tol) {
    OpenMM_VariableLangevinIntegrator_setErrorTolerance(target, tol);
};
void OPENMM_VARIABLELANGEVININTEGRATOR_SETERRORTOLERANCE(OpenMM_VariableLangevinIntegrator*& target, double const& tol) {
    OpenMM_VariableLangevinIntegrator_setErrorTolerance(target, tol);
};
int openmm_variablelangevinintegrator_getrandomnumberseed_(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getRandomNumberSeed(target);
};
int OPENMM_VARIABLELANGEVININTEGRATOR_GETRANDOMNUMBERSEED(const OpenMM_VariableLangevinIntegrator*& target) {
    return OpenMM_VariableLangevinIntegrator_getRandomNumberSeed(target);
};
void openmm_variablelangevinintegrator_setrandomnumberseed_(OpenMM_VariableLangevinIntegrator*& target, int const& seed) {
    OpenMM_VariableLangevinIntegrator_setRandomNumberSeed(target, seed);
};
void OPENMM_VARIABLELANGEVININTEGRATOR_SETRANDOMNUMBERSEED(OpenMM_VariableLangevinIntegrator*& target, int const& seed) {
    OpenMM_VariableLangevinIntegrator_setRandomNumberSeed(target, seed);
};
void openmm_variablelangevinintegrator_step_(OpenMM_VariableLangevinIntegrator*& target, int const& steps) {
    OpenMM_VariableLangevinIntegrator_step(target, steps);
};
void OPENMM_VARIABLELANGEVININTEGRATOR_STEP(OpenMM_VariableLangevinIntegrator*& target, int const& steps) {
    OpenMM_VariableLangevinIntegrator_step(target, steps);
};
void openmm_variablelangevinintegrator_stepto_(OpenMM_VariableLangevinIntegrator*& target, double const& time) {
    OpenMM_VariableLangevinIntegrator_stepTo(target, time);
};
void OPENMM_VARIABLELANGEVININTEGRATOR_STEPTO(OpenMM_VariableLangevinIntegrator*& target, double const& time) {
    OpenMM_VariableLangevinIntegrator_stepTo(target, time);
};


/* OpenMM::GBVIForce*/
void openmm_gbviforce_create_(OpenMM_GBVIForce*& result) {
    result = OpenMM_GBVIForce_create();
}
void OPENMM_GBVIFORCE_CREATE(OpenMM_GBVIForce*& result) {
    result = OpenMM_GBVIForce_create();
}
void openmm_gbviforce_destroy_(OpenMM_GBVIForce*& destroy) {
    OpenMM_GBVIForce_destroy(destroy);
    destroy = 0;
}
void OPENMM_GBVIFORCE_DESTROY(OpenMM_GBVIForce*& destroy) {
    OpenMM_GBVIForce_destroy(destroy);
    destroy = 0;
}
int openmm_gbviforce_getnumparticles_(const OpenMM_GBVIForce*& target) {
    return OpenMM_GBVIForce_getNumParticles(target);
};
int OPENMM_GBVIFORCE_GETNUMPARTICLES(const OpenMM_GBVIForce*& target) {
    return OpenMM_GBVIForce_getNumParticles(target);
};
int openmm_gbviforce_addparticle_(OpenMM_GBVIForce*& target, double const& charge, double const& radius, double const& gamma) {
    return OpenMM_GBVIForce_addParticle(target, charge, radius, gamma);
};
int OPENMM_GBVIFORCE_ADDPARTICLE(OpenMM_GBVIForce*& target, double const& charge, double const& radius, double const& gamma) {
    return OpenMM_GBVIForce_addParticle(target, charge, radius, gamma);
};
void openmm_gbviforce_getparticleparameters_(const OpenMM_GBVIForce*& target, int const& index, double* charge, double* radius, double* gamma) {
    OpenMM_GBVIForce_getParticleParameters(target, index, charge, radius, gamma);
};
void OPENMM_GBVIFORCE_GETPARTICLEPARAMETERS(const OpenMM_GBVIForce*& target, int const& index, double* charge, double* radius, double* gamma) {
    OpenMM_GBVIForce_getParticleParameters(target, index, charge, radius, gamma);
};
void openmm_gbviforce_setparticleparameters_(OpenMM_GBVIForce*& target, int const& index, double const& charge, double const& radius, double const& gamma) {
    OpenMM_GBVIForce_setParticleParameters(target, index, charge, radius, gamma);
};
void OPENMM_GBVIFORCE_SETPARTICLEPARAMETERS(OpenMM_GBVIForce*& target, int const& index, double const& charge, double const& radius, double const& gamma) {
    OpenMM_GBVIForce_setParticleParameters(target, index, charge, radius, gamma);
};
double openmm_gbviforce_getsolventdielectric_(const OpenMM_GBVIForce*& target) {
    return OpenMM_GBVIForce_getSolventDielectric(target);
};
double OPENMM_GBVIFORCE_GETSOLVENTDIELECTRIC(const OpenMM_GBVIForce*& target) {
    return OpenMM_GBVIForce_getSolventDielectric(target);
};
void openmm_gbviforce_setsolventdielectric_(OpenMM_GBVIForce*& target, double const& dielectric) {
    OpenMM_GBVIForce_setSolventDielectric(target, dielectric);
};
void OPENMM_GBVIFORCE_SETSOLVENTDIELECTRIC(OpenMM_GBVIForce*& target, double const& dielectric) {
    OpenMM_GBVIForce_setSolventDielectric(target, dielectric);
};
double openmm_gbviforce_getsolutedielectric_(const OpenMM_GBVIForce*& target) {
    return OpenMM_GBVIForce_getSoluteDielectric(target);
};
double OPENMM_GBVIFORCE_GETSOLUTEDIELECTRIC(const OpenMM_GBVIForce*& target) {
    return OpenMM_GBVIForce_getSoluteDielectric(target);
};
void openmm_gbviforce_setsolutedielectric_(OpenMM_GBVIForce*& target, double const& dielectric) {
    OpenMM_GBVIForce_setSoluteDielectric(target, dielectric);
};
void OPENMM_GBVIFORCE_SETSOLUTEDIELECTRIC(OpenMM_GBVIForce*& target, double const& dielectric) {
    OpenMM_GBVIForce_setSoluteDielectric(target, dielectric);
};


/* OpenMM::Context*/
void openmm_context_create_(OpenMM_Context*& result, OpenMM_System*& system, OpenMM_Integrator*& integrator) {
    result = OpenMM_Context_create(system, integrator);
}
void OPENMM_CONTEXT_CREATE(OpenMM_Context*& result, OpenMM_System*& system, OpenMM_Integrator*& integrator) {
    result = OpenMM_Context_create(system, integrator);
}
void openmm_context_create_2_(OpenMM_Context*& result, OpenMM_System*& system, OpenMM_Integrator*& integrator, OpenMM_Platform*& platform) {
    result = OpenMM_Context_create_2(system, integrator, platform);
}
void OPENMM_CONTEXT_CREATE_2(OpenMM_Context*& result, OpenMM_System*& system, OpenMM_Integrator*& integrator, OpenMM_Platform*& platform) {
    result = OpenMM_Context_create_2(system, integrator, platform);
}
void openmm_context_destroy_(OpenMM_Context*& destroy) {
    OpenMM_Context_destroy(destroy);
    destroy = 0;
}
void OPENMM_CONTEXT_DESTROY(OpenMM_Context*& destroy) {
    OpenMM_Context_destroy(destroy);
    destroy = 0;
}
void openmm_context_getsystem_(OpenMM_Context*& target, OpenMM_System*& result) {
    result = OpenMM_Context_getSystem(target);
};
void OPENMM_CONTEXT_GETSYSTEM(OpenMM_Context*& target, OpenMM_System*& result) {
    result = OpenMM_Context_getSystem(target);
};
void openmm_context_getintegrator_(OpenMM_Context*& target, OpenMM_Integrator*& result) {
    result = OpenMM_Context_getIntegrator(target);
};
void OPENMM_CONTEXT_GETINTEGRATOR(OpenMM_Context*& target, OpenMM_Integrator*& result) {
    result = OpenMM_Context_getIntegrator(target);
};
void openmm_context_getplatform_(OpenMM_Context*& target, OpenMM_Platform*& result) {
    result = OpenMM_Context_getPlatform(target);
};
void OPENMM_CONTEXT_GETPLATFORM(OpenMM_Context*& target, OpenMM_Platform*& result) {
    result = OpenMM_Context_getPlatform(target);
};
void openmm_context_settime_(OpenMM_Context*& target, double const& time) {
    OpenMM_Context_setTime(target, time);
};
void OPENMM_CONTEXT_SETTIME(OpenMM_Context*& target, double const& time) {
    OpenMM_Context_setTime(target, time);
};
void openmm_context_setpositions_(OpenMM_Context*& target, const OpenMM_Vec3Array*& positions) {
    OpenMM_Context_setPositions(target, positions);
};
void OPENMM_CONTEXT_SETPOSITIONS(OpenMM_Context*& target, const OpenMM_Vec3Array*& positions) {
    OpenMM_Context_setPositions(target, positions);
};
void openmm_context_setvelocities_(OpenMM_Context*& target, const OpenMM_Vec3Array*& velocities) {
    OpenMM_Context_setVelocities(target, velocities);
};
void OPENMM_CONTEXT_SETVELOCITIES(OpenMM_Context*& target, const OpenMM_Vec3Array*& velocities) {
    OpenMM_Context_setVelocities(target, velocities);
};
double openmm_context_getparameter_(OpenMM_Context*& target, const char* name, int name_length) {
    return OpenMM_Context_getParameter(target, string(name, name_length).c_str());
};
double OPENMM_CONTEXT_GETPARAMETER(OpenMM_Context*& target, const char* name, int name_length) {
    return OpenMM_Context_getParameter(target, string(name, name_length).c_str());
};
void openmm_context_setparameter_(OpenMM_Context*& target, const char* name, double const& value, int name_length) {
    OpenMM_Context_setParameter(target, string(name, name_length).c_str(), value);
};
void OPENMM_CONTEXT_SETPARAMETER(OpenMM_Context*& target, const char* name, double const& value, int name_length) {
    OpenMM_Context_setParameter(target, string(name, name_length).c_str(), value);
};
void openmm_context_reinitialize_(OpenMM_Context*& target) {
    OpenMM_Context_reinitialize(target);
};
void OPENMM_CONTEXT_REINITIALIZE(OpenMM_Context*& target) {
    OpenMM_Context_reinitialize(target);
};


/* OpenMM::GBSAOBCForce*/
void openmm_gbsaobcforce_create_(OpenMM_GBSAOBCForce*& result) {
    result = OpenMM_GBSAOBCForce_create();
}
void OPENMM_GBSAOBCFORCE_CREATE(OpenMM_GBSAOBCForce*& result) {
    result = OpenMM_GBSAOBCForce_create();
}
void openmm_gbsaobcforce_destroy_(OpenMM_GBSAOBCForce*& destroy) {
    OpenMM_GBSAOBCForce_destroy(destroy);
    destroy = 0;
}
void OPENMM_GBSAOBCFORCE_DESTROY(OpenMM_GBSAOBCForce*& destroy) {
    OpenMM_GBSAOBCForce_destroy(destroy);
    destroy = 0;
}
int openmm_gbsaobcforce_getnumparticles_(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getNumParticles(target);
};
int OPENMM_GBSAOBCFORCE_GETNUMPARTICLES(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getNumParticles(target);
};
int openmm_gbsaobcforce_addparticle_(OpenMM_GBSAOBCForce*& target, double const& charge, double const& radius, double const& scalingFactor) {
    return OpenMM_GBSAOBCForce_addParticle(target, charge, radius, scalingFactor);
};
int OPENMM_GBSAOBCFORCE_ADDPARTICLE(OpenMM_GBSAOBCForce*& target, double const& charge, double const& radius, double const& scalingFactor) {
    return OpenMM_GBSAOBCForce_addParticle(target, charge, radius, scalingFactor);
};
void openmm_gbsaobcforce_getparticleparameters_(const OpenMM_GBSAOBCForce*& target, int const& index, double* charge, double* radius, double* scalingFactor) {
    OpenMM_GBSAOBCForce_getParticleParameters(target, index, charge, radius, scalingFactor);
};
void OPENMM_GBSAOBCFORCE_GETPARTICLEPARAMETERS(const OpenMM_GBSAOBCForce*& target, int const& index, double* charge, double* radius, double* scalingFactor) {
    OpenMM_GBSAOBCForce_getParticleParameters(target, index, charge, radius, scalingFactor);
};
void openmm_gbsaobcforce_setparticleparameters_(OpenMM_GBSAOBCForce*& target, int const& index, double const& charge, double const& radius, double const& scalingFactor) {
    OpenMM_GBSAOBCForce_setParticleParameters(target, index, charge, radius, scalingFactor);
};
void OPENMM_GBSAOBCFORCE_SETPARTICLEPARAMETERS(OpenMM_GBSAOBCForce*& target, int const& index, double const& charge, double const& radius, double const& scalingFactor) {
    OpenMM_GBSAOBCForce_setParticleParameters(target, index, charge, radius, scalingFactor);
};
double openmm_gbsaobcforce_getsolventdielectric_(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSolventDielectric(target);
};
double OPENMM_GBSAOBCFORCE_GETSOLVENTDIELECTRIC(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSolventDielectric(target);
};
void openmm_gbsaobcforce_setsolventdielectric_(OpenMM_GBSAOBCForce*& target, double const& dielectric) {
    OpenMM_GBSAOBCForce_setSolventDielectric(target, dielectric);
};
void OPENMM_GBSAOBCFORCE_SETSOLVENTDIELECTRIC(OpenMM_GBSAOBCForce*& target, double const& dielectric) {
    OpenMM_GBSAOBCForce_setSolventDielectric(target, dielectric);
};
double openmm_gbsaobcforce_getsolutedielectric_(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSoluteDielectric(target);
};
double OPENMM_GBSAOBCFORCE_GETSOLUTEDIELECTRIC(const OpenMM_GBSAOBCForce*& target) {
    return OpenMM_GBSAOBCForce_getSoluteDielectric(target);
};
void openmm_gbsaobcforce_setsolutedielectric_(OpenMM_GBSAOBCForce*& target, double const& dielectric) {
    OpenMM_GBSAOBCForce_setSoluteDielectric(target, dielectric);
};
void OPENMM_GBSAOBCFORCE_SETSOLUTEDIELECTRIC(OpenMM_GBSAOBCForce*& target, double const& dielectric) {
    OpenMM_GBSAOBCForce_setSoluteDielectric(target, dielectric);
};


/* OpenMM::VariableVerletIntegrator*/
void openmm_variableverletintegrator_create_(OpenMM_VariableVerletIntegrator*& result, double const& errorTol) {
    result = OpenMM_VariableVerletIntegrator_create(errorTol);
}
void OPENMM_VARIABLEVERLETINTEGRATOR_CREATE(OpenMM_VariableVerletIntegrator*& result, double const& errorTol) {
    result = OpenMM_VariableVerletIntegrator_create(errorTol);
}
void openmm_variableverletintegrator_destroy_(OpenMM_VariableVerletIntegrator*& destroy) {
    OpenMM_VariableVerletIntegrator_destroy(destroy);
    destroy = 0;
}
void OPENMM_VARIABLEVERLETINTEGRATOR_DESTROY(OpenMM_VariableVerletIntegrator*& destroy) {
    OpenMM_VariableVerletIntegrator_destroy(destroy);
    destroy = 0;
}
double openmm_variableverletintegrator_geterrortolerance_(const OpenMM_VariableVerletIntegrator*& target) {
    return OpenMM_VariableVerletIntegrator_getErrorTolerance(target);
};
double OPENMM_VARIABLEVERLETINTEGRATOR_GETERRORTOLERANCE(const OpenMM_VariableVerletIntegrator*& target) {
    return OpenMM_VariableVerletIntegrator_getErrorTolerance(target);
};
void openmm_variableverletintegrator_seterrortolerance_(OpenMM_VariableVerletIntegrator*& target, double const& tol) {
    OpenMM_VariableVerletIntegrator_setErrorTolerance(target, tol);
};
void OPENMM_VARIABLEVERLETINTEGRATOR_SETERRORTOLERANCE(OpenMM_VariableVerletIntegrator*& target, double const& tol) {
    OpenMM_VariableVerletIntegrator_setErrorTolerance(target, tol);
};
void openmm_variableverletintegrator_step_(OpenMM_VariableVerletIntegrator*& target, int const& steps) {
    OpenMM_VariableVerletIntegrator_step(target, steps);
};
void OPENMM_VARIABLEVERLETINTEGRATOR_STEP(OpenMM_VariableVerletIntegrator*& target, int const& steps) {
    OpenMM_VariableVerletIntegrator_step(target, steps);
};
void openmm_variableverletintegrator_stepto_(OpenMM_VariableVerletIntegrator*& target, double const& time) {
    OpenMM_VariableVerletIntegrator_stepTo(target, time);
};
void OPENMM_VARIABLEVERLETINTEGRATOR_STEPTO(OpenMM_VariableVerletIntegrator*& target, double const& time) {
    OpenMM_VariableVerletIntegrator_stepTo(target, time);
};


/* OpenMM::CMMotionRemover*/
void openmm_cmmotionremover_create_(OpenMM_CMMotionRemover*& result, int const& frequency) {
    result = OpenMM_CMMotionRemover_create(frequency);
}
void OPENMM_CMMOTIONREMOVER_CREATE(OpenMM_CMMotionRemover*& result, int const& frequency) {
    result = OpenMM_CMMotionRemover_create(frequency);
}
void openmm_cmmotionremover_destroy_(OpenMM_CMMotionRemover*& destroy) {
    OpenMM_CMMotionRemover_destroy(destroy);
    destroy = 0;
}
void OPENMM_CMMOTIONREMOVER_DESTROY(OpenMM_CMMotionRemover*& destroy) {
    OpenMM_CMMotionRemover_destroy(destroy);
    destroy = 0;
}
int openmm_cmmotionremover_getfrequency_(const OpenMM_CMMotionRemover*& target) {
    return OpenMM_CMMotionRemover_getFrequency(target);
};
int OPENMM_CMMOTIONREMOVER_GETFREQUENCY(const OpenMM_CMMotionRemover*& target) {
    return OpenMM_CMMotionRemover_getFrequency(target);
};
void openmm_cmmotionremover_setfrequency_(OpenMM_CMMotionRemover*& target, int const& freq) {
    OpenMM_CMMotionRemover_setFrequency(target, freq);
};
void OPENMM_CMMOTIONREMOVER_SETFREQUENCY(OpenMM_CMMotionRemover*& target, int const& freq) {
    OpenMM_CMMotionRemover_setFrequency(target, freq);
};


/* OpenMM::VerletIntegrator*/
void openmm_verletintegrator_create_(OpenMM_VerletIntegrator*& result, double const& stepSize) {
    result = OpenMM_VerletIntegrator_create(stepSize);
}
void OPENMM_VERLETINTEGRATOR_CREATE(OpenMM_VerletIntegrator*& result, double const& stepSize) {
    result = OpenMM_VerletIntegrator_create(stepSize);
}
void openmm_verletintegrator_destroy_(OpenMM_VerletIntegrator*& destroy) {
    OpenMM_VerletIntegrator_destroy(destroy);
    destroy = 0;
}
void OPENMM_VERLETINTEGRATOR_DESTROY(OpenMM_VerletIntegrator*& destroy) {
    OpenMM_VerletIntegrator_destroy(destroy);
    destroy = 0;
}
void openmm_verletintegrator_step_(OpenMM_VerletIntegrator*& target, int const& steps) {
    OpenMM_VerletIntegrator_step(target, steps);
};
void OPENMM_VERLETINTEGRATOR_STEP(OpenMM_VerletIntegrator*& target, int const& steps) {
    OpenMM_VerletIntegrator_step(target, steps);
};


/* OpenMM::RBTorsionForce*/
void openmm_rbtorsionforce_create_(OpenMM_RBTorsionForce*& result) {
    result = OpenMM_RBTorsionForce_create();
}
void OPENMM_RBTORSIONFORCE_CREATE(OpenMM_RBTorsionForce*& result) {
    result = OpenMM_RBTorsionForce_create();
}
void openmm_rbtorsionforce_destroy_(OpenMM_RBTorsionForce*& destroy) {
    OpenMM_RBTorsionForce_destroy(destroy);
    destroy = 0;
}
void OPENMM_RBTORSIONFORCE_DESTROY(OpenMM_RBTorsionForce*& destroy) {
    OpenMM_RBTorsionForce_destroy(destroy);
    destroy = 0;
}
int openmm_rbtorsionforce_getnumtorsions_(const OpenMM_RBTorsionForce*& target) {
    return OpenMM_RBTorsionForce_getNumTorsions(target);
};
int OPENMM_RBTORSIONFORCE_GETNUMTORSIONS(const OpenMM_RBTorsionForce*& target) {
    return OpenMM_RBTorsionForce_getNumTorsions(target);
};
int openmm_rbtorsionforce_addtorsion_(OpenMM_RBTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5) {
    return OpenMM_RBTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
};
int OPENMM_RBTORSIONFORCE_ADDTORSION(OpenMM_RBTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5) {
    return OpenMM_RBTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
};
void openmm_rbtorsionforce_gettorsionparameters_(const OpenMM_RBTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, double* c0, double* c1, double* c2, double* c3, double* c4, double* c5) {
    OpenMM_RBTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
};
void OPENMM_RBTORSIONFORCE_GETTORSIONPARAMETERS(const OpenMM_RBTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, double* c0, double* c1, double* c2, double* c3, double* c4, double* c5) {
    OpenMM_RBTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
};
void openmm_rbtorsionforce_settorsionparameters_(OpenMM_RBTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5) {
    OpenMM_RBTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
};
void OPENMM_RBTORSIONFORCE_SETTORSIONPARAMETERS(OpenMM_RBTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, double const& c0, double const& c1, double const& c2, double const& c3, double const& c4, double const& c5) {
    OpenMM_RBTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
};


/* OpenMM::LangevinIntegrator*/
void openmm_langevinintegrator_create_(OpenMM_LangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_LangevinIntegrator_create(temperature, frictionCoeff, stepSize);
}
void OPENMM_LANGEVININTEGRATOR_CREATE(OpenMM_LangevinIntegrator*& result, double const& temperature, double const& frictionCoeff, double const& stepSize) {
    result = OpenMM_LangevinIntegrator_create(temperature, frictionCoeff, stepSize);
}
void openmm_langevinintegrator_destroy_(OpenMM_LangevinIntegrator*& destroy) {
    OpenMM_LangevinIntegrator_destroy(destroy);
    destroy = 0;
}
void OPENMM_LANGEVININTEGRATOR_DESTROY(OpenMM_LangevinIntegrator*& destroy) {
    OpenMM_LangevinIntegrator_destroy(destroy);
    destroy = 0;
}
double openmm_langevinintegrator_gettemperature_(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getTemperature(target);
};
double OPENMM_LANGEVININTEGRATOR_GETTEMPERATURE(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getTemperature(target);
};
void openmm_langevinintegrator_settemperature_(OpenMM_LangevinIntegrator*& target, double const& temp) {
    OpenMM_LangevinIntegrator_setTemperature(target, temp);
};
void OPENMM_LANGEVININTEGRATOR_SETTEMPERATURE(OpenMM_LangevinIntegrator*& target, double const& temp) {
    OpenMM_LangevinIntegrator_setTemperature(target, temp);
};
double openmm_langevinintegrator_getfriction_(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getFriction(target);
};
double OPENMM_LANGEVININTEGRATOR_GETFRICTION(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getFriction(target);
};
void openmm_langevinintegrator_setfriction_(OpenMM_LangevinIntegrator*& target, double const& coeff) {
    OpenMM_LangevinIntegrator_setFriction(target, coeff);
};
void OPENMM_LANGEVININTEGRATOR_SETFRICTION(OpenMM_LangevinIntegrator*& target, double const& coeff) {
    OpenMM_LangevinIntegrator_setFriction(target, coeff);
};
int openmm_langevinintegrator_getrandomnumberseed_(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getRandomNumberSeed(target);
};
int OPENMM_LANGEVININTEGRATOR_GETRANDOMNUMBERSEED(const OpenMM_LangevinIntegrator*& target) {
    return OpenMM_LangevinIntegrator_getRandomNumberSeed(target);
};
void openmm_langevinintegrator_setrandomnumberseed_(OpenMM_LangevinIntegrator*& target, int const& seed) {
    OpenMM_LangevinIntegrator_setRandomNumberSeed(target, seed);
};
void OPENMM_LANGEVININTEGRATOR_SETRANDOMNUMBERSEED(OpenMM_LangevinIntegrator*& target, int const& seed) {
    OpenMM_LangevinIntegrator_setRandomNumberSeed(target, seed);
};
void openmm_langevinintegrator_step_(OpenMM_LangevinIntegrator*& target, int const& steps) {
    OpenMM_LangevinIntegrator_step(target, steps);
};
void OPENMM_LANGEVININTEGRATOR_STEP(OpenMM_LangevinIntegrator*& target, int const& steps) {
    OpenMM_LangevinIntegrator_step(target, steps);
};


/* OpenMM::Force*/
void openmm_force_destroy_(OpenMM_Force*& destroy) {
    OpenMM_Force_destroy(destroy);
    destroy = 0;
}
void OPENMM_FORCE_DESTROY(OpenMM_Force*& destroy) {
    OpenMM_Force_destroy(destroy);
    destroy = 0;
}


/* OpenMM::HarmonicAngleForce*/
void openmm_harmonicangleforce_create_(OpenMM_HarmonicAngleForce*& result) {
    result = OpenMM_HarmonicAngleForce_create();
}
void OPENMM_HARMONICANGLEFORCE_CREATE(OpenMM_HarmonicAngleForce*& result) {
    result = OpenMM_HarmonicAngleForce_create();
}
void openmm_harmonicangleforce_destroy_(OpenMM_HarmonicAngleForce*& destroy) {
    OpenMM_HarmonicAngleForce_destroy(destroy);
    destroy = 0;
}
void OPENMM_HARMONICANGLEFORCE_DESTROY(OpenMM_HarmonicAngleForce*& destroy) {
    OpenMM_HarmonicAngleForce_destroy(destroy);
    destroy = 0;
}
int openmm_harmonicangleforce_getnumangles_(const OpenMM_HarmonicAngleForce*& target) {
    return OpenMM_HarmonicAngleForce_getNumAngles(target);
};
int OPENMM_HARMONICANGLEFORCE_GETNUMANGLES(const OpenMM_HarmonicAngleForce*& target) {
    return OpenMM_HarmonicAngleForce_getNumAngles(target);
};
int openmm_harmonicangleforce_addangle_(OpenMM_HarmonicAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, double const& angle, double const& k) {
    return OpenMM_HarmonicAngleForce_addAngle(target, particle1, particle2, particle3, angle, k);
};
int OPENMM_HARMONICANGLEFORCE_ADDANGLE(OpenMM_HarmonicAngleForce*& target, int const& particle1, int const& particle2, int const& particle3, double const& angle, double const& k) {
    return OpenMM_HarmonicAngleForce_addAngle(target, particle1, particle2, particle3, angle, k);
};
void openmm_harmonicangleforce_getangleparameters_(const OpenMM_HarmonicAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, double* angle, double* k) {
    OpenMM_HarmonicAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, angle, k);
};
void OPENMM_HARMONICANGLEFORCE_GETANGLEPARAMETERS(const OpenMM_HarmonicAngleForce*& target, int const& index, int* particle1, int* particle2, int* particle3, double* angle, double* k) {
    OpenMM_HarmonicAngleForce_getAngleParameters(target, index, particle1, particle2, particle3, angle, k);
};
void openmm_harmonicangleforce_setangleparameters_(OpenMM_HarmonicAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, double const& angle, double const& k) {
    OpenMM_HarmonicAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, angle, k);
};
void OPENMM_HARMONICANGLEFORCE_SETANGLEPARAMETERS(OpenMM_HarmonicAngleForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, double const& angle, double const& k) {
    OpenMM_HarmonicAngleForce_setAngleParameters(target, index, particle1, particle2, particle3, angle, k);
};


/* OpenMM::AndersenThermostat*/
void openmm_andersenthermostat_create_(OpenMM_AndersenThermostat*& result, double const& defaultTemperature, double const& defaultCollisionFrequency) {
    result = OpenMM_AndersenThermostat_create(defaultTemperature, defaultCollisionFrequency);
}
void OPENMM_ANDERSENTHERMOSTAT_CREATE(OpenMM_AndersenThermostat*& result, double const& defaultTemperature, double const& defaultCollisionFrequency) {
    result = OpenMM_AndersenThermostat_create(defaultTemperature, defaultCollisionFrequency);
}
void openmm_andersenthermostat_destroy_(OpenMM_AndersenThermostat*& destroy) {
    OpenMM_AndersenThermostat_destroy(destroy);
    destroy = 0;
}
void OPENMM_ANDERSENTHERMOSTAT_DESTROY(OpenMM_AndersenThermostat*& destroy) {
    OpenMM_AndersenThermostat_destroy(destroy);
    destroy = 0;
}
void openmm_andersenthermostat_temperature_(char* result, int result_length) {
    const char* result_chars = OpenMM_AndersenThermostat_Temperature();
    strncpy(result, result_chars, result_length);
 
};
void OPENMM_ANDERSENTHERMOSTAT_TEMPERATURE(char* result, int result_length) {
    const char* result_chars = OpenMM_AndersenThermostat_Temperature();
    strncpy(result, result_chars, result_length);
 
};
void openmm_andersenthermostat_collisionfrequency_(char* result, int result_length) {
    const char* result_chars = OpenMM_AndersenThermostat_CollisionFrequency();
    strncpy(result, result_chars, result_length);
 
};
void OPENMM_ANDERSENTHERMOSTAT_COLLISIONFREQUENCY(char* result, int result_length) {
    const char* result_chars = OpenMM_AndersenThermostat_CollisionFrequency();
    strncpy(result, result_chars, result_length);
 
};
double openmm_andersenthermostat_getdefaulttemperature_(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getDefaultTemperature(target);
};
double OPENMM_ANDERSENTHERMOSTAT_GETDEFAULTTEMPERATURE(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getDefaultTemperature(target);
};
double openmm_andersenthermostat_getdefaultcollisionfrequency_(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getDefaultCollisionFrequency(target);
};
double OPENMM_ANDERSENTHERMOSTAT_GETDEFAULTCOLLISIONFREQUENCY(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getDefaultCollisionFrequency(target);
};
int openmm_andersenthermostat_getrandomnumberseed_(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getRandomNumberSeed(target);
};
int OPENMM_ANDERSENTHERMOSTAT_GETRANDOMNUMBERSEED(const OpenMM_AndersenThermostat*& target) {
    return OpenMM_AndersenThermostat_getRandomNumberSeed(target);
};
void openmm_andersenthermostat_setrandomnumberseed_(OpenMM_AndersenThermostat*& target, int const& seed) {
    OpenMM_AndersenThermostat_setRandomNumberSeed(target, seed);
};
void OPENMM_ANDERSENTHERMOSTAT_SETRANDOMNUMBERSEED(OpenMM_AndersenThermostat*& target, int const& seed) {
    OpenMM_AndersenThermostat_setRandomNumberSeed(target, seed);
};


/* OpenMM::Platform*/
void openmm_platform_destroy_(OpenMM_Platform*& destroy) {
    OpenMM_Platform_destroy(destroy);
    destroy = 0;
}
void OPENMM_PLATFORM_DESTROY(OpenMM_Platform*& destroy) {
    OpenMM_Platform_destroy(destroy);
    destroy = 0;
}
void openmm_platform_getname_(const OpenMM_Platform*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getName(target);
    strncpy(result, result_chars, result_length);
 
};
void OPENMM_PLATFORM_GETNAME(const OpenMM_Platform*& target, char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getName(target);
    strncpy(result, result_chars, result_length);
 
};
double openmm_platform_getspeed_(const OpenMM_Platform*& target) {
    return OpenMM_Platform_getSpeed(target);
};
double OPENMM_PLATFORM_GETSPEED(const OpenMM_Platform*& target) {
    return OpenMM_Platform_getSpeed(target);
};
void openmm_platform_supportsdoubleprecision_(const OpenMM_Platform*& target, OpenMM_Boolean& result) {
    result = OpenMM_Platform_supportsDoublePrecision(target);
};
void OPENMM_PLATFORM_SUPPORTSDOUBLEPRECISION(const OpenMM_Platform*& target, OpenMM_Boolean& result) {
    result = OpenMM_Platform_supportsDoublePrecision(target);
};
void openmm_platform_getpropertynames_(OpenMM_Platform*& target, const OpenMM_StringArray*& result) {
    result = OpenMM_Platform_getPropertyNames(target);
};
void OPENMM_PLATFORM_GETPROPERTYNAMES(OpenMM_Platform*& target, const OpenMM_StringArray*& result) {
    result = OpenMM_Platform_getPropertyNames(target);
};
void openmm_platform_getpropertyvalue_(const OpenMM_Platform*& target, const OpenMM_Context*& context, const char* property, char* result, int property_length, int result_length) {
    const char* result_chars = OpenMM_Platform_getPropertyValue(target, context, string(property, property_length).c_str());
    strncpy(result, result_chars, result_length);
 
};
void OPENMM_PLATFORM_GETPROPERTYVALUE(const OpenMM_Platform*& target, const OpenMM_Context*& context, const char* property, char* result, int property_length, int result_length) {
    const char* result_chars = OpenMM_Platform_getPropertyValue(target, context, string(property, property_length).c_str());
    strncpy(result, result_chars, result_length);
 
};
void openmm_platform_setpropertyvalue_(const OpenMM_Platform*& target, OpenMM_Context*& context, const char* property, const char* value, int property_length, int value_length) {
    OpenMM_Platform_setPropertyValue(target, context, string(property, property_length).c_str(), string(value, value_length).c_str());
};
void OPENMM_PLATFORM_SETPROPERTYVALUE(const OpenMM_Platform*& target, OpenMM_Context*& context, const char* property, const char* value, int property_length, int value_length) {
    OpenMM_Platform_setPropertyValue(target, context, string(property, property_length).c_str(), string(value, value_length).c_str());
};
void openmm_platform_getpropertydefaultvalue_(const OpenMM_Platform*& target, const char* property, char* result, int property_length, int result_length) {
    const char* result_chars = OpenMM_Platform_getPropertyDefaultValue(target, string(property, property_length).c_str());
    strncpy(result, result_chars, result_length);
 
};
void OPENMM_PLATFORM_GETPROPERTYDEFAULTVALUE(const OpenMM_Platform*& target, const char* property, char* result, int property_length, int result_length) {
    const char* result_chars = OpenMM_Platform_getPropertyDefaultValue(target, string(property, property_length).c_str());
    strncpy(result, result_chars, result_length);
 
};
void openmm_platform_setpropertydefaultvalue_(OpenMM_Platform*& target, const char* property, const char* value, int property_length, int value_length) {
    OpenMM_Platform_setPropertyDefaultValue(target, string(property, property_length).c_str(), string(value, value_length).c_str());
};
void OPENMM_PLATFORM_SETPROPERTYDEFAULTVALUE(OpenMM_Platform*& target, const char* property, const char* value, int property_length, int value_length) {
    OpenMM_Platform_setPropertyDefaultValue(target, string(property, property_length).c_str(), string(value, value_length).c_str());
};
void openmm_platform_contextcreated_(const OpenMM_Platform*& target, OpenMM_ContextImpl* context) {
    OpenMM_Platform_contextCreated(target, context);
};
void OPENMM_PLATFORM_CONTEXTCREATED(const OpenMM_Platform*& target, OpenMM_ContextImpl* context) {
    OpenMM_Platform_contextCreated(target, context);
};
void openmm_platform_contextdestroyed_(const OpenMM_Platform*& target, OpenMM_ContextImpl* context) {
    OpenMM_Platform_contextDestroyed(target, context);
};
void OPENMM_PLATFORM_CONTEXTDESTROYED(const OpenMM_Platform*& target, OpenMM_ContextImpl* context) {
    OpenMM_Platform_contextDestroyed(target, context);
};
void openmm_platform_supportskernels_(const OpenMM_Platform*& target, const OpenMM_StringArray*& kernelNames, OpenMM_Boolean& result) {
    result = OpenMM_Platform_supportsKernels(target, kernelNames);
};
void OPENMM_PLATFORM_SUPPORTSKERNELS(const OpenMM_Platform*& target, const OpenMM_StringArray*& kernelNames, OpenMM_Boolean& result) {
    result = OpenMM_Platform_supportsKernels(target, kernelNames);
};
void openmm_platform_registerplatform_(OpenMM_Platform*& platform) {
    OpenMM_Platform_registerPlatform(platform);
};
void OPENMM_PLATFORM_REGISTERPLATFORM(OpenMM_Platform*& platform) {
    OpenMM_Platform_registerPlatform(platform);
};
int openmm_platform_getnumplatforms_() {
    return OpenMM_Platform_getNumPlatforms();
};
int OPENMM_PLATFORM_GETNUMPLATFORMS() {
    return OpenMM_Platform_getNumPlatforms();
};
void openmm_platform_getplatform_(int const& index, OpenMM_Platform*& result) {
    result = OpenMM_Platform_getPlatform(index);
};
void OPENMM_PLATFORM_GETPLATFORM(int const& index, OpenMM_Platform*& result) {
    result = OpenMM_Platform_getPlatform(index);
};
void openmm_platform_findplatform_(const OpenMM_StringArray*& kernelNames, OpenMM_Platform*& result) {
    result = OpenMM_Platform_findPlatform(kernelNames);
};
void OPENMM_PLATFORM_FINDPLATFORM(const OpenMM_StringArray*& kernelNames, OpenMM_Platform*& result) {
    result = OpenMM_Platform_findPlatform(kernelNames);
};
void openmm_platform_loadpluginlibrary_(const char* file, int file_length) {
    OpenMM_Platform_loadPluginLibrary(string(file, file_length).c_str());
};
void OPENMM_PLATFORM_LOADPLUGINLIBRARY(const char* file, int file_length) {
    OpenMM_Platform_loadPluginLibrary(string(file, file_length).c_str());
};
void openmm_platform_getdefaultpluginsdirectory_(char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getDefaultPluginsDirectory();
    strncpy(result, result_chars, result_length);
 
};
void OPENMM_PLATFORM_GETDEFAULTPLUGINSDIRECTORY(char* result, int result_length) {
    const char* result_chars = OpenMM_Platform_getDefaultPluginsDirectory();
    strncpy(result, result_chars, result_length);
 
};


/* OpenMM::State*/
void openmm_state_destroy_(OpenMM_State*& destroy) {
    OpenMM_State_destroy(destroy);
    destroy = 0;
}
void OPENMM_STATE_DESTROY(OpenMM_State*& destroy) {
    OpenMM_State_destroy(destroy);
    destroy = 0;
}
double openmm_state_gettime_(const OpenMM_State*& target) {
    return OpenMM_State_getTime(target);
};
double OPENMM_STATE_GETTIME(const OpenMM_State*& target) {
    return OpenMM_State_getTime(target);
};
void openmm_state_getpositions_(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getPositions(target);
};
void OPENMM_STATE_GETPOSITIONS(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getPositions(target);
};
void openmm_state_getvelocities_(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getVelocities(target);
};
void OPENMM_STATE_GETVELOCITIES(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getVelocities(target);
};
void openmm_state_getforces_(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getForces(target);
};
void OPENMM_STATE_GETFORCES(const OpenMM_State*& target, const OpenMM_Vec3Array*& result) {
    result = OpenMM_State_getForces(target);
};
double openmm_state_getkineticenergy_(const OpenMM_State*& target) {
    return OpenMM_State_getKineticEnergy(target);
};
double OPENMM_STATE_GETKINETICENERGY(const OpenMM_State*& target) {
    return OpenMM_State_getKineticEnergy(target);
};
double openmm_state_getpotentialenergy_(const OpenMM_State*& target) {
    return OpenMM_State_getPotentialEnergy(target);
};
double OPENMM_STATE_GETPOTENTIALENERGY(const OpenMM_State*& target) {
    return OpenMM_State_getPotentialEnergy(target);
};
void openmm_state_getparameters_(const OpenMM_State*& target, const OpenMM_ParameterArray*& result) {
    result = OpenMM_State_getParameters(target);
};
void OPENMM_STATE_GETPARAMETERS(const OpenMM_State*& target, const OpenMM_ParameterArray*& result) {
    result = OpenMM_State_getParameters(target);
};


/* OpenMM::PeriodicTorsionForce*/
void openmm_periodictorsionforce_create_(OpenMM_PeriodicTorsionForce*& result) {
    result = OpenMM_PeriodicTorsionForce_create();
}
void OPENMM_PERIODICTORSIONFORCE_CREATE(OpenMM_PeriodicTorsionForce*& result) {
    result = OpenMM_PeriodicTorsionForce_create();
}
void openmm_periodictorsionforce_destroy_(OpenMM_PeriodicTorsionForce*& destroy) {
    OpenMM_PeriodicTorsionForce_destroy(destroy);
    destroy = 0;
}
void OPENMM_PERIODICTORSIONFORCE_DESTROY(OpenMM_PeriodicTorsionForce*& destroy) {
    OpenMM_PeriodicTorsionForce_destroy(destroy);
    destroy = 0;
}
int openmm_periodictorsionforce_getnumtorsions_(const OpenMM_PeriodicTorsionForce*& target) {
    return OpenMM_PeriodicTorsionForce_getNumTorsions(target);
};
int OPENMM_PERIODICTORSIONFORCE_GETNUMTORSIONS(const OpenMM_PeriodicTorsionForce*& target) {
    return OpenMM_PeriodicTorsionForce_getNumTorsions(target);
};
int openmm_periodictorsionforce_addtorsion_(OpenMM_PeriodicTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& periodicity, double const& phase, double const& k) {
    return OpenMM_PeriodicTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, periodicity, phase, k);
};
int OPENMM_PERIODICTORSIONFORCE_ADDTORSION(OpenMM_PeriodicTorsionForce*& target, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& periodicity, double const& phase, double const& k) {
    return OpenMM_PeriodicTorsionForce_addTorsion(target, particle1, particle2, particle3, particle4, periodicity, phase, k);
};
void openmm_periodictorsionforce_gettorsionparameters_(const OpenMM_PeriodicTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, int* periodicity, double* phase, double* k) {
    OpenMM_PeriodicTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k);
};
void OPENMM_PERIODICTORSIONFORCE_GETTORSIONPARAMETERS(const OpenMM_PeriodicTorsionForce*& target, int const& index, int* particle1, int* particle2, int* particle3, int* particle4, int* periodicity, double* phase, double* k) {
    OpenMM_PeriodicTorsionForce_getTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k);
};
void openmm_periodictorsionforce_settorsionparameters_(OpenMM_PeriodicTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& periodicity, double const& phase, double const& k) {
    OpenMM_PeriodicTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k);
};
void OPENMM_PERIODICTORSIONFORCE_SETTORSIONPARAMETERS(OpenMM_PeriodicTorsionForce*& target, int const& index, int const& particle1, int const& particle2, int const& particle3, int const& particle4, int const& periodicity, double const& phase, double const& k) {
    OpenMM_PeriodicTorsionForce_setTorsionParameters(target, index, particle1, particle2, particle3, particle4, periodicity, phase, k);
};


/* OpenMM::Integrator*/
void openmm_integrator_destroy_(OpenMM_Integrator*& destroy) {
    OpenMM_Integrator_destroy(destroy);
    destroy = 0;
}
void OPENMM_INTEGRATOR_DESTROY(OpenMM_Integrator*& destroy) {
    OpenMM_Integrator_destroy(destroy);
    destroy = 0;
}
double openmm_integrator_getstepsize_(const OpenMM_Integrator*& target) {
    return OpenMM_Integrator_getStepSize(target);
};
double OPENMM_INTEGRATOR_GETSTEPSIZE(const OpenMM_Integrator*& target) {
    return OpenMM_Integrator_getStepSize(target);
};
void openmm_integrator_setstepsize_(OpenMM_Integrator*& target, double const& size) {
    OpenMM_Integrator_setStepSize(target, size);
};
void OPENMM_INTEGRATOR_SETSTEPSIZE(OpenMM_Integrator*& target, double const& size) {
    OpenMM_Integrator_setStepSize(target, size);
};
double openmm_integrator_getconstrainttolerance_(const OpenMM_Integrator*& target) {
    return OpenMM_Integrator_getConstraintTolerance(target);
};
double OPENMM_INTEGRATOR_GETCONSTRAINTTOLERANCE(const OpenMM_Integrator*& target) {
    return OpenMM_Integrator_getConstraintTolerance(target);
};
void openmm_integrator_setconstrainttolerance_(OpenMM_Integrator*& target, double const& tol) {
    OpenMM_Integrator_setConstraintTolerance(target, tol);
};
void OPENMM_INTEGRATOR_SETCONSTRAINTTOLERANCE(OpenMM_Integrator*& target, double const& tol) {
    OpenMM_Integrator_setConstraintTolerance(target, tol);
};
void openmm_integrator_step_(OpenMM_Integrator*& target, int const& steps) {
    OpenMM_Integrator_step(target, steps);
};
void OPENMM_INTEGRATOR_STEP(OpenMM_Integrator*& target, int const& steps) {
    OpenMM_Integrator_step(target, steps);
};


/* OpenMM::System*/
void openmm_system_create_(OpenMM_System*& result) {
    result = OpenMM_System_create();
}
void OPENMM_SYSTEM_CREATE(OpenMM_System*& result) {
    result = OpenMM_System_create();
}
void openmm_system_destroy_(OpenMM_System*& destroy) {
    OpenMM_System_destroy(destroy);
    destroy = 0;
}
void OPENMM_SYSTEM_DESTROY(OpenMM_System*& destroy) {
    OpenMM_System_destroy(destroy);
    destroy = 0;
}
int openmm_system_getnumparticles_(const OpenMM_System*& target) {
    return OpenMM_System_getNumParticles(target);
};
int OPENMM_SYSTEM_GETNUMPARTICLES(const OpenMM_System*& target) {
    return OpenMM_System_getNumParticles(target);
};
int openmm_system_addparticle_(OpenMM_System*& target, double const& mass) {
    return OpenMM_System_addParticle(target, mass);
};
int OPENMM_SYSTEM_ADDPARTICLE(OpenMM_System*& target, double const& mass) {
    return OpenMM_System_addParticle(target, mass);
};
double openmm_system_getparticlemass_(const OpenMM_System*& target, int const& index) {
    return OpenMM_System_getParticleMass(target, index);
};
double OPENMM_SYSTEM_GETPARTICLEMASS(const OpenMM_System*& target, int const& index) {
    return OpenMM_System_getParticleMass(target, index);
};
void openmm_system_setparticlemass_(OpenMM_System*& target, int const& index, double const& mass) {
    OpenMM_System_setParticleMass(target, index, mass);
};
void OPENMM_SYSTEM_SETPARTICLEMASS(OpenMM_System*& target, int const& index, double const& mass) {
    OpenMM_System_setParticleMass(target, index, mass);
};
int openmm_system_getnumconstraints_(const OpenMM_System*& target) {
    return OpenMM_System_getNumConstraints(target);
};
int OPENMM_SYSTEM_GETNUMCONSTRAINTS(const OpenMM_System*& target) {
    return OpenMM_System_getNumConstraints(target);
};
int openmm_system_addconstraint_(OpenMM_System*& target, int const& particle1, int const& particle2, double const& distance) {
    return OpenMM_System_addConstraint(target, particle1, particle2, distance);
};
int OPENMM_SYSTEM_ADDCONSTRAINT(OpenMM_System*& target, int const& particle1, int const& particle2, double const& distance) {
    return OpenMM_System_addConstraint(target, particle1, particle2, distance);
};
void openmm_system_getconstraintparameters_(const OpenMM_System*& target, int const& index, int* particle1, int* particle2, double* distance) {
    OpenMM_System_getConstraintParameters(target, index, particle1, particle2, distance);
};
void OPENMM_SYSTEM_GETCONSTRAINTPARAMETERS(const OpenMM_System*& target, int const& index, int* particle1, int* particle2, double* distance) {
    OpenMM_System_getConstraintParameters(target, index, particle1, particle2, distance);
};
void openmm_system_setconstraintparameters_(OpenMM_System*& target, int const& index, int const& particle1, int const& particle2, double const& distance) {
    OpenMM_System_setConstraintParameters(target, index, particle1, particle2, distance);
};
void OPENMM_SYSTEM_SETCONSTRAINTPARAMETERS(OpenMM_System*& target, int const& index, int const& particle1, int const& particle2, double const& distance) {
    OpenMM_System_setConstraintParameters(target, index, particle1, particle2, distance);
};
int openmm_system_addforce_(OpenMM_System*& target, OpenMM_Force*& force) {
    return OpenMM_System_addForce(target, force);
};
int OPENMM_SYSTEM_ADDFORCE(OpenMM_System*& target, OpenMM_Force*& force) {
    return OpenMM_System_addForce(target, force);
};
int openmm_system_getnumforces_(const OpenMM_System*& target) {
    return OpenMM_System_getNumForces(target);
};
int OPENMM_SYSTEM_GETNUMFORCES(const OpenMM_System*& target) {
    return OpenMM_System_getNumForces(target);
};
void openmm_system_getforce_(OpenMM_System*& target, int const& index, OpenMM_Force*& result) {
    result = OpenMM_System_getForce(target, index);
};
void OPENMM_SYSTEM_GETFORCE(OpenMM_System*& target, int const& index, OpenMM_Force*& result) {
    result = OpenMM_System_getForce(target, index);
};


#if defined(__cplusplus)
}
#endif
