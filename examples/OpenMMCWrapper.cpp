
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
OpenMM_Vec3 OpenMM_Vec3_scale(const OpenMM_Vec3 vec, double scale) {
    OpenMM_Vec3 result = {vec.x*scale, vec.y*scale, vec.z*scale};
    return result;
}

/* OpenMM_Vec3Array */
OpenMM_Vec3Array* OpenMM_Vec3Array_create(int size) {
    return reinterpret_cast<OpenMM_Vec3Array*>(new vector<Vec3>(size));
}
void OpenMM_Vec3Array_destroy(OpenMM_Vec3Array* array) {
    delete reinterpret_cast<vector<Vec3>*>(array);
}
int OpenMM_Vec3Array_getSize(const OpenMM_Vec3Array* array) {
    return reinterpret_cast<const vector<Vec3>*>(array)->size();
}
void OpenMM_Vec3Array_resize(OpenMM_Vec3Array* array, int size) {
    reinterpret_cast<vector<Vec3>*>(array)->resize(size);
}
void OpenMM_Vec3Array_append(OpenMM_Vec3Array* array, const OpenMM_Vec3 vec) {
    reinterpret_cast<vector<Vec3>*>(array)->push_back(Vec3(vec.x, vec.y, vec.z));
}
void OpenMM_Vec3Array_set(OpenMM_Vec3Array* array, int index, const OpenMM_Vec3 vec) {
    (*reinterpret_cast<vector<Vec3>*>(array))[index] = Vec3(vec.x, vec.y, vec.z);
}
const OpenMM_Vec3* OpenMM_Vec3Array_get(const OpenMM_Vec3Array* array, int index) {
    return reinterpret_cast<const OpenMM_Vec3*>((&(*reinterpret_cast<const vector<Vec3>*>(array))[index]));
}

/* OpenMM_StringArray */
OpenMM_StringArray* OpenMM_StringArray_create(int size) {
    return reinterpret_cast<OpenMM_StringArray*>(new vector<string>(size));
}
void OpenMM_StringArray_destroy(OpenMM_StringArray* array) {
    delete reinterpret_cast<vector<string>*>(array);
}
int OpenMM_StringArray_getSize(const OpenMM_StringArray* array) {
    return reinterpret_cast<const vector<string>*>(array)->size();
}
void OpenMM_StringArray_resize(OpenMM_StringArray* array, int size) {
    reinterpret_cast<vector<string>*>(array)->resize(size);
}
void OpenMM_StringArray_append(OpenMM_StringArray* array, const char* str) {
    reinterpret_cast<vector<string>*>(array)->push_back(string(str));
}
void OpenMM_StringArray_set(OpenMM_StringArray* array, int index, const char* str) {
    (*reinterpret_cast<vector<string>*>(array))[index] = string(str);
}
const char* OpenMM_StringArray_get(const OpenMM_StringArray* array, int index) {
    return (*reinterpret_cast<const vector<string>*>(array))[index].c_str();
}

/* OpenMM_BondArray */
OpenMM_BondArray* OpenMM_BondArray_create(int size) {
    return reinterpret_cast<OpenMM_BondArray*>(new vector<pair<int, int> >(size));
}
void OpenMM_BondArray_destroy(OpenMM_BondArray* array) {
    delete reinterpret_cast<vector<pair<int, int> >*>(array);
}
int OpenMM_BondArray_getSize(const OpenMM_BondArray* array) {
    return reinterpret_cast<const vector<pair<int, int> >*>(array)->size();
}
void OpenMM_BondArray_resize(OpenMM_BondArray* array, int size) {
    reinterpret_cast<vector<pair<int, int> >*>(array)->resize(size);
}
void OpenMM_BondArray_append(OpenMM_BondArray* array, int particle1, int particle2) {
    reinterpret_cast<vector<pair<int, int> >*>(array)->push_back(pair<int, int>(particle1, particle2));
}
void OpenMM_BondArray_set(OpenMM_BondArray* array, int index, int particle1, int particle2) {
    (*reinterpret_cast<vector<pair<int, int> >*>(array))[index] = pair<int, int>(particle1, particle2);
}
void OpenMM_BondArray_get(const OpenMM_BondArray* array, int index, int* particle1, int* particle2) {
    pair<int, int> particles = (*reinterpret_cast<const vector<pair<int, int> >*>(array))[index];
    *particle1 = particles.first;
    *particle2 = particles.second;
}

/* OpenMM_ParameterArray */
int OpenMM_ParameterArray_getSize(const OpenMM_ParameterArray* array) {
    return reinterpret_cast<const map<string, double>*>(array)->size();
}
double OpenMM_ParameterArray_get(const OpenMM_ParameterArray* array, const char* name) {
    const map<string, double>* params = reinterpret_cast<const map<string, double>*>(array);
    const map<string, double>::const_iterator iter = params->find(string(name));
    if (iter == params->end())
        throw OpenMMException("OpenMM_ParameterArray_get: No such parameter");
    return iter->second;
}

/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
OpenMM_State* OpenMM_Context_getState(const OpenMM_Context* target, int types) {
    State result = reinterpret_cast<const Context*>(target)->getState(types);
    return reinterpret_cast<OpenMM_State*>(new State(result));
};
OpenMM_StringArray* OpenMM_Platform_loadPluginsFromDirectory(const char* directory) {
    vector<string> result = Platform::loadPluginsFromDirectory(string(directory));
    return reinterpret_cast<OpenMM_StringArray*>(new vector<string>(result));
};

 
 

/* OpenMM::HarmonicBondForce*/
OpenMM_HarmonicBondForce* OpenMM_HarmonicBondForce_create() {
    return reinterpret_cast<OpenMM_HarmonicBondForce*>(new HarmonicBondForce());
}

void OpenMM_HarmonicBondForce_destroy(OpenMM_HarmonicBondForce* target) {
    delete reinterpret_cast<HarmonicBondForce*>(target);
}
 
 
int OpenMM_HarmonicBondForce_getNumBonds(const OpenMM_HarmonicBondForce* target) {
    int result = reinterpret_cast<const HarmonicBondForce*>(target)->getNumBonds();
    return result;
};

int OpenMM_HarmonicBondForce_addBond(OpenMM_HarmonicBondForce* target, int particle1, int particle2, double length, double k) {
    int result = reinterpret_cast<HarmonicBondForce*>(target)->addBond(particle1, particle2, length, k);
    return result;
};

void OpenMM_HarmonicBondForce_getBondParameters(const OpenMM_HarmonicBondForce* target, int index, int* particle1, int* particle2, double* length, double* k) {
    reinterpret_cast<const HarmonicBondForce*>(target)->getBondParameters(index, *reinterpret_cast<int* >(particle1), *reinterpret_cast<int* >(particle2), *reinterpret_cast<double* >(length), *reinterpret_cast<double* >(k));
};

void OpenMM_HarmonicBondForce_setBondParameters(OpenMM_HarmonicBondForce* target, int index, int particle1, int particle2, double length, double k) {
    reinterpret_cast<HarmonicBondForce*>(target)->setBondParameters(index, particle1, particle2, length, k);
};


/* OpenMM::BrownianIntegrator*/
OpenMM_BrownianIntegrator* OpenMM_BrownianIntegrator_create(double temperature, double frictionCoeff, double stepSize) {
    return reinterpret_cast<OpenMM_BrownianIntegrator*>(new BrownianIntegrator(temperature, frictionCoeff, stepSize));
}

void OpenMM_BrownianIntegrator_destroy(OpenMM_BrownianIntegrator* target) {
    delete reinterpret_cast<BrownianIntegrator*>(target);
}
 
 
double OpenMM_BrownianIntegrator_getTemperature(const OpenMM_BrownianIntegrator* target) {
    double result = reinterpret_cast<const BrownianIntegrator*>(target)->getTemperature();
    return result;
};

void OpenMM_BrownianIntegrator_setTemperature(OpenMM_BrownianIntegrator* target, double temp) {
    reinterpret_cast<BrownianIntegrator*>(target)->setTemperature(temp);
};

double OpenMM_BrownianIntegrator_getFriction(const OpenMM_BrownianIntegrator* target) {
    double result = reinterpret_cast<const BrownianIntegrator*>(target)->getFriction();
    return result;
};

void OpenMM_BrownianIntegrator_setFriction(OpenMM_BrownianIntegrator* target, double coeff) {
    reinterpret_cast<BrownianIntegrator*>(target)->setFriction(coeff);
};

int OpenMM_BrownianIntegrator_getRandomNumberSeed(const OpenMM_BrownianIntegrator* target) {
    int result = reinterpret_cast<const BrownianIntegrator*>(target)->getRandomNumberSeed();
    return result;
};

void OpenMM_BrownianIntegrator_setRandomNumberSeed(OpenMM_BrownianIntegrator* target, int seed) {
    reinterpret_cast<BrownianIntegrator*>(target)->setRandomNumberSeed(seed);
};

void OpenMM_BrownianIntegrator_step(OpenMM_BrownianIntegrator* target, int steps) {
    reinterpret_cast<BrownianIntegrator*>(target)->step(steps);
};


/* OpenMM::OpenMMException*/
OpenMM_OpenMMException* OpenMM_OpenMMException_create(const char* message) {
    return reinterpret_cast<OpenMM_OpenMMException*>(new OpenMMException(string(message)));
}

void OpenMM_OpenMMException_destroy(OpenMM_OpenMMException* target) {
    delete reinterpret_cast<OpenMMException*>(target);
}
 
 
const char* OpenMM_OpenMMException_what(const OpenMM_OpenMMException* target) {
    const char* result = reinterpret_cast<const OpenMMException*>(target)->what();
    return reinterpret_cast<const char*>(result);
};


/* OpenMM::NonbondedForce*/
OpenMM_NonbondedForce* OpenMM_NonbondedForce_create() {
    return reinterpret_cast<OpenMM_NonbondedForce*>(new NonbondedForce());
}

void OpenMM_NonbondedForce_destroy(OpenMM_NonbondedForce* target) {
    delete reinterpret_cast<NonbondedForce*>(target);
}
 
 
int OpenMM_NonbondedForce_getNumParticles(const OpenMM_NonbondedForce* target) {
    int result = reinterpret_cast<const NonbondedForce*>(target)->getNumParticles();
    return result;
};

int OpenMM_NonbondedForce_getNumExceptions(const OpenMM_NonbondedForce* target) {
    int result = reinterpret_cast<const NonbondedForce*>(target)->getNumExceptions();
    return result;
};

OpenMM_NonbondedForce_NonbondedMethod OpenMM_NonbondedForce_getNonbondedMethod(const OpenMM_NonbondedForce* target) {
    NonbondedForce::NonbondedMethod result = reinterpret_cast<const NonbondedForce*>(target)->getNonbondedMethod();
    return static_cast<OpenMM_NonbondedForce_NonbondedMethod>(result);
};

void OpenMM_NonbondedForce_setNonbondedMethod(OpenMM_NonbondedForce* target, OpenMM_NonbondedForce_NonbondedMethod method) {
    reinterpret_cast<NonbondedForce*>(target)->setNonbondedMethod(static_cast<NonbondedForce::NonbondedMethod >(method));
};

double OpenMM_NonbondedForce_getCutoffDistance(const OpenMM_NonbondedForce* target) {
    double result = reinterpret_cast<const NonbondedForce*>(target)->getCutoffDistance();
    return result;
};

void OpenMM_NonbondedForce_setCutoffDistance(OpenMM_NonbondedForce* target, double distance) {
    reinterpret_cast<NonbondedForce*>(target)->setCutoffDistance(distance);
};

double OpenMM_NonbondedForce_getReactionFieldDielectric(const OpenMM_NonbondedForce* target) {
    double result = reinterpret_cast<const NonbondedForce*>(target)->getReactionFieldDielectric();
    return result;
};

void OpenMM_NonbondedForce_setReactionFieldDielectric(OpenMM_NonbondedForce* target, double dielectric) {
    reinterpret_cast<NonbondedForce*>(target)->setReactionFieldDielectric(dielectric);
};

double OpenMM_NonbondedForce_getEwaldErrorTolerance(const OpenMM_NonbondedForce* target) {
    double result = reinterpret_cast<const NonbondedForce*>(target)->getEwaldErrorTolerance();
    return result;
};

void OpenMM_NonbondedForce_setEwaldErrorTolerance(OpenMM_NonbondedForce* target, double tol) {
    reinterpret_cast<NonbondedForce*>(target)->setEwaldErrorTolerance(tol);
};

void OpenMM_NonbondedForce_getPeriodicBoxVectors(const OpenMM_NonbondedForce* target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c) {
    reinterpret_cast<const NonbondedForce*>(target)->getPeriodicBoxVectors(*reinterpret_cast<Vec3* >(a), *reinterpret_cast<Vec3* >(b), *reinterpret_cast<Vec3* >(c));
};

void OpenMM_NonbondedForce_setPeriodicBoxVectors(OpenMM_NonbondedForce* target, OpenMM_Vec3 a, OpenMM_Vec3 b, OpenMM_Vec3 c) {
    reinterpret_cast<NonbondedForce*>(target)->setPeriodicBoxVectors(Vec3(a.x, a.y, a.z), Vec3(b.x, b.y, b.z), Vec3(c.x, c.y, c.z));
};

int OpenMM_NonbondedForce_addParticle(OpenMM_NonbondedForce* target, double charge, double sigma, double epsilon) {
    int result = reinterpret_cast<NonbondedForce*>(target)->addParticle(charge, sigma, epsilon);
    return result;
};

void OpenMM_NonbondedForce_getParticleParameters(const OpenMM_NonbondedForce* target, int index, double* charge, double* sigma, double* epsilon) {
    reinterpret_cast<const NonbondedForce*>(target)->getParticleParameters(index, *reinterpret_cast<double* >(charge), *reinterpret_cast<double* >(sigma), *reinterpret_cast<double* >(epsilon));
};

void OpenMM_NonbondedForce_setParticleParameters(OpenMM_NonbondedForce* target, int index, double charge, double sigma, double epsilon) {
    reinterpret_cast<NonbondedForce*>(target)->setParticleParameters(index, charge, sigma, epsilon);
};

int OpenMM_NonbondedForce_addException(OpenMM_NonbondedForce* target, int particle1, int particle2, double chargeProd, double sigma, double epsilon, OpenMM_Boolean replace) {
    int result = reinterpret_cast<NonbondedForce*>(target)->addException(particle1, particle2, chargeProd, sigma, epsilon, (replace != OpenMM_False));
    return result;
};

void OpenMM_NonbondedForce_getExceptionParameters(const OpenMM_NonbondedForce* target, int index, int* particle1, int* particle2, double* chargeProd, double* sigma, double* epsilon) {
    reinterpret_cast<const NonbondedForce*>(target)->getExceptionParameters(index, *reinterpret_cast<int* >(particle1), *reinterpret_cast<int* >(particle2), *reinterpret_cast<double* >(chargeProd), *reinterpret_cast<double* >(sigma), *reinterpret_cast<double* >(epsilon));
};

void OpenMM_NonbondedForce_setExceptionParameters(OpenMM_NonbondedForce* target, int index, int particle1, int particle2, double chargeProd, double sigma, double epsilon) {
    reinterpret_cast<NonbondedForce*>(target)->setExceptionParameters(index, particle1, particle2, chargeProd, sigma, epsilon);
};

void OpenMM_NonbondedForce_createExceptionsFromBonds(OpenMM_NonbondedForce* target, const OpenMM_BondArray* bonds, double coulomb14Scale, double lj14Scale) {
    reinterpret_cast<NonbondedForce*>(target)->createExceptionsFromBonds(*reinterpret_cast<const vector<pair<int, int> >* >(bonds), coulomb14Scale, lj14Scale);
};


/* OpenMM::VariableLangevinIntegrator*/
OpenMM_VariableLangevinIntegrator* OpenMM_VariableLangevinIntegrator_create(double temperature, double frictionCoeff, double errorTol) {
    return reinterpret_cast<OpenMM_VariableLangevinIntegrator*>(new VariableLangevinIntegrator(temperature, frictionCoeff, errorTol));
}

void OpenMM_VariableLangevinIntegrator_destroy(OpenMM_VariableLangevinIntegrator* target) {
    delete reinterpret_cast<VariableLangevinIntegrator*>(target);
}
 
 
double OpenMM_VariableLangevinIntegrator_getTemperature(const OpenMM_VariableLangevinIntegrator* target) {
    double result = reinterpret_cast<const VariableLangevinIntegrator*>(target)->getTemperature();
    return result;
};

void OpenMM_VariableLangevinIntegrator_setTemperature(OpenMM_VariableLangevinIntegrator* target, double temp) {
    reinterpret_cast<VariableLangevinIntegrator*>(target)->setTemperature(temp);
};

double OpenMM_VariableLangevinIntegrator_getFriction(const OpenMM_VariableLangevinIntegrator* target) {
    double result = reinterpret_cast<const VariableLangevinIntegrator*>(target)->getFriction();
    return result;
};

void OpenMM_VariableLangevinIntegrator_setFriction(OpenMM_VariableLangevinIntegrator* target, double coeff) {
    reinterpret_cast<VariableLangevinIntegrator*>(target)->setFriction(coeff);
};

double OpenMM_VariableLangevinIntegrator_getErrorTolerance(const OpenMM_VariableLangevinIntegrator* target) {
    double result = reinterpret_cast<const VariableLangevinIntegrator*>(target)->getErrorTolerance();
    return result;
};

void OpenMM_VariableLangevinIntegrator_setErrorTolerance(OpenMM_VariableLangevinIntegrator* target, double tol) {
    reinterpret_cast<VariableLangevinIntegrator*>(target)->setErrorTolerance(tol);
};

int OpenMM_VariableLangevinIntegrator_getRandomNumberSeed(const OpenMM_VariableLangevinIntegrator* target) {
    int result = reinterpret_cast<const VariableLangevinIntegrator*>(target)->getRandomNumberSeed();
    return result;
};

void OpenMM_VariableLangevinIntegrator_setRandomNumberSeed(OpenMM_VariableLangevinIntegrator* target, int seed) {
    reinterpret_cast<VariableLangevinIntegrator*>(target)->setRandomNumberSeed(seed);
};

void OpenMM_VariableLangevinIntegrator_step(OpenMM_VariableLangevinIntegrator* target, int steps) {
    reinterpret_cast<VariableLangevinIntegrator*>(target)->step(steps);
};

void OpenMM_VariableLangevinIntegrator_stepTo(OpenMM_VariableLangevinIntegrator* target, double time) {
    reinterpret_cast<VariableLangevinIntegrator*>(target)->stepTo(time);
};


/* OpenMM::GBVIForce*/
OpenMM_GBVIForce* OpenMM_GBVIForce_create() {
    return reinterpret_cast<OpenMM_GBVIForce*>(new GBVIForce());
}

void OpenMM_GBVIForce_destroy(OpenMM_GBVIForce* target) {
    delete reinterpret_cast<GBVIForce*>(target);
}
 
 
int OpenMM_GBVIForce_getNumParticles(const OpenMM_GBVIForce* target) {
    int result = reinterpret_cast<const GBVIForce*>(target)->getNumParticles();
    return result;
};

int OpenMM_GBVIForce_addParticle(OpenMM_GBVIForce* target, double charge, double radius, double gamma) {
    int result = reinterpret_cast<GBVIForce*>(target)->addParticle(charge, radius, gamma);
    return result;
};

void OpenMM_GBVIForce_getParticleParameters(const OpenMM_GBVIForce* target, int index, double* charge, double* radius, double* gamma) {
    reinterpret_cast<const GBVIForce*>(target)->getParticleParameters(index, *reinterpret_cast<double* >(charge), *reinterpret_cast<double* >(radius), *reinterpret_cast<double* >(gamma));
};

void OpenMM_GBVIForce_setParticleParameters(OpenMM_GBVIForce* target, int index, double charge, double radius, double gamma) {
    reinterpret_cast<GBVIForce*>(target)->setParticleParameters(index, charge, radius, gamma);
};

double OpenMM_GBVIForce_getSolventDielectric(const OpenMM_GBVIForce* target) {
    double result = reinterpret_cast<const GBVIForce*>(target)->getSolventDielectric();
    return result;
};

void OpenMM_GBVIForce_setSolventDielectric(OpenMM_GBVIForce* target, double dielectric) {
    reinterpret_cast<GBVIForce*>(target)->setSolventDielectric(dielectric);
};

double OpenMM_GBVIForce_getSoluteDielectric(const OpenMM_GBVIForce* target) {
    double result = reinterpret_cast<const GBVIForce*>(target)->getSoluteDielectric();
    return result;
};

void OpenMM_GBVIForce_setSoluteDielectric(OpenMM_GBVIForce* target, double dielectric) {
    reinterpret_cast<GBVIForce*>(target)->setSoluteDielectric(dielectric);
};


/* OpenMM::Context*/
OpenMM_Context* OpenMM_Context_create(OpenMM_System* system, OpenMM_Integrator* integrator) {
    return reinterpret_cast<OpenMM_Context*>(new Context(*reinterpret_cast<System* >(system), *reinterpret_cast<Integrator* >(integrator)));
}

OpenMM_Context* OpenMM_Context_create_2(OpenMM_System* system, OpenMM_Integrator* integrator, OpenMM_Platform* platform) {
    return reinterpret_cast<OpenMM_Context*>(new Context(*reinterpret_cast<System* >(system), *reinterpret_cast<Integrator* >(integrator), *reinterpret_cast<Platform* >(platform)));
}

void OpenMM_Context_destroy(OpenMM_Context* target) {
    delete reinterpret_cast<Context*>(target);
}
 
 
OpenMM_System* OpenMM_Context_getSystem(OpenMM_Context* target) {
    System* result = &reinterpret_cast<Context*>(target)->getSystem();
    return reinterpret_cast<OpenMM_System*>(result);
};

OpenMM_Integrator* OpenMM_Context_getIntegrator(OpenMM_Context* target) {
    Integrator* result = &reinterpret_cast<Context*>(target)->getIntegrator();
    return reinterpret_cast<OpenMM_Integrator*>(result);
};

OpenMM_Platform* OpenMM_Context_getPlatform(OpenMM_Context* target) {
    Platform* result = &reinterpret_cast<Context*>(target)->getPlatform();
    return reinterpret_cast<OpenMM_Platform*>(result);
};

void OpenMM_Context_setTime(OpenMM_Context* target, double time) {
    reinterpret_cast<Context*>(target)->setTime(time);
};

void OpenMM_Context_setPositions(OpenMM_Context* target, const OpenMM_Vec3Array* positions) {
    reinterpret_cast<Context*>(target)->setPositions(*reinterpret_cast<const vector<Vec3>* >(positions));
};

void OpenMM_Context_setVelocities(OpenMM_Context* target, const OpenMM_Vec3Array* velocities) {
    reinterpret_cast<Context*>(target)->setVelocities(*reinterpret_cast<const vector<Vec3>* >(velocities));
};

double OpenMM_Context_getParameter(OpenMM_Context* target, const char* name) {
    double result = reinterpret_cast<Context*>(target)->getParameter(string(name));
    return result;
};

void OpenMM_Context_setParameter(OpenMM_Context* target, const char* name, double value) {
    reinterpret_cast<Context*>(target)->setParameter(string(name), value);
};

void OpenMM_Context_reinitialize(OpenMM_Context* target) {
    reinterpret_cast<Context*>(target)->reinitialize();
};


/* OpenMM::GBSAOBCForce*/
OpenMM_GBSAOBCForce* OpenMM_GBSAOBCForce_create() {
    return reinterpret_cast<OpenMM_GBSAOBCForce*>(new GBSAOBCForce());
}

void OpenMM_GBSAOBCForce_destroy(OpenMM_GBSAOBCForce* target) {
    delete reinterpret_cast<GBSAOBCForce*>(target);
}
 
 
int OpenMM_GBSAOBCForce_getNumParticles(const OpenMM_GBSAOBCForce* target) {
    int result = reinterpret_cast<const GBSAOBCForce*>(target)->getNumParticles();
    return result;
};

int OpenMM_GBSAOBCForce_addParticle(OpenMM_GBSAOBCForce* target, double charge, double radius, double scalingFactor) {
    int result = reinterpret_cast<GBSAOBCForce*>(target)->addParticle(charge, radius, scalingFactor);
    return result;
};

void OpenMM_GBSAOBCForce_getParticleParameters(const OpenMM_GBSAOBCForce* target, int index, double* charge, double* radius, double* scalingFactor) {
    reinterpret_cast<const GBSAOBCForce*>(target)->getParticleParameters(index, *reinterpret_cast<double* >(charge), *reinterpret_cast<double* >(radius), *reinterpret_cast<double* >(scalingFactor));
};

void OpenMM_GBSAOBCForce_setParticleParameters(OpenMM_GBSAOBCForce* target, int index, double charge, double radius, double scalingFactor) {
    reinterpret_cast<GBSAOBCForce*>(target)->setParticleParameters(index, charge, radius, scalingFactor);
};

double OpenMM_GBSAOBCForce_getSolventDielectric(const OpenMM_GBSAOBCForce* target) {
    double result = reinterpret_cast<const GBSAOBCForce*>(target)->getSolventDielectric();
    return result;
};

void OpenMM_GBSAOBCForce_setSolventDielectric(OpenMM_GBSAOBCForce* target, double dielectric) {
    reinterpret_cast<GBSAOBCForce*>(target)->setSolventDielectric(dielectric);
};

double OpenMM_GBSAOBCForce_getSoluteDielectric(const OpenMM_GBSAOBCForce* target) {
    double result = reinterpret_cast<const GBSAOBCForce*>(target)->getSoluteDielectric();
    return result;
};

void OpenMM_GBSAOBCForce_setSoluteDielectric(OpenMM_GBSAOBCForce* target, double dielectric) {
    reinterpret_cast<GBSAOBCForce*>(target)->setSoluteDielectric(dielectric);
};


/* OpenMM::VariableVerletIntegrator*/
OpenMM_VariableVerletIntegrator* OpenMM_VariableVerletIntegrator_create(double errorTol) {
    return reinterpret_cast<OpenMM_VariableVerletIntegrator*>(new VariableVerletIntegrator(errorTol));
}

void OpenMM_VariableVerletIntegrator_destroy(OpenMM_VariableVerletIntegrator* target) {
    delete reinterpret_cast<VariableVerletIntegrator*>(target);
}
 
 
double OpenMM_VariableVerletIntegrator_getErrorTolerance(const OpenMM_VariableVerletIntegrator* target) {
    double result = reinterpret_cast<const VariableVerletIntegrator*>(target)->getErrorTolerance();
    return result;
};

void OpenMM_VariableVerletIntegrator_setErrorTolerance(OpenMM_VariableVerletIntegrator* target, double tol) {
    reinterpret_cast<VariableVerletIntegrator*>(target)->setErrorTolerance(tol);
};

void OpenMM_VariableVerletIntegrator_step(OpenMM_VariableVerletIntegrator* target, int steps) {
    reinterpret_cast<VariableVerletIntegrator*>(target)->step(steps);
};

void OpenMM_VariableVerletIntegrator_stepTo(OpenMM_VariableVerletIntegrator* target, double time) {
    reinterpret_cast<VariableVerletIntegrator*>(target)->stepTo(time);
};


/* OpenMM::CMMotionRemover*/
OpenMM_CMMotionRemover* OpenMM_CMMotionRemover_create(int frequency) {
    return reinterpret_cast<OpenMM_CMMotionRemover*>(new CMMotionRemover(frequency));
}

void OpenMM_CMMotionRemover_destroy(OpenMM_CMMotionRemover* target) {
    delete reinterpret_cast<CMMotionRemover*>(target);
}
 
 
int OpenMM_CMMotionRemover_getFrequency(const OpenMM_CMMotionRemover* target) {
    int result = reinterpret_cast<const CMMotionRemover*>(target)->getFrequency();
    return result;
};

void OpenMM_CMMotionRemover_setFrequency(OpenMM_CMMotionRemover* target, int freq) {
    reinterpret_cast<CMMotionRemover*>(target)->setFrequency(freq);
};


/* OpenMM::VerletIntegrator*/
OpenMM_VerletIntegrator* OpenMM_VerletIntegrator_create(double stepSize) {
    return reinterpret_cast<OpenMM_VerletIntegrator*>(new VerletIntegrator(stepSize));
}

void OpenMM_VerletIntegrator_destroy(OpenMM_VerletIntegrator* target) {
    delete reinterpret_cast<VerletIntegrator*>(target);
}
 
 
void OpenMM_VerletIntegrator_step(OpenMM_VerletIntegrator* target, int steps) {
    reinterpret_cast<VerletIntegrator*>(target)->step(steps);
};


/* OpenMM::RBTorsionForce*/
OpenMM_RBTorsionForce* OpenMM_RBTorsionForce_create() {
    return reinterpret_cast<OpenMM_RBTorsionForce*>(new RBTorsionForce());
}

void OpenMM_RBTorsionForce_destroy(OpenMM_RBTorsionForce* target) {
    delete reinterpret_cast<RBTorsionForce*>(target);
}
 
 
int OpenMM_RBTorsionForce_getNumTorsions(const OpenMM_RBTorsionForce* target) {
    int result = reinterpret_cast<const RBTorsionForce*>(target)->getNumTorsions();
    return result;
};

int OpenMM_RBTorsionForce_addTorsion(OpenMM_RBTorsionForce* target, int particle1, int particle2, int particle3, int particle4, double c0, double c1, double c2, double c3, double c4, double c5) {
    int result = reinterpret_cast<RBTorsionForce*>(target)->addTorsion(particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
    return result;
};

void OpenMM_RBTorsionForce_getTorsionParameters(const OpenMM_RBTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, double* c0, double* c1, double* c2, double* c3, double* c4, double* c5) {
    reinterpret_cast<const RBTorsionForce*>(target)->getTorsionParameters(index, *reinterpret_cast<int* >(particle1), *reinterpret_cast<int* >(particle2), *reinterpret_cast<int* >(particle3), *reinterpret_cast<int* >(particle4), *reinterpret_cast<double* >(c0), *reinterpret_cast<double* >(c1), *reinterpret_cast<double* >(c2), *reinterpret_cast<double* >(c3), *reinterpret_cast<double* >(c4), *reinterpret_cast<double* >(c5));
};

void OpenMM_RBTorsionForce_setTorsionParameters(OpenMM_RBTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, double c0, double c1, double c2, double c3, double c4, double c5) {
    reinterpret_cast<RBTorsionForce*>(target)->setTorsionParameters(index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
};


/* OpenMM::LangevinIntegrator*/
OpenMM_LangevinIntegrator* OpenMM_LangevinIntegrator_create(double temperature, double frictionCoeff, double stepSize) {
    return reinterpret_cast<OpenMM_LangevinIntegrator*>(new LangevinIntegrator(temperature, frictionCoeff, stepSize));
}

void OpenMM_LangevinIntegrator_destroy(OpenMM_LangevinIntegrator* target) {
    delete reinterpret_cast<LangevinIntegrator*>(target);
}
 
 
double OpenMM_LangevinIntegrator_getTemperature(const OpenMM_LangevinIntegrator* target) {
    double result = reinterpret_cast<const LangevinIntegrator*>(target)->getTemperature();
    return result;
};

void OpenMM_LangevinIntegrator_setTemperature(OpenMM_LangevinIntegrator* target, double temp) {
    reinterpret_cast<LangevinIntegrator*>(target)->setTemperature(temp);
};

double OpenMM_LangevinIntegrator_getFriction(const OpenMM_LangevinIntegrator* target) {
    double result = reinterpret_cast<const LangevinIntegrator*>(target)->getFriction();
    return result;
};

void OpenMM_LangevinIntegrator_setFriction(OpenMM_LangevinIntegrator* target, double coeff) {
    reinterpret_cast<LangevinIntegrator*>(target)->setFriction(coeff);
};

int OpenMM_LangevinIntegrator_getRandomNumberSeed(const OpenMM_LangevinIntegrator* target) {
    int result = reinterpret_cast<const LangevinIntegrator*>(target)->getRandomNumberSeed();
    return result;
};

void OpenMM_LangevinIntegrator_setRandomNumberSeed(OpenMM_LangevinIntegrator* target, int seed) {
    reinterpret_cast<LangevinIntegrator*>(target)->setRandomNumberSeed(seed);
};

void OpenMM_LangevinIntegrator_step(OpenMM_LangevinIntegrator* target, int steps) {
    reinterpret_cast<LangevinIntegrator*>(target)->step(steps);
};


/* OpenMM::Force*/
void OpenMM_Force_destroy(OpenMM_Force* target) {
    delete reinterpret_cast<Force*>(target);
}
 
 

/* OpenMM::HarmonicAngleForce*/
OpenMM_HarmonicAngleForce* OpenMM_HarmonicAngleForce_create() {
    return reinterpret_cast<OpenMM_HarmonicAngleForce*>(new HarmonicAngleForce());
}

void OpenMM_HarmonicAngleForce_destroy(OpenMM_HarmonicAngleForce* target) {
    delete reinterpret_cast<HarmonicAngleForce*>(target);
}
 
 
int OpenMM_HarmonicAngleForce_getNumAngles(const OpenMM_HarmonicAngleForce* target) {
    int result = reinterpret_cast<const HarmonicAngleForce*>(target)->getNumAngles();
    return result;
};

int OpenMM_HarmonicAngleForce_addAngle(OpenMM_HarmonicAngleForce* target, int particle1, int particle2, int particle3, double angle, double k) {
    int result = reinterpret_cast<HarmonicAngleForce*>(target)->addAngle(particle1, particle2, particle3, angle, k);
    return result;
};

void OpenMM_HarmonicAngleForce_getAngleParameters(const OpenMM_HarmonicAngleForce* target, int index, int* particle1, int* particle2, int* particle3, double* angle, double* k) {
    reinterpret_cast<const HarmonicAngleForce*>(target)->getAngleParameters(index, *reinterpret_cast<int* >(particle1), *reinterpret_cast<int* >(particle2), *reinterpret_cast<int* >(particle3), *reinterpret_cast<double* >(angle), *reinterpret_cast<double* >(k));
};

void OpenMM_HarmonicAngleForce_setAngleParameters(OpenMM_HarmonicAngleForce* target, int index, int particle1, int particle2, int particle3, double angle, double k) {
    reinterpret_cast<HarmonicAngleForce*>(target)->setAngleParameters(index, particle1, particle2, particle3, angle, k);
};


/* OpenMM::AndersenThermostat*/
OpenMM_AndersenThermostat* OpenMM_AndersenThermostat_create(double defaultTemperature, double defaultCollisionFrequency) {
    return reinterpret_cast<OpenMM_AndersenThermostat*>(new AndersenThermostat(defaultTemperature, defaultCollisionFrequency));
}

void OpenMM_AndersenThermostat_destroy(OpenMM_AndersenThermostat* target) {
    delete reinterpret_cast<AndersenThermostat*>(target);
}
 
 
const char* OpenMM_AndersenThermostat_Temperature() {
    const string* result = &AndersenThermostat::Temperature();
    return result->c_str();
};

const char* OpenMM_AndersenThermostat_CollisionFrequency() {
    const string* result = &AndersenThermostat::CollisionFrequency();
    return result->c_str();
};

double OpenMM_AndersenThermostat_getDefaultTemperature(const OpenMM_AndersenThermostat* target) {
    double result = reinterpret_cast<const AndersenThermostat*>(target)->getDefaultTemperature();
    return result;
};

double OpenMM_AndersenThermostat_getDefaultCollisionFrequency(const OpenMM_AndersenThermostat* target) {
    double result = reinterpret_cast<const AndersenThermostat*>(target)->getDefaultCollisionFrequency();
    return result;
};

int OpenMM_AndersenThermostat_getRandomNumberSeed(const OpenMM_AndersenThermostat* target) {
    int result = reinterpret_cast<const AndersenThermostat*>(target)->getRandomNumberSeed();
    return result;
};

void OpenMM_AndersenThermostat_setRandomNumberSeed(OpenMM_AndersenThermostat* target, int seed) {
    reinterpret_cast<AndersenThermostat*>(target)->setRandomNumberSeed(seed);
};


/* OpenMM::Platform*/
void OpenMM_Platform_destroy(OpenMM_Platform* target) {
    delete reinterpret_cast<Platform*>(target);
}
 
 
const char* OpenMM_Platform_getName(const OpenMM_Platform* target) {
    const string* result = &reinterpret_cast<const Platform*>(target)->getName();
    return result->c_str();
};

double OpenMM_Platform_getSpeed(const OpenMM_Platform* target) {
    double result = reinterpret_cast<const Platform*>(target)->getSpeed();
    return result;
};

OpenMM_Boolean OpenMM_Platform_supportsDoublePrecision(const OpenMM_Platform* target) {
    bool result = reinterpret_cast<const Platform*>(target)->supportsDoublePrecision();
    return (result ? OpenMM_True : OpenMM_False);
};

const OpenMM_StringArray* OpenMM_Platform_getPropertyNames(OpenMM_Platform* target) {
    const vector<string>* result = &reinterpret_cast<Platform*>(target)->getPropertyNames();
    return reinterpret_cast<const OpenMM_StringArray*>(result);
};

const char* OpenMM_Platform_getPropertyValue(const OpenMM_Platform* target, const OpenMM_Context* context, const char* property) {
    const string* result = &reinterpret_cast<const Platform*>(target)->getPropertyValue(*reinterpret_cast<const Context* >(context), string(property));
    return result->c_str();
};

void OpenMM_Platform_setPropertyValue(const OpenMM_Platform* target, OpenMM_Context* context, const char* property, const char* value) {
    reinterpret_cast<const Platform*>(target)->setPropertyValue(*reinterpret_cast<Context* >(context), string(property), string(value));
};

const char* OpenMM_Platform_getPropertyDefaultValue(const OpenMM_Platform* target, const char* property) {
    const string* result = &reinterpret_cast<const Platform*>(target)->getPropertyDefaultValue(string(property));
    return result->c_str();
};

void OpenMM_Platform_setPropertyDefaultValue(OpenMM_Platform* target, const char* property, const char* value) {
    reinterpret_cast<Platform*>(target)->setPropertyDefaultValue(string(property), string(value));
};

void OpenMM_Platform_contextCreated(const OpenMM_Platform* target, OpenMM_ContextImpl* context) {
    reinterpret_cast<const Platform*>(target)->contextCreated(*reinterpret_cast<ContextImpl* >(context));
};

void OpenMM_Platform_contextDestroyed(const OpenMM_Platform* target, OpenMM_ContextImpl* context) {
    reinterpret_cast<const Platform*>(target)->contextDestroyed(*reinterpret_cast<ContextImpl* >(context));
};

OpenMM_Boolean OpenMM_Platform_supportsKernels(const OpenMM_Platform* target, const OpenMM_StringArray* kernelNames) {
    bool result = reinterpret_cast<const Platform*>(target)->supportsKernels(*reinterpret_cast<const vector<string>* >(kernelNames));
    return (result ? OpenMM_True : OpenMM_False);
};

void OpenMM_Platform_registerPlatform(OpenMM_Platform* platform) {
    Platform::registerPlatform(reinterpret_cast<Platform* >(platform));
};

int OpenMM_Platform_getNumPlatforms() {
    int result = Platform::getNumPlatforms();
    return result;
};

OpenMM_Platform* OpenMM_Platform_getPlatform(int index) {
    Platform* result = &Platform::getPlatform(index);
    return reinterpret_cast<OpenMM_Platform*>(result);
};

OpenMM_Platform* OpenMM_Platform_findPlatform(const OpenMM_StringArray* kernelNames) {
    Platform* result = &Platform::findPlatform(*reinterpret_cast<const vector<string>* >(kernelNames));
    return reinterpret_cast<OpenMM_Platform*>(result);
};

void OpenMM_Platform_loadPluginLibrary(const char* file) {
    Platform::loadPluginLibrary(string(file));
};

const char* OpenMM_Platform_getDefaultPluginsDirectory() {
    const string* result = &Platform::getDefaultPluginsDirectory();
    return result->c_str();
};


/* OpenMM::State*/
void OpenMM_State_destroy(OpenMM_State* target) {
    delete reinterpret_cast<State*>(target);
}
 
 
double OpenMM_State_getTime(const OpenMM_State* target) {
    double result = reinterpret_cast<const State*>(target)->getTime();
    return result;
};

const OpenMM_Vec3Array* OpenMM_State_getPositions(const OpenMM_State* target) {
    const vector<Vec3>* result = &reinterpret_cast<const State*>(target)->getPositions();
    return reinterpret_cast<const OpenMM_Vec3Array*>(result);
};

const OpenMM_Vec3Array* OpenMM_State_getVelocities(const OpenMM_State* target) {
    const vector<Vec3>* result = &reinterpret_cast<const State*>(target)->getVelocities();
    return reinterpret_cast<const OpenMM_Vec3Array*>(result);
};

const OpenMM_Vec3Array* OpenMM_State_getForces(const OpenMM_State* target) {
    const vector<Vec3>* result = &reinterpret_cast<const State*>(target)->getForces();
    return reinterpret_cast<const OpenMM_Vec3Array*>(result);
};

double OpenMM_State_getKineticEnergy(const OpenMM_State* target) {
    double result = reinterpret_cast<const State*>(target)->getKineticEnergy();
    return result;
};

double OpenMM_State_getPotentialEnergy(const OpenMM_State* target) {
    double result = reinterpret_cast<const State*>(target)->getPotentialEnergy();
    return result;
};

const OpenMM_ParameterArray* OpenMM_State_getParameters(const OpenMM_State* target) {
    const map<string, double>* result = &reinterpret_cast<const State*>(target)->getParameters();
    return reinterpret_cast<const OpenMM_ParameterArray*>(result);
};


/* OpenMM::PeriodicTorsionForce*/
OpenMM_PeriodicTorsionForce* OpenMM_PeriodicTorsionForce_create() {
    return reinterpret_cast<OpenMM_PeriodicTorsionForce*>(new PeriodicTorsionForce());
}

void OpenMM_PeriodicTorsionForce_destroy(OpenMM_PeriodicTorsionForce* target) {
    delete reinterpret_cast<PeriodicTorsionForce*>(target);
}
 
 
int OpenMM_PeriodicTorsionForce_getNumTorsions(const OpenMM_PeriodicTorsionForce* target) {
    int result = reinterpret_cast<const PeriodicTorsionForce*>(target)->getNumTorsions();
    return result;
};

int OpenMM_PeriodicTorsionForce_addTorsion(OpenMM_PeriodicTorsionForce* target, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    int result = reinterpret_cast<PeriodicTorsionForce*>(target)->addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, k);
    return result;
};

void OpenMM_PeriodicTorsionForce_getTorsionParameters(const OpenMM_PeriodicTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, int* periodicity, double* phase, double* k) {
    reinterpret_cast<const PeriodicTorsionForce*>(target)->getTorsionParameters(index, *reinterpret_cast<int* >(particle1), *reinterpret_cast<int* >(particle2), *reinterpret_cast<int* >(particle3), *reinterpret_cast<int* >(particle4), *reinterpret_cast<int* >(periodicity), *reinterpret_cast<double* >(phase), *reinterpret_cast<double* >(k));
};

void OpenMM_PeriodicTorsionForce_setTorsionParameters(OpenMM_PeriodicTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    reinterpret_cast<PeriodicTorsionForce*>(target)->setTorsionParameters(index, particle1, particle2, particle3, particle4, periodicity, phase, k);
};


/* OpenMM::Integrator*/
void OpenMM_Integrator_destroy(OpenMM_Integrator* target) {
    delete reinterpret_cast<Integrator*>(target);
}
 
 
double OpenMM_Integrator_getStepSize(const OpenMM_Integrator* target) {
    double result = reinterpret_cast<const Integrator*>(target)->getStepSize();
    return result;
};

void OpenMM_Integrator_setStepSize(OpenMM_Integrator* target, double size) {
    reinterpret_cast<Integrator*>(target)->setStepSize(size);
};

double OpenMM_Integrator_getConstraintTolerance(const OpenMM_Integrator* target) {
    double result = reinterpret_cast<const Integrator*>(target)->getConstraintTolerance();
    return result;
};

void OpenMM_Integrator_setConstraintTolerance(OpenMM_Integrator* target, double tol) {
    reinterpret_cast<Integrator*>(target)->setConstraintTolerance(tol);
};

void OpenMM_Integrator_step(OpenMM_Integrator* target, int steps) {
    reinterpret_cast<Integrator*>(target)->step(steps);
};


/* OpenMM::System*/
OpenMM_System* OpenMM_System_create() {
    return reinterpret_cast<OpenMM_System*>(new System());
}

void OpenMM_System_destroy(OpenMM_System* target) {
    delete reinterpret_cast<System*>(target);
}
 
 
int OpenMM_System_getNumParticles(const OpenMM_System* target) {
    int result = reinterpret_cast<const System*>(target)->getNumParticles();
    return result;
};

int OpenMM_System_addParticle(OpenMM_System* target, double mass) {
    int result = reinterpret_cast<System*>(target)->addParticle(mass);
    return result;
};

double OpenMM_System_getParticleMass(const OpenMM_System* target, int index) {
    double result = reinterpret_cast<const System*>(target)->getParticleMass(index);
    return result;
};

void OpenMM_System_setParticleMass(OpenMM_System* target, int index, double mass) {
    reinterpret_cast<System*>(target)->setParticleMass(index, mass);
};

int OpenMM_System_getNumConstraints(const OpenMM_System* target) {
    int result = reinterpret_cast<const System*>(target)->getNumConstraints();
    return result;
};

int OpenMM_System_addConstraint(OpenMM_System* target, int particle1, int particle2, double distance) {
    int result = reinterpret_cast<System*>(target)->addConstraint(particle1, particle2, distance);
    return result;
};

void OpenMM_System_getConstraintParameters(const OpenMM_System* target, int index, int* particle1, int* particle2, double* distance) {
    reinterpret_cast<const System*>(target)->getConstraintParameters(index, *reinterpret_cast<int* >(particle1), *reinterpret_cast<int* >(particle2), *reinterpret_cast<double* >(distance));
};

void OpenMM_System_setConstraintParameters(OpenMM_System* target, int index, int particle1, int particle2, double distance) {
    reinterpret_cast<System*>(target)->setConstraintParameters(index, particle1, particle2, distance);
};

int OpenMM_System_addForce(OpenMM_System* target, OpenMM_Force* force) {
    int result = reinterpret_cast<System*>(target)->addForce(reinterpret_cast<Force* >(force));
    return result;
};

int OpenMM_System_getNumForces(const OpenMM_System* target) {
    int result = reinterpret_cast<const System*>(target)->getNumForces();
    return result;
};

OpenMM_Force* OpenMM_System_getForce(OpenMM_System* target, int index) {
    Force* result = &reinterpret_cast<System*>(target)->getForce(index);
    return reinterpret_cast<OpenMM_Force*>(result);
};

 
#if defined(__cplusplus)
}
#endif
