
#ifndef OPENMM_CWRAPPER_H_
#define OPENMM_CWRAPPER_H_

/* Global Constants */
 
static const double OpenMM_KJPerKcal = 4.184;
static const double OpenMM_KcalPerKJ = 0.2390057361376673;
static const double OpenMM_PsPerFs = 0.001;
static const double OpenMM_AngstromsPerNm = 10;
static const double OpenMM_FsPerPs = 1000;
static const double OpenMM_RadiansPerDegree = 0.017453292519943295;
static const double OpenMM_NmPerAngstrom = 0.1;
static const double OpenMM_SigmaPerVdwRadius = 1.7817974362806785;
static const double OpenMM_VdwRadiusPerSigma = 0.5612310241546865;
static const double OpenMM_DegreesPerRadian = 57.29577951308232;

/* Type Declarations */
 
typedef struct OpenMM_HarmonicBondForce_struct OpenMM_HarmonicBondForce;
typedef struct OpenMM_BrownianIntegrator_struct OpenMM_BrownianIntegrator;
typedef struct OpenMM_OpenMMException_struct OpenMM_OpenMMException;
typedef struct OpenMM_NonbondedForce_struct OpenMM_NonbondedForce;
typedef struct OpenMM_VariableLangevinIntegrator_struct OpenMM_VariableLangevinIntegrator;
typedef struct OpenMM_GBVIForce_struct OpenMM_GBVIForce;
typedef struct OpenMM_Context_struct OpenMM_Context;
typedef struct OpenMM_GBSAOBCForce_struct OpenMM_GBSAOBCForce;
typedef struct OpenMM_VariableVerletIntegrator_struct OpenMM_VariableVerletIntegrator;
typedef struct OpenMM_CMMotionRemover_struct OpenMM_CMMotionRemover;
typedef struct OpenMM_VerletIntegrator_struct OpenMM_VerletIntegrator;
typedef struct OpenMM_ContextImpl_struct OpenMM_ContextImpl;
typedef struct OpenMM_RBTorsionForce_struct OpenMM_RBTorsionForce;
typedef struct OpenMM_LangevinIntegrator_struct OpenMM_LangevinIntegrator;
typedef struct OpenMM_Force_struct OpenMM_Force;
typedef struct OpenMM_HarmonicAngleForce_struct OpenMM_HarmonicAngleForce;
typedef struct OpenMM_AndersenThermostat_struct OpenMM_AndersenThermostat;
typedef struct OpenMM_ForceImpl_struct OpenMM_ForceImpl;
typedef struct OpenMM_Platform_struct OpenMM_Platform;
typedef struct OpenMM_State_struct OpenMM_State;
typedef struct OpenMM_PeriodicTorsionForce_struct OpenMM_PeriodicTorsionForce;
typedef struct OpenMM_Integrator_struct OpenMM_Integrator;
typedef struct OpenMM_System_struct OpenMM_System;
typedef struct OpenMM_Vec3Array_struct OpenMM_Vec3Array;
typedef struct OpenMM_StringArray_struct OpenMM_StringArray;
typedef struct OpenMM_BondArray_struct OpenMM_BondArray;
typedef struct OpenMM_ParameterArray_struct OpenMM_ParameterArray;
typedef struct {double x, y, z;} OpenMM_Vec3;

/* This struct collects all the runtime object pointers together to 
 * facilitate use of an opaque handle in high-level C or Fortran code 
 * that doesn't have (or want) access to OpenMM declarations. This
 * does not have an equivalent in the OpenMM C++ API. */
typedef struct OpenMM_RuntimeObjects_s {
    OpenMM_System*      system;
    OpenMM_Integrator*  integrator;
    OpenMM_Context*     context;
} OpenMM_RuntimeObjects;

typedef enum {OpenMM_False = 0, OpenMM_True = 1} OpenMM_Boolean;

#if defined(__cplusplus)
extern "C" {
#endif

/* OpenMM_Vec3 */
extern OpenMM_Vec3 OpenMM_Vec3_scale(const OpenMM_Vec3 vec, double scale);

/* OpenMM_Vec3Array */
extern OpenMM_Vec3Array* OpenMM_Vec3Array_create(int size);
extern void OpenMM_Vec3Array_destroy(OpenMM_Vec3Array* array);
extern int OpenMM_Vec3Array_getSize(const OpenMM_Vec3Array* array);
extern void OpenMM_Vec3Array_resize(OpenMM_Vec3Array* array, int size);
extern void OpenMM_Vec3Array_append(OpenMM_Vec3Array* array, const OpenMM_Vec3 vec);
extern void OpenMM_Vec3Array_set(OpenMM_Vec3Array* array, int index, const OpenMM_Vec3 vec);
extern const OpenMM_Vec3* OpenMM_Vec3Array_get(const OpenMM_Vec3Array* array, int index);

/* OpenMM_StringArray */
extern OpenMM_StringArray* OpenMM_StringArray_create(int size);
extern void OpenMM_StringArray_destroy(OpenMM_StringArray* array);
extern int OpenMM_StringArray_getSize(const OpenMM_StringArray* array);
extern void OpenMM_StringArray_resize(OpenMM_StringArray* array, int size);
extern void OpenMM_StringArray_append(OpenMM_StringArray* array, const char* string);
extern void OpenMM_StringArray_set(OpenMM_StringArray* array, int index, const char* string);
extern const char* OpenMM_StringArray_get(const OpenMM_StringArray* array, int index);

/* OpenMM_BondArray */
extern OpenMM_BondArray* OpenMM_BondArray_create(int size);
extern void OpenMM_BondArray_destroy(OpenMM_BondArray* array);
extern int OpenMM_BondArray_getSize(const OpenMM_BondArray* array);
extern void OpenMM_BondArray_resize(OpenMM_BondArray* array, int size);
extern void OpenMM_BondArray_append(OpenMM_BondArray* array, int particle1, int particle2);
extern void OpenMM_BondArray_set(OpenMM_BondArray* array, int index, int particle1, int particle2);
extern void OpenMM_BondArray_get(const OpenMM_BondArray* array, int index, int* particle1, int* particle2);

/* OpenMM_ParameterArray */
extern int OpenMM_ParameterArray_getSize(const OpenMM_ParameterArray* array);
extern double OpenMM_ParameterArray_get(const OpenMM_ParameterArray* array, const char* name);

/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
extern OpenMM_State* OpenMM_Context_getState(const OpenMM_Context* target, int types);
extern OpenMM_StringArray* OpenMM_Platform_loadPluginsFromDirectory(const char* directory);

 
 

/* OpenMM::HarmonicBondForce*/
extern OpenMM_HarmonicBondForce* OpenMM_HarmonicBondForce_create();
extern void OpenMM_HarmonicBondForce_destroy(OpenMM_HarmonicBondForce* target);
extern int OpenMM_HarmonicBondForce_getNumBonds(const OpenMM_HarmonicBondForce* target);
extern int OpenMM_HarmonicBondForce_addBond(OpenMM_HarmonicBondForce* target, int particle1, int particle2, double length, double k);
extern void OpenMM_HarmonicBondForce_getBondParameters(const OpenMM_HarmonicBondForce* target, int index, int* particle1, int* particle2, double* length, double* k);
extern void OpenMM_HarmonicBondForce_setBondParameters(OpenMM_HarmonicBondForce* target, int index, int particle1, int particle2, double length, double k);

/* OpenMM::BrownianIntegrator*/
extern OpenMM_BrownianIntegrator* OpenMM_BrownianIntegrator_create(double temperature, double frictionCoeff, double stepSize);
extern void OpenMM_BrownianIntegrator_destroy(OpenMM_BrownianIntegrator* target);
extern double OpenMM_BrownianIntegrator_getTemperature(const OpenMM_BrownianIntegrator* target);
extern void OpenMM_BrownianIntegrator_setTemperature(OpenMM_BrownianIntegrator* target, double temp);
extern double OpenMM_BrownianIntegrator_getFriction(const OpenMM_BrownianIntegrator* target);
extern void OpenMM_BrownianIntegrator_setFriction(OpenMM_BrownianIntegrator* target, double coeff);
extern int OpenMM_BrownianIntegrator_getRandomNumberSeed(const OpenMM_BrownianIntegrator* target);
extern void OpenMM_BrownianIntegrator_setRandomNumberSeed(OpenMM_BrownianIntegrator* target, int seed);
extern void OpenMM_BrownianIntegrator_step(OpenMM_BrownianIntegrator* target, int steps);

/* OpenMM::OpenMMException*/
extern OpenMM_OpenMMException* OpenMM_OpenMMException_create(const char* message);
extern void OpenMM_OpenMMException_destroy(OpenMM_OpenMMException* target);
extern const char* OpenMM_OpenMMException_what(const OpenMM_OpenMMException* target);

/* OpenMM::NonbondedForce*/
typedef enum {
  OpenMM_NonbondedForce_NoCutoff = 0, OpenMM_NonbondedForce_CutoffNonPeriodic = 1, OpenMM_NonbondedForce_CutoffPeriodic = 2, OpenMM_NonbondedForce_Ewald = 3, OpenMM_NonbondedForce_PME = 4
} OpenMM_NonbondedForce_NonbondedMethod;
extern OpenMM_NonbondedForce* OpenMM_NonbondedForce_create();
extern void OpenMM_NonbondedForce_destroy(OpenMM_NonbondedForce* target);
extern int OpenMM_NonbondedForce_getNumParticles(const OpenMM_NonbondedForce* target);
extern int OpenMM_NonbondedForce_getNumExceptions(const OpenMM_NonbondedForce* target);
extern OpenMM_NonbondedForce_NonbondedMethod OpenMM_NonbondedForce_getNonbondedMethod(const OpenMM_NonbondedForce* target);
extern void OpenMM_NonbondedForce_setNonbondedMethod(OpenMM_NonbondedForce* target, OpenMM_NonbondedForce_NonbondedMethod method);
extern double OpenMM_NonbondedForce_getCutoffDistance(const OpenMM_NonbondedForce* target);
extern void OpenMM_NonbondedForce_setCutoffDistance(OpenMM_NonbondedForce* target, double distance);
extern double OpenMM_NonbondedForce_getReactionFieldDielectric(const OpenMM_NonbondedForce* target);
extern void OpenMM_NonbondedForce_setReactionFieldDielectric(OpenMM_NonbondedForce* target, double dielectric);
extern double OpenMM_NonbondedForce_getEwaldErrorTolerance(const OpenMM_NonbondedForce* target);
extern void OpenMM_NonbondedForce_setEwaldErrorTolerance(OpenMM_NonbondedForce* target, double tol);
extern void OpenMM_NonbondedForce_getPeriodicBoxVectors(const OpenMM_NonbondedForce* target, OpenMM_Vec3* a, OpenMM_Vec3* b, OpenMM_Vec3* c);
extern void OpenMM_NonbondedForce_setPeriodicBoxVectors(OpenMM_NonbondedForce* target, OpenMM_Vec3 a, OpenMM_Vec3 b, OpenMM_Vec3 c);
extern int OpenMM_NonbondedForce_addParticle(OpenMM_NonbondedForce* target, double charge, double sigma, double epsilon);
extern void OpenMM_NonbondedForce_getParticleParameters(const OpenMM_NonbondedForce* target, int index, double* charge, double* sigma, double* epsilon);
extern void OpenMM_NonbondedForce_setParticleParameters(OpenMM_NonbondedForce* target, int index, double charge, double sigma, double epsilon);
extern int OpenMM_NonbondedForce_addException(OpenMM_NonbondedForce* target, int particle1, int particle2, double chargeProd, double sigma, double epsilon, OpenMM_Boolean replace);
extern void OpenMM_NonbondedForce_getExceptionParameters(const OpenMM_NonbondedForce* target, int index, int* particle1, int* particle2, double* chargeProd, double* sigma, double* epsilon);
extern void OpenMM_NonbondedForce_setExceptionParameters(OpenMM_NonbondedForce* target, int index, int particle1, int particle2, double chargeProd, double sigma, double epsilon);
extern void OpenMM_NonbondedForce_createExceptionsFromBonds(OpenMM_NonbondedForce* target, const OpenMM_BondArray* bonds, double coulomb14Scale, double lj14Scale);

/* OpenMM::VariableLangevinIntegrator*/
extern OpenMM_VariableLangevinIntegrator* OpenMM_VariableLangevinIntegrator_create(double temperature, double frictionCoeff, double errorTol);
extern void OpenMM_VariableLangevinIntegrator_destroy(OpenMM_VariableLangevinIntegrator* target);
extern double OpenMM_VariableLangevinIntegrator_getTemperature(const OpenMM_VariableLangevinIntegrator* target);
extern void OpenMM_VariableLangevinIntegrator_setTemperature(OpenMM_VariableLangevinIntegrator* target, double temp);
extern double OpenMM_VariableLangevinIntegrator_getFriction(const OpenMM_VariableLangevinIntegrator* target);
extern void OpenMM_VariableLangevinIntegrator_setFriction(OpenMM_VariableLangevinIntegrator* target, double coeff);
extern double OpenMM_VariableLangevinIntegrator_getErrorTolerance(const OpenMM_VariableLangevinIntegrator* target);
extern void OpenMM_VariableLangevinIntegrator_setErrorTolerance(OpenMM_VariableLangevinIntegrator* target, double tol);
extern int OpenMM_VariableLangevinIntegrator_getRandomNumberSeed(const OpenMM_VariableLangevinIntegrator* target);
extern void OpenMM_VariableLangevinIntegrator_setRandomNumberSeed(OpenMM_VariableLangevinIntegrator* target, int seed);
extern void OpenMM_VariableLangevinIntegrator_step(OpenMM_VariableLangevinIntegrator* target, int steps);
extern void OpenMM_VariableLangevinIntegrator_stepTo(OpenMM_VariableLangevinIntegrator* target, double time);

/* OpenMM::GBVIForce*/
extern OpenMM_GBVIForce* OpenMM_GBVIForce_create();
extern void OpenMM_GBVIForce_destroy(OpenMM_GBVIForce* target);
extern int OpenMM_GBVIForce_getNumParticles(const OpenMM_GBVIForce* target);
extern int OpenMM_GBVIForce_addParticle(OpenMM_GBVIForce* target, double charge, double radius, double gamma);
extern void OpenMM_GBVIForce_getParticleParameters(const OpenMM_GBVIForce* target, int index, double* charge, double* radius, double* gamma);
extern void OpenMM_GBVIForce_setParticleParameters(OpenMM_GBVIForce* target, int index, double charge, double radius, double gamma);
extern double OpenMM_GBVIForce_getSolventDielectric(const OpenMM_GBVIForce* target);
extern void OpenMM_GBVIForce_setSolventDielectric(OpenMM_GBVIForce* target, double dielectric);
extern double OpenMM_GBVIForce_getSoluteDielectric(const OpenMM_GBVIForce* target);
extern void OpenMM_GBVIForce_setSoluteDielectric(OpenMM_GBVIForce* target, double dielectric);

/* OpenMM::Context*/
extern OpenMM_Context* OpenMM_Context_create(OpenMM_System* system, OpenMM_Integrator* integrator);
extern OpenMM_Context* OpenMM_Context_create_2(OpenMM_System* system, OpenMM_Integrator* integrator, OpenMM_Platform* platform);
extern void OpenMM_Context_destroy(OpenMM_Context* target);
extern OpenMM_System* OpenMM_Context_getSystem(OpenMM_Context* target);
extern OpenMM_Integrator* OpenMM_Context_getIntegrator(OpenMM_Context* target);
extern OpenMM_Platform* OpenMM_Context_getPlatform(OpenMM_Context* target);
extern void OpenMM_Context_setTime(OpenMM_Context* target, double time);
extern void OpenMM_Context_setPositions(OpenMM_Context* target, const OpenMM_Vec3Array* positions);
extern void OpenMM_Context_setVelocities(OpenMM_Context* target, const OpenMM_Vec3Array* velocities);
extern double OpenMM_Context_getParameter(OpenMM_Context* target, const char* name);
extern void OpenMM_Context_setParameter(OpenMM_Context* target, const char* name, double value);
extern void OpenMM_Context_reinitialize(OpenMM_Context* target);

/* OpenMM::GBSAOBCForce*/
extern OpenMM_GBSAOBCForce* OpenMM_GBSAOBCForce_create();
extern void OpenMM_GBSAOBCForce_destroy(OpenMM_GBSAOBCForce* target);
extern int OpenMM_GBSAOBCForce_getNumParticles(const OpenMM_GBSAOBCForce* target);
extern int OpenMM_GBSAOBCForce_addParticle(OpenMM_GBSAOBCForce* target, double charge, double radius, double scalingFactor);
extern void OpenMM_GBSAOBCForce_getParticleParameters(const OpenMM_GBSAOBCForce* target, int index, double* charge, double* radius, double* scalingFactor);
extern void OpenMM_GBSAOBCForce_setParticleParameters(OpenMM_GBSAOBCForce* target, int index, double charge, double radius, double scalingFactor);
extern double OpenMM_GBSAOBCForce_getSolventDielectric(const OpenMM_GBSAOBCForce* target);
extern void OpenMM_GBSAOBCForce_setSolventDielectric(OpenMM_GBSAOBCForce* target, double dielectric);
extern double OpenMM_GBSAOBCForce_getSoluteDielectric(const OpenMM_GBSAOBCForce* target);
extern void OpenMM_GBSAOBCForce_setSoluteDielectric(OpenMM_GBSAOBCForce* target, double dielectric);

/* OpenMM::VariableVerletIntegrator*/
extern OpenMM_VariableVerletIntegrator* OpenMM_VariableVerletIntegrator_create(double errorTol);
extern void OpenMM_VariableVerletIntegrator_destroy(OpenMM_VariableVerletIntegrator* target);
extern double OpenMM_VariableVerletIntegrator_getErrorTolerance(const OpenMM_VariableVerletIntegrator* target);
extern void OpenMM_VariableVerletIntegrator_setErrorTolerance(OpenMM_VariableVerletIntegrator* target, double tol);
extern void OpenMM_VariableVerletIntegrator_step(OpenMM_VariableVerletIntegrator* target, int steps);
extern void OpenMM_VariableVerletIntegrator_stepTo(OpenMM_VariableVerletIntegrator* target, double time);

/* OpenMM::CMMotionRemover*/
extern OpenMM_CMMotionRemover* OpenMM_CMMotionRemover_create(int frequency);
extern void OpenMM_CMMotionRemover_destroy(OpenMM_CMMotionRemover* target);
extern int OpenMM_CMMotionRemover_getFrequency(const OpenMM_CMMotionRemover* target);
extern void OpenMM_CMMotionRemover_setFrequency(OpenMM_CMMotionRemover* target, int freq);

/* OpenMM::VerletIntegrator*/
extern OpenMM_VerletIntegrator* OpenMM_VerletIntegrator_create(double stepSize);
extern void OpenMM_VerletIntegrator_destroy(OpenMM_VerletIntegrator* target);
extern void OpenMM_VerletIntegrator_step(OpenMM_VerletIntegrator* target, int steps);

/* OpenMM::RBTorsionForce*/
extern OpenMM_RBTorsionForce* OpenMM_RBTorsionForce_create();
extern void OpenMM_RBTorsionForce_destroy(OpenMM_RBTorsionForce* target);
extern int OpenMM_RBTorsionForce_getNumTorsions(const OpenMM_RBTorsionForce* target);
extern int OpenMM_RBTorsionForce_addTorsion(OpenMM_RBTorsionForce* target, int particle1, int particle2, int particle3, int particle4, double c0, double c1, double c2, double c3, double c4, double c5);
extern void OpenMM_RBTorsionForce_getTorsionParameters(const OpenMM_RBTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, double* c0, double* c1, double* c2, double* c3, double* c4, double* c5);
extern void OpenMM_RBTorsionForce_setTorsionParameters(OpenMM_RBTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, double c0, double c1, double c2, double c3, double c4, double c5);

/* OpenMM::LangevinIntegrator*/
extern OpenMM_LangevinIntegrator* OpenMM_LangevinIntegrator_create(double temperature, double frictionCoeff, double stepSize);
extern void OpenMM_LangevinIntegrator_destroy(OpenMM_LangevinIntegrator* target);
extern double OpenMM_LangevinIntegrator_getTemperature(const OpenMM_LangevinIntegrator* target);
extern void OpenMM_LangevinIntegrator_setTemperature(OpenMM_LangevinIntegrator* target, double temp);
extern double OpenMM_LangevinIntegrator_getFriction(const OpenMM_LangevinIntegrator* target);
extern void OpenMM_LangevinIntegrator_setFriction(OpenMM_LangevinIntegrator* target, double coeff);
extern int OpenMM_LangevinIntegrator_getRandomNumberSeed(const OpenMM_LangevinIntegrator* target);
extern void OpenMM_LangevinIntegrator_setRandomNumberSeed(OpenMM_LangevinIntegrator* target, int seed);
extern void OpenMM_LangevinIntegrator_step(OpenMM_LangevinIntegrator* target, int steps);

/* OpenMM::Force*/
extern void OpenMM_Force_destroy(OpenMM_Force* target);

/* OpenMM::HarmonicAngleForce*/
extern OpenMM_HarmonicAngleForce* OpenMM_HarmonicAngleForce_create();
extern void OpenMM_HarmonicAngleForce_destroy(OpenMM_HarmonicAngleForce* target);
extern int OpenMM_HarmonicAngleForce_getNumAngles(const OpenMM_HarmonicAngleForce* target);
extern int OpenMM_HarmonicAngleForce_addAngle(OpenMM_HarmonicAngleForce* target, int particle1, int particle2, int particle3, double angle, double k);
extern void OpenMM_HarmonicAngleForce_getAngleParameters(const OpenMM_HarmonicAngleForce* target, int index, int* particle1, int* particle2, int* particle3, double* angle, double* k);
extern void OpenMM_HarmonicAngleForce_setAngleParameters(OpenMM_HarmonicAngleForce* target, int index, int particle1, int particle2, int particle3, double angle, double k);

/* OpenMM::AndersenThermostat*/
extern OpenMM_AndersenThermostat* OpenMM_AndersenThermostat_create(double defaultTemperature, double defaultCollisionFrequency);
extern void OpenMM_AndersenThermostat_destroy(OpenMM_AndersenThermostat* target);
extern const char* OpenMM_AndersenThermostat_Temperature();
extern const char* OpenMM_AndersenThermostat_CollisionFrequency();
extern double OpenMM_AndersenThermostat_getDefaultTemperature(const OpenMM_AndersenThermostat* target);
extern double OpenMM_AndersenThermostat_getDefaultCollisionFrequency(const OpenMM_AndersenThermostat* target);
extern int OpenMM_AndersenThermostat_getRandomNumberSeed(const OpenMM_AndersenThermostat* target);
extern void OpenMM_AndersenThermostat_setRandomNumberSeed(OpenMM_AndersenThermostat* target, int seed);

/* OpenMM::Platform*/
extern void OpenMM_Platform_destroy(OpenMM_Platform* target);
extern const char* OpenMM_Platform_getName(const OpenMM_Platform* target);
extern double OpenMM_Platform_getSpeed(const OpenMM_Platform* target);
extern OpenMM_Boolean OpenMM_Platform_supportsDoublePrecision(const OpenMM_Platform* target);
extern const OpenMM_StringArray* OpenMM_Platform_getPropertyNames(OpenMM_Platform* target);
extern const char* OpenMM_Platform_getPropertyValue(const OpenMM_Platform* target, const OpenMM_Context* context, const char* property);
extern void OpenMM_Platform_setPropertyValue(const OpenMM_Platform* target, OpenMM_Context* context, const char* property, const char* value);
extern const char* OpenMM_Platform_getPropertyDefaultValue(const OpenMM_Platform* target, const char* property);
extern void OpenMM_Platform_setPropertyDefaultValue(OpenMM_Platform* target, const char* property, const char* value);
extern void OpenMM_Platform_contextCreated(const OpenMM_Platform* target, OpenMM_ContextImpl* context);
extern void OpenMM_Platform_contextDestroyed(const OpenMM_Platform* target, OpenMM_ContextImpl* context);
extern OpenMM_Boolean OpenMM_Platform_supportsKernels(const OpenMM_Platform* target, const OpenMM_StringArray* kernelNames);
extern void OpenMM_Platform_registerPlatform(OpenMM_Platform* platform);
extern int OpenMM_Platform_getNumPlatforms();
extern OpenMM_Platform* OpenMM_Platform_getPlatform(int index);
extern OpenMM_Platform* OpenMM_Platform_findPlatform(const OpenMM_StringArray* kernelNames);
extern void OpenMM_Platform_loadPluginLibrary(const char* file);
extern const char* OpenMM_Platform_getDefaultPluginsDirectory();

/* OpenMM::State*/
typedef enum {
  OpenMM_State_Positions = 1, OpenMM_State_Velocities = 2, OpenMM_State_Forces = 4, OpenMM_State_Energy = 8, OpenMM_State_Parameters = 16
} OpenMM_State_DataType;
extern void OpenMM_State_destroy(OpenMM_State* target);
extern double OpenMM_State_getTime(const OpenMM_State* target);
extern const OpenMM_Vec3Array* OpenMM_State_getPositions(const OpenMM_State* target);
extern const OpenMM_Vec3Array* OpenMM_State_getVelocities(const OpenMM_State* target);
extern const OpenMM_Vec3Array* OpenMM_State_getForces(const OpenMM_State* target);
extern double OpenMM_State_getKineticEnergy(const OpenMM_State* target);
extern double OpenMM_State_getPotentialEnergy(const OpenMM_State* target);
extern const OpenMM_ParameterArray* OpenMM_State_getParameters(const OpenMM_State* target);

/* OpenMM::PeriodicTorsionForce*/
extern OpenMM_PeriodicTorsionForce* OpenMM_PeriodicTorsionForce_create();
extern void OpenMM_PeriodicTorsionForce_destroy(OpenMM_PeriodicTorsionForce* target);
extern int OpenMM_PeriodicTorsionForce_getNumTorsions(const OpenMM_PeriodicTorsionForce* target);
extern int OpenMM_PeriodicTorsionForce_addTorsion(OpenMM_PeriodicTorsionForce* target, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k);
extern void OpenMM_PeriodicTorsionForce_getTorsionParameters(const OpenMM_PeriodicTorsionForce* target, int index, int* particle1, int* particle2, int* particle3, int* particle4, int* periodicity, double* phase, double* k);
extern void OpenMM_PeriodicTorsionForce_setTorsionParameters(OpenMM_PeriodicTorsionForce* target, int index, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k);

/* OpenMM::Integrator*/
extern void OpenMM_Integrator_destroy(OpenMM_Integrator* target);
extern double OpenMM_Integrator_getStepSize(const OpenMM_Integrator* target);
extern void OpenMM_Integrator_setStepSize(OpenMM_Integrator* target, double size);
extern double OpenMM_Integrator_getConstraintTolerance(const OpenMM_Integrator* target);
extern void OpenMM_Integrator_setConstraintTolerance(OpenMM_Integrator* target, double tol);
extern void OpenMM_Integrator_step(OpenMM_Integrator* target, int steps);

/* OpenMM::System*/
extern OpenMM_System* OpenMM_System_create();
extern void OpenMM_System_destroy(OpenMM_System* target);
extern int OpenMM_System_getNumParticles(const OpenMM_System* target);
extern int OpenMM_System_addParticle(OpenMM_System* target, double mass);
extern double OpenMM_System_getParticleMass(const OpenMM_System* target, int index);
extern void OpenMM_System_setParticleMass(OpenMM_System* target, int index, double mass);
extern int OpenMM_System_getNumConstraints(const OpenMM_System* target);
extern int OpenMM_System_addConstraint(OpenMM_System* target, int particle1, int particle2, double distance);
extern void OpenMM_System_getConstraintParameters(const OpenMM_System* target, int index, int* particle1, int* particle2, double* distance);
extern void OpenMM_System_setConstraintParameters(OpenMM_System* target, int index, int particle1, int particle2, double distance);
extern int OpenMM_System_addForce(OpenMM_System* target, OpenMM_Force* force);
extern int OpenMM_System_getNumForces(const OpenMM_System* target);
extern OpenMM_Force* OpenMM_System_getForce(OpenMM_System* target, int index);

#if defined(__cplusplus)
}
#endif

#endif /*OPENMM_CWRAPPER_H_*/
