/* --------------------------------------------------------------------------
 *     OpenMM(tm) PROTOTYPE C wrapper function declarations (June 2009)
 * --------------------------------------------------------------------------
 * This header should be included by a C main program that would like to
 * access the OpenMM API through the C wrappers. Please note that this is an
 * experimental prototype, not an official part of OpenMM; it is just an 
 * example of how the C++ API can be wrapped for access from C.
 *
 * This set of wrappers is incomplete. If you add more, please send them
 * to us. Improvements in substance and style would also be greatly 
 * appreciated. If you have ideas (or better code) please post to the OpenMM 
 * forum on simtk.org/home/openmm or if you're shy you can email Michael 
 * Sherman at msherman@stanford.edu.
 *
 * Note: this header must be includable in both ANSI C and C++ code, because
 * the function declarations must be common to both the users and the
 * implementing code, which is in C++.
 * -------------------------------------------------------------------------- */

#ifndef OPENMM_CWRAPPER_H_
#define OPENMM_CWRAPPER_H_

/* Declare incomplete types corresponding to each of the OpenMM objects that
 * we want to make available. This allows us to have unique pointer types
 * for each object to maintain some semblance of type safety.
 */

/* These first three types represent the three OpenMM runtime objects that
 * must persist from call to call during an OpenMM-powered simulation. 
 * OpenMM_Integrator is the generic type of all Integrator objects. */
typedef struct OpenMM_System_s      OpenMM_System;
typedef struct OpenMM_Integrator_s  OpenMM_Integrator;
typedef struct OpenMM_Context_s     OpenMM_Context;

/* This is the generic type of all Force objects. */
typedef struct OpenMM_Force_s       OpenMM_Force;

/* This struct collects all the runtime object pointers together to 
 * facilitate use of an opaque handle in high-level C or Fortran code 
 * that doesn't have (or want) access to OpenMM declarations. This
 * does not have an equivalent in the OpenMM C++ API. */
typedef struct OpenMM_RuntimeObjects_s {
    OpenMM_System*      system;
    OpenMM_Integrator*  integrator;
    OpenMM_Context*     context;
} OpenMM_RuntimeObjects;

typedef double                          OpenMM_Vec3[3];
typedef struct OpenMM_Vec3Array_s       OpenMM_Vec3Array;
typedef struct OpenMM_BondArray_s   OpenMM_BondArray;
typedef struct OpenMM_String_s          OpenMM_String;

/*
 * OpenMM_Integrator is the generic type for all integrators. Cast your
 * specific pointer type to the generic one for communication with functions
 * that take type OpenMM_Integrator.
 */
typedef struct OpenMM_VerletIntegrator_s        OpenMM_VerletIntegrator;
typedef struct OpenMM_LangevinIntegrator_s      OpenMM_LangevinIntegrator;

/*
 * OpenMM_Force is the generic type for all Force objects. Create the
 * concrete Force object and then cast it to the generic type.
 */
typedef struct OpenMM_NonbondedForce_s          OpenMM_NonbondedForce;
typedef struct OpenMM_GBSAOBCForce_s            OpenMM_GBSAOBCForce;
typedef struct OpenMM_HarmonicBondForce_s       OpenMM_HarmonicBondForce;
typedef struct OpenMM_HarmonicAngleForce_s      OpenMM_HarmonicAngleForce;
typedef struct OpenMM_PeriodicTorsionForce_s    OpenMM_PeriodicTorsionForce;

typedef enum {
    OpenMM_NonbondedForce_NoCutoff            = 0,
    OpenMM_NonbondedForce_CutoffNonPeriodic   = 1,
    OpenMM_NonbondedForce_CutoffPeriodic      = 2,
    OpenMM_NonbondedForce_Ewald               = 3
} OpenMM_NonbondedForce_NonbondedMethod;


typedef struct OpenMM_State_s               OpenMM_State;
typedef enum {
    OpenMM_State_Positions   = 1, 
    OpenMM_State_Velocities  = 2,
    OpenMM_State_Forces      = 4,
    OpenMM_State_Energy      = 8,
    OpenMM_State_Parameters  = 16
} OpenMM_State_DataType;

/* Conversion constants from openmm/Units.h */

/*
 * The number of nanometers in an Angstrom.
 */
static const double OpenMM_NmPerAngstrom = 0.1;
/*
 * The number of Angstroms in a nanometer.
 */
static const double OpenMM_AngstromsPerNm = 10.0;
/*
 * The number of picoseconds in a femtosecond.
 */
static const double OpenMM_PsPerFs = 0.001;
/*
 * The number of femtoseconds in a picosecond.
 */
static const double OpenMM_FsPerPs = 1000.0;
/*
 * The number of kJ in a kcal.
 */
static const double OpenMM_KJPerKcal = 4.184;
/*
 * The number of kcal in a kJ.
 */
static const double OpenMM_KcalPerKJ = 1.0/4.184;
/*
 * The number of radians in a degree.
 */
static const double OpenMM_RadiansPerDegree = 3.1415926535897932385/180.0;
/*
 * The number of degrees in a radian.
 */
static const double OpenMM_DegreesPerRadian = 180.0/3.1415926535897932385;
/*
 * This is the conversion factor that takes you from a van der Waals radius
 * (defined as 1/2 the minimum energy separation) to the related Lennard Jones 
 * "sigma" parameter (defined as the zero crossing separation). The value
 * is 2*pow(2, -1/6).
 */
static const double OpenMM_SigmaPerVdwRadius = 1.78179743628068;


#if defined(__cplusplus)
extern "C" {
#endif

/* OpenMM_Vec3Array */
extern OpenMM_Vec3Array*    OpenMM_Vec3Array_create(int n);
extern int                  OpenMM_Vec3Array_size(const OpenMM_Vec3Array*);
extern void                 OpenMM_Vec3Array_resize(OpenMM_Vec3Array*, int n);
extern void                 OpenMM_Vec3Array_destroy(OpenMM_Vec3Array*);
extern void                 OpenMM_Vec3Array_append(OpenMM_Vec3Array*, const double[3]);
extern void                 OpenMM_Vec3Array_get(const OpenMM_Vec3Array*, int i, double[3]);
extern void                 OpenMM_Vec3Array_getScaled(const OpenMM_Vec3Array*, int i, double s, double[3]);
extern void                 OpenMM_Vec3Array_set(OpenMM_Vec3Array*, int i, const double[3]);
extern void                 OpenMM_Vec3Array_setScaled(OpenMM_Vec3Array*, int i, const double[3], double s);
extern void                 OpenMM_Vec3_scale(const double in[3], double s, double out[3]);

/* OpenMM_BondArray */
extern OpenMM_BondArray*    OpenMM_BondArray_create(int n);
extern int                  OpenMM_BondArray_size(const OpenMM_BondArray*);
extern void                 OpenMM_BondArray_resize(OpenMM_BondArray*, int n);
extern void                 OpenMM_BondArray_destroy(OpenMM_BondArray*);
extern void                 OpenMM_BondArray_append(OpenMM_BondArray*, int p1, int p2);
extern void                 OpenMM_BondArray_get(const OpenMM_BondArray*, int i, int* p1, int* p2);
extern void                 OpenMM_BondArray_set(OpenMM_BondArray*, int i, int p1, int p2);

/* OpenMM_String */
extern OpenMM_String*       OpenMM_String_create(const char* init);
extern void                 OpenMM_String_destroy(OpenMM_String*);
extern int                  OpenMM_String_length(const OpenMM_String*);
extern const char*          OpenMM_String_getAsC(const OpenMM_String*);
extern void                 OpenMM_String_get(const OpenMM_String*, char* buf, int buflen);
extern void                 OpenMM_String_set(OpenMM_String*, const char* buf);

/* OpenMM::Platform */
extern void OpenMM_Platform_loadPluginsFromDirectory(const char*);
extern const char* OpenMM_Platform_getDefaultPluginsDirectory();

/* OpenMM::System */
extern OpenMM_System* OpenMM_System_create();
extern void           OpenMM_System_destroy (OpenMM_System*);

extern int    OpenMM_System_addParticle(OpenMM_System*, double mass);
extern void   OpenMM_System_setParticleMass(OpenMM_System*, int ix, double mass);
extern double OpenMM_System_getParticleMass(const OpenMM_System*, int ix);

extern int  OpenMM_System_addConstraint(OpenMM_System*, int p1, int p2, double distance);
extern void OpenMM_System_setConstraintParameters(OpenMM_System*, int ix, 
                                                  int p1, int p2, double distance);
extern void OpenMM_System_getConstraintParameters(const OpenMM_System*, int ix, 
                                                  int* p1, int* p2, double* distance);

extern int                 OpenMM_System_addForce(OpenMM_System*, OpenMM_Force*);
extern OpenMM_Force*       OpenMM_System_updForce(OpenMM_System*, int ix);
extern const OpenMM_Force* OpenMM_System_getForce(const OpenMM_System*, int ix);

extern int OpenMM_System_getNumParticles(const OpenMM_System*);
extern int OpenMM_System_getNumConstraints(const OpenMM_System*);
extern int OpenMM_System_getNumForces(const OpenMM_System*);

/* OpenMM::NonbondedForce */
extern OpenMM_NonbondedForce* OpenMM_NonbondedForce_create();
extern void OpenMM_NonbondedForce_destroy              (OpenMM_NonbondedForce*);
extern void OpenMM_NonbondedForce_setNonbondedMethod   (OpenMM_NonbondedForce*, 
                                                        OpenMM_NonbondedForce_NonbondedMethod);
extern OpenMM_NonbondedForce_NonbondedMethod
            OpenMM_NonbondedForce_getNonbondedMethod   (const OpenMM_NonbondedForce*);
extern void OpenMM_NonbondedForce_setCutoffDistance    (OpenMM_NonbondedForce*, double);
extern double OpenMM_NonbondedForce_getCutoffDistance  (const OpenMM_NonbondedForce*);
extern void OpenMM_NonbondedForce_setPeriodicBoxVectors(OpenMM_NonbondedForce*, 
                                                        const OpenMM_Vec3,const OpenMM_Vec3,const OpenMM_Vec3);
extern void OpenMM_NonbondedForce_getPeriodicBoxVectors(const OpenMM_NonbondedForce*, 
                                                        OpenMM_Vec3, OpenMM_Vec3, OpenMM_Vec3);
extern int  OpenMM_NonbondedForce_addParticle(OpenMM_NonbondedForce*,
                                              double charge,
                                              double sigmaInNm,
                                              double vdwEnergyInKJ);
extern void OpenMM_NonbondedForce_setParticleParameters(OpenMM_NonbondedForce*, int index,
                                                        double charge,
                                                        double sigmaInNm,
                                                        double vdwEnergyInKJ);
extern void OpenMM_NonbondedForce_getParticleParameters(const OpenMM_NonbondedForce*, int index,
                                                        double* charge,
                                                        double* sigmaInNm,
                                                        double* vdwEnergyInKJ);
extern int OpenMM_NonbondedForce_getNumParticles(const OpenMM_NonbondedForce*);
extern int OpenMM_NonbondedForce_getNumExceptions(const OpenMM_NonbondedForce*);
extern int OpenMM_NonbondedForce_addException(OpenMM_NonbondedForce*, int p1, int p2, 
                                              double chargeProd, double sigma, double epsilon);
extern void OpenMM_NonbondedForce_getExceptionParameters(const OpenMM_NonbondedForce*, int index, int* p1, int* p2, 
                                                         double* chargeProd, double* sigma, double* epsilon);
extern void OpenMM_NonbondedForce_setExceptionParameters(OpenMM_NonbondedForce*, int index, int p1, int p2, 
                                                         double chargeProd, double sigma, double epsilon);
extern void OpenMM_NonbondedForces_createExceptionsFromBonds(OpenMM_NonbondedForce*, const OpenMM_BondArray*,
                                                             double coulomb14Scale, double lj14Scale);


/* OpenMM::GBSAOBCForce */
extern OpenMM_GBSAOBCForce* OpenMM_GBSAOBCForce_create();
extern void OpenMM_GBSAOBCForce_destroy             (OpenMM_GBSAOBCForce*);
extern void OpenMM_GBSAOBCForce_setSolventDielectric(OpenMM_GBSAOBCForce*, double);
extern void OpenMM_GBSAOBCForce_setSoluteDielectric (OpenMM_GBSAOBCForce*, double);
extern int  OpenMM_GBSAOBCForce_addParticle(OpenMM_GBSAOBCForce*,
                                            double charge,
                                            double radiusInNm,
                                            double scalingFactor);

/* OpenMM::HarmonicBondForce */
extern OpenMM_HarmonicBondForce* OpenMM_HarmonicBondForce_create();
extern void OpenMM_HarmonicBondForce_destroy(OpenMM_HarmonicBondForce*);
extern int OpenMM_HarmonicBondForce_getNumBonds(const OpenMM_HarmonicBondForce*);
extern int OpenMM_HarmonicBondForce_addBond(OpenMM_HarmonicBondForce*, int p1, int p2, double len, double k);
extern void OpenMM_HarmonicBondForce_getBondParameters(const OpenMM_HarmonicBondForce*, int ix, int* p1, int* p2, double* len, double* k);
extern void OpenMM_HarmonicBondForce_setBondParameters(OpenMM_HarmonicBondForce*, int ix, int p1, int p2, double len, double k);

/* OpenMM::HarmonicAngleForce */
extern OpenMM_HarmonicAngleForce* OpenMM_HarmonicAngleForce_create();
extern void OpenMM_HarmonicAngleForce_destroy(OpenMM_HarmonicAngleForce*);
extern int OpenMM_HarmonicAngleForce_getNumAngles(const OpenMM_HarmonicAngleForce*);
extern int OpenMM_HarmonicAngleForce_addAngle(OpenMM_HarmonicAngleForce*, int p1, int p2, int p3, double angle, double k);
extern void OpenMM_HarmonicAngleForce_getAngleParameters(const OpenMM_HarmonicAngleForce*, int ix, int* p1, int* p2, int* p3, double* angle, double* k);
extern void OpenMM_HarmonicAngleForce_setAngleParameters(OpenMM_HarmonicAngleForce*, int ix, int p1, int p2, int p3, double angle, double k);

/* OpenMM::PeriodicTorsionForce */
extern OpenMM_PeriodicTorsionForce* OpenMM_PeriodicTorsionForce_create();
extern void OpenMM_PeriodicTorsionForce_destroy(OpenMM_PeriodicTorsionForce*);
extern int OpenMM_PeriodicTorsionForce_getNumTorsions(const OpenMM_PeriodicTorsionForce*);
extern int OpenMM_PeriodicTorsionForce_addTorsion(OpenMM_PeriodicTorsionForce*, int p1, int p2, int p3, int p4,
                                                int periodicity, double phase, double k);
extern void OpenMM_PeriodicTorsionForce_getTorsionParameters(const OpenMM_PeriodicTorsionForce*, int ix, int* p1, int* p2, int* p3, int* p4,
                                                           int* periodicity, double* phase, double* k);
extern void OpenMM_PeriodicTorsionForce_setTorsionParameters(OpenMM_PeriodicTorsionForce*, int ix, int p1, int p2, int p3, int p4,
                                                           int periodicity, double phase, double k);

/* OpenMM::Integrator */
extern void OpenMM_Integrator_step(OpenMM_Integrator*, int numSteps);
extern void OpenMM_Integrator_destroy(OpenMM_Integrator*);
/* OpenMM::VerletIntegrator */
extern OpenMM_VerletIntegrator* OpenMM_VerletIntegrator_create(double stepSzInPs);
extern void                     OpenMM_VerletIntegrator_destroy(OpenMM_VerletIntegrator*);
extern void                     OpenMM_VerletIntegrator_step(OpenMM_VerletIntegrator*, int numSteps);
/* OpenMM::LangevinIntegrator */
extern OpenMM_LangevinIntegrator* OpenMM_LangevinIntegrator_create(double temperature, double frictionInPerPs, double stepSzInPs);
extern void                       OpenMM_LangevinIntegrator_destroy(OpenMM_LangevinIntegrator*);
extern void                       OpenMM_LangevinIntegrator_step(OpenMM_LangevinIntegrator*, int numSteps);

/* OpenMM::Context */
extern OpenMM_Context*  OpenMM_Context_create(OpenMM_System*, OpenMM_Integrator*);
extern void             OpenMM_Context_destroy(OpenMM_Context*);
extern void             OpenMM_Context_setPositions(OpenMM_Context*, const OpenMM_Vec3Array*);
extern void             OpenMM_Context_setVelocities(OpenMM_Context*, const OpenMM_Vec3Array*);
extern OpenMM_State*    OpenMM_Context_createState(const OpenMM_Context*, int types);
extern const char*      OpenMM_Context_getPlatformName(const OpenMM_Context*);

/* OpenMM::State */
extern void     OpenMM_State_destroy(OpenMM_State*);
extern double   OpenMM_State_getTime(const OpenMM_State*);
extern double   OpenMM_State_getPotentialEnergy(const OpenMM_State*);
extern double   OpenMM_State_getKineticEnergy(const OpenMM_State*);
extern const OpenMM_Vec3Array* OpenMM_State_getPositions(const OpenMM_State*);
extern const OpenMM_Vec3Array* OpenMM_State_getVelocities(const OpenMM_State*);

/* OpenMM_Runtime_Objects */
extern OpenMM_RuntimeObjects* OpenMM_RuntimeObjects_create();
extern void                   OpenMM_RuntimeObjects_clear(OpenMM_RuntimeObjects*);
extern void                   OpenMM_RuntimeObjects_destroy(OpenMM_RuntimeObjects*);
extern void                   OpenMM_RuntimeObjects_setSystem(OpenMM_RuntimeObjects*, OpenMM_System*);
extern void                   OpenMM_RuntimeObjects_setIntegrator(OpenMM_RuntimeObjects*, OpenMM_Integrator*);
extern void                   OpenMM_RuntimeObjects_setContext(OpenMM_RuntimeObjects*, OpenMM_Context*);
extern OpenMM_System*         OpenMM_RuntimeObjects_getSystem(OpenMM_RuntimeObjects*);
extern OpenMM_Integrator*     OpenMM_RuntimeObjects_getIntegrator(OpenMM_RuntimeObjects*);
extern OpenMM_Context*        OpenMM_RuntimeObjects_getContext(OpenMM_RuntimeObjects*);


#if defined(__cplusplus)
}
#endif

#endif /*OPENMM_CWRAPPER_H_*/




