/* --------------------------------------------------------------------------
 *       OpenMM(tm) example C wrapper function declarations (June 2009)
 * --------------------------------------------------------------------------
 * This header should be included by a C main program that would like to
 * access the OpenMM API through the C wrappers. Please note that this is not
 * an official part of OpenMM; it is just an example of how the C++ API can
 * be wrapped for access from C.
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
 * must persist from call to call during an OpenMM-powered simulation. */
typedef struct OpenMM_System_s      OpenMM_System;
typedef struct OpenMM_Integrator_s  OpenMM_Integrator;
typedef struct OpenMM_Context_s     OpenMM_Context;

/* This struct collects all the runtime object pointers together to 
 * facilitate use of an opaque handle in high-level C or Fortran code 
 * that doesn't have (or want) access to OpenMM declarations. This
 * does not have an equivalent in the OpenMM C++ API. */
typedef struct OpenMM_RuntimeObjects_s {
    OpenMM_System*      system;
    OpenMM_Integrator*  integrator;
    OpenMM_Context*     context;
} OpenMM_RuntimeObjects;

typedef double                      OpenMM_Vec3[3];
typedef struct OpenMM_Vec3Array_s   OpenMM_Vec3Array;
typedef struct OpenMM_String_s      OpenMM_String;

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
typedef struct OpenMM_Force_s               OpenMM_Force;
typedef struct OpenMM_NonbondedForce_s          OpenMM_NonbondedForce;
typedef struct OpenMM_GBSAOBCForce_s            OpenMM_GBSAOBCForce;

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
extern void OpenMM_System_destroy (OpenMM_System*);
extern void OpenMM_System_addForce(OpenMM_System*, OpenMM_Force*);
extern void OpenMM_System_addParticle(OpenMM_System*, double mass);


/* OpenMM::NonbondedForce */
extern OpenMM_NonbondedForce* OpenMM_NonbondedForce_create();
extern void OpenMM_NonbondedForce_destroy              (OpenMM_NonbondedForce*);
extern void OpenMM_NonbondedForce_setNonbondedMethod   (OpenMM_NonbondedForce*, 
                                                        OpenMM_NonbondedForce_NonbondedMethod);
extern void OpenMM_NonbondedForce_setCutoffDistance    (OpenMM_NonbondedForce*, double);
extern void OpenMM_NonbondedForce_setPeriodicBoxVectors(OpenMM_NonbondedForce*, 
                                                        const OpenMM_Vec3,const OpenMM_Vec3,const OpenMM_Vec3);
extern void OpenMM_NonbondedForce_addParticle(OpenMM_NonbondedForce*,
                                              double charge,
                                              double sigmaInNm,
                                              double vdwEnergyInKJ);

/* OpenMM::GBSAOBCForce */
extern OpenMM_GBSAOBCForce* OpenMM_GBSAOBCForce_create();
extern void OpenMM_GBSAOBCForce_destroy             (OpenMM_GBSAOBCForce*);
extern void OpenMM_GBSAOBCForce_setSolventDielectric(OpenMM_GBSAOBCForce*, double);
extern void OpenMM_GBSAOBCForce_setSoluteDielectric (OpenMM_GBSAOBCForce*, double);
extern void OpenMM_GBSAOBCForce_addParticle(OpenMM_GBSAOBCForce*,
                                            double charge,
                                            double radiusInNm,
                                            double scalingFactor);

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
extern double           OpenMM_Context_getTime(OpenMM_Context*);

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




