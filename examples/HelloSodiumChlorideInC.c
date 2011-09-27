/* -----------------------------------------------------------------------------
 *           OpenMM(tm) HelloSodiumChloride example in C (June 2009)
 * -----------------------------------------------------------------------------
 * This is a complete, self-contained "hello world" example demonstrating 
 * GPU-accelerated constant temperature simulation of a very simple system with
 * just nonbonded forces, consisting of several sodium (Na+) and chloride (Cl-) 
 * ions in implicit solvent. A multi-frame PDB file is written to stdout which 
 * can be read by VMD or other visualization tool to produce an animation of the 
 * resulting trajectory.
 *
 * Pay particular attention to the handling of units in this example. Incorrect
 * handling of units is a very common error; this example shows how you can
 * continue to work with Amber-style units like Angstroms, kCals, and van der
 * Waals radii while correctly communicating with OpenMM in nm, kJ, and sigma.
 *
 * This example is written entirely in ANSI C, using the OpenMM C bindings.
 * -------------------------------------------------------------------------- */


#include <stdio.h>
#include <stdlib.h>

/* --------------------------------------------------------------------------
 *                   MODELING AND SIMULATION PARAMETERS
 * -------------------------------------------------------------------------- */
static const double Temperature         = 300;    /*Kelvins */
static const double FrictionInPerPs     = 91.;    /*collisions per ps*/
static const double SolventDielectric   = 80.;    /*typical for water    */
static const double SoluteDielectric    = 2.;     /*typical for protein  */

static const double StepSizeInFs        = 2;      /*integration step size (fs)  */
static const double ReportIntervalInFs  = 50;     /*how often for PDB frame (fs)*/
static const double SimulationTimeInPs  = 100;    /*total simulation time (ps)  */

static const int    WantEnergy          = 1;


/* --------------------------------------------------------------------------
 *                          ATOM AND FORCE FIELD DATA
 * --------------------------------------------------------------------------
 * This is not part of OpenMM; just a struct we can use to collect atom 
 * parameters for this example. Normally atom parameters would come from the 
 * force field's parameterization file. We're going to use data in Angstrom and 
 * Kilocalorie units and show how to safely convert to OpenMM's internal unit 
 * system which uses nanometers and kilojoules.
 */
typedef struct MyAtomInfo_s { 
    const char* pdb; 
    double      mass, charge, vdwRadiusInAng, vdwEnergyInKcal,
                gbsaRadiusInAng, gbsaScaleFactor;
    double      initPosInAng[3];
    double      posInAng[3]; /*leave room for runtime state info*/
} MyAtomInfo;
static MyAtomInfo atoms[] = {
/* pdb   mass  charge  vdwRad vdwEnergy   gbsaRad gbsaScale  initPos */
{" NA ", 22.99,  1,    1.8680, 0.00277,    1.992,   0.8,     8, 0,  0},
{" CL ", 35.45, -1,    2.4700, 0.1000,     1.735,   0.8,    -8, 0,  0},
{" NA ", 22.99,  1,    1.8680, 0.00277,    1.992,   0.8,     0, 9,  0},
{" CL ", 35.45, -1,    2.4700, 0.1000,     1.735,   0.8,     0,-9,  0},
{" NA ", 22.99,  1,    1.8680, 0.00277,    1.992,   0.8,     0, 0,-10},
{" CL ", 35.45, -1,    2.4700, 0.1000,     1.735,   0.8,     0, 0, 10},
{""} // end of list
};


/* --------------------------------------------------------------------------
 *                           INTERFACE TO OpenMM
 * --------------------------------------------------------------------------
 * These four functions and an opaque structure are used to interface our main
 * program with OpenMM without the main program having any direct interaction
 * with the OpenMM API. This is a clean approach for interfacing with any MD
 * code, although the details of the interface routines will differ.
 */
typedef struct MyOpenMMData_s MyOpenMMData;
static MyOpenMMData* myInitializeOpenMM(const MyAtomInfo atoms[],
                                        double temperature,
                                        double frictionInPs,
                                        double solventDielectric,
                                        double soluteDielectric,
                                        double stepSizeInFs, 
                                        const char** platformName);
static void          myStepWithOpenMM(MyOpenMMData*, int numSteps);
static void          myGetOpenMMState(MyOpenMMData*, int wantEnergy,
                                      double* time, double* energy, 
                                      MyAtomInfo atoms[]);
static void          myTerminateOpenMM(MyOpenMMData*);

/* --------------------------------------------------------------------------
 *                               PDB FILE WRITER
 * --------------------------------------------------------------------------
 * Given state data, output a single frame (pdb "model") of the trajectory.
 */
static void
myWritePDBFrame(int frameNum, double timeInPs, double energyInKcal, 
                const MyAtomInfo atoms[]) 
{
    int n;

    /* Write out in PDB format. */
    printf("MODEL     %d\n", frameNum);
    printf("REMARK 250 time=%.3f ps; energy=%.3f kcal/mole\n", 
           timeInPs, energyInKcal);
    for (n=0; *atoms[n].pdb; ++n)
        printf("ATOM  %5d %4s SLT     1    %8.3f%8.3f%8.3f  1.00  0.00\n", 
            n+1, atoms[n].pdb, 
            atoms[n].posInAng[0], atoms[n].posInAng[1], atoms[n].posInAng[2]);
    printf("ENDMDL\n");
}


/* --------------------------------------------------------------------------
 *                                MAIN PROGRAM
 * -------------------------------------------------------------------------- */
int main() {
    const int NumReports     = (int)(SimulationTimeInPs*1000 / ReportIntervalInFs + 0.5);
    const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);
    int frame;

    /* TODO: what about thrown exceptions? */
    double        time, energy;
    const char*   platformName;

    /* Set up OpenMM data structures; returns OpenMM Platform name. */
    MyOpenMMData* omm = myInitializeOpenMM(atoms, Temperature, FrictionInPerPs,
                                           SolventDielectric, SoluteDielectric,
                                           StepSizeInFs, &platformName);

    /* Run the simulation:
     *  (1) Write the first line of the PDB file and the initial configuration.
     *  (2) Run silently entirely within OpenMM between reporting intervals.
     *  (3) Write a PDB frame when the time comes. */
    printf("REMARK  Using OpenMM platform %s\n", platformName);
    myGetOpenMMState(omm, WantEnergy, &time, &energy, atoms);
    myWritePDBFrame(0, time, energy, atoms);

    for (frame=1; frame <= NumReports; ++frame) {
        myStepWithOpenMM(omm, NumSilentSteps);
        myGetOpenMMState(omm, WantEnergy, &time, &energy, atoms);
        myWritePDBFrame(frame, time, energy, atoms);
    } 

    /* Clean up OpenMM data structures. */
    myTerminateOpenMM(omm);

    return 0; /* Normal return from main. */
}


/* --------------------------------------------------------------------------
 *                           OpenMM-USING CODE
 * --------------------------------------------------------------------------
 * The OpenMM C-wrapped API is visible only at this point and below. Normally 
 * this would be in a separate compilation module; we're including it here for 
 * simplicity. We suggest that you write them in C++ if possible; in fact you
 * can use the implementation from the C++ version of this example if you 
 * want. However, the methods are reimplemented in C below in case you prefer.
 */
#include "OpenMMCWrapper.h"


struct MyOpenMMData_s {
    OpenMM_System*      system;
    OpenMM_Context*     context;
    OpenMM_Integrator*  integrator;
};

/* --------------------------------------------------------------------------
 *                      INITIALIZE OpenMM DATA STRUCTURES
 * --------------------------------------------------------------------------
 * We take these actions here:
 * (1) Load any available OpenMM plugins, e.g. Cuda and Brook.
 * (2) Allocate a MyOpenMMData structure to hang on to OpenMM data structures
 *     in a manner which is opaque to the caller.
 * (3) Fill the OpenMM::System with the force field parameters we want to
 *     use and the particular set of atoms to be simulated.
 * (4) Create an Integrator and a Context associating the Integrator with
 *     the System.
 * (5) Select the OpenMM platform to be used.
 * (6) Return an opaque pointer to the MyOpenMMData struct and the name 
 *     of the Platform in use.
 *
 * Note that this function must understand the calling MD code's molecule and
 * force field data structures so will need to be customized for each MD code.
 */
static MyOpenMMData* 
myInitializeOpenMM( const MyAtomInfo    atoms[],
                    double              temperature,
                    double              frictionInPerPs,
                    double              solventDielectric,
                    double              soluteDielectric,
                    double              stepSizeInFs, 
                    const char**        platformName) 
{
    /* Allocate space to hold OpenMM objects while we're using them. */
    MyOpenMMData* omm = (MyOpenMMData*)malloc(sizeof(struct MyOpenMMData_s));
    /* These are temporary OpenMM objects used and discarded here. */
    OpenMM_Vec3Array*       initialPosInNm;
    OpenMM_StringArray*     pluginList;
    OpenMM_NonbondedForce*  nonbond;
    OpenMM_GBSAOBCForce*    gbsa;
    OpenMM_Platform*        platform;
    int n;

    /* Load all available OpenMM plugins from their default location. */
    pluginList = OpenMM_Platform_loadPluginsFromDirectory
       (OpenMM_Platform_getDefaultPluginsDirectory());
    OpenMM_StringArray_destroy(pluginList);

    /* Create a System and Force objects within the System. Retain a reference
     * to each force object so we can fill in the forces. Note: the OpenMM
     * System takes ownership of the force objects; don't delete them yourself. */
    omm->system = OpenMM_System_create();
    nonbond     = OpenMM_NonbondedForce_create();
    gbsa        = OpenMM_GBSAOBCForce_create();
    OpenMM_System_addForce(omm->system, (OpenMM_Force*)nonbond);
    OpenMM_System_addForce(omm->system, (OpenMM_Force*)gbsa);

    /* Specify dielectrics for GBSA implicit solvation. */
    OpenMM_GBSAOBCForce_setSolventDielectric(gbsa, solventDielectric);
    OpenMM_GBSAOBCForce_setSoluteDielectric(gbsa, soluteDielectric);

    /* Specify the atoms and their properties:
     *  (1) System needs to know the masses.
     *  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
     *  (3) GBSA needs charge, radius, and scale factor.
     *  (4) Collect default positions for initializing the simulation later. */
    initialPosInNm = OpenMM_Vec3Array_create(0);
    for (n=0; *atoms[n].pdb; ++n) {
        const MyAtomInfo* atom = &atoms[n];
        OpenMM_Vec3 posInNm;

        OpenMM_System_addParticle(omm->system, atom->mass);

        OpenMM_NonbondedForce_addParticle(nonbond,
                                atom->charge,
                                atom->vdwRadiusInAng  * OpenMM_NmPerAngstrom 
                                                      * OpenMM_SigmaPerVdwRadius,
                                atom->vdwEnergyInKcal * OpenMM_KJPerKcal);

        OpenMM_GBSAOBCForce_addParticle(gbsa,
                                atom->charge, 
                                atom->gbsaRadiusInAng * OpenMM_NmPerAngstrom,
                                atom->gbsaScaleFactor);

        /* Convert the initial position to nm and append to the array. */
        posInNm = OpenMM_Vec3_scale(*(const OpenMM_Vec3*)atom->initPosInAng, OpenMM_NmPerAngstrom);
        OpenMM_Vec3Array_append(initialPosInNm, posInNm);
    }

    /* Choose an Integrator for advancing time, and a Context connecting the
     * System with the Integrator for simulation. Let the Context choose the
     * best available Platform. Initialize the configuration from the default
     * positions we collected above. Initial velocities will be zero but could
     * have been set here. */
    omm->integrator = (OpenMM_Integrator*)OpenMM_LangevinIntegrator_create(
                                            temperature, frictionInPerPs, 
                                            stepSizeInFs * OpenMM_PsPerFs);
    omm->context    = OpenMM_Context_create(omm->system, omm->integrator);
    OpenMM_Context_setPositions(omm->context, initialPosInNm);

    platform = OpenMM_Context_getPlatform(omm->context);
    *platformName = OpenMM_Platform_getName(platform);
    return omm;
}




/* --------------------------------------------------------------------------
 *                    COPY STATE BACK TO CPU FROM OPENMM
 * -------------------------------------------------------------------------- */
static void
myGetOpenMMState(MyOpenMMData* omm, int wantEnergy, 
                 double* timeInPs, double* energyInKcal,
                 MyAtomInfo atoms[])
{
    OpenMM_State*           state;
    const OpenMM_Vec3Array* posArrayInNm;
    int                     infoMask;
    int n;

    infoMask = OpenMM_State_Positions;
    if (wantEnergy) {
        infoMask += OpenMM_State_Velocities; /*for kinetic energy (cheap)*/
        infoMask += OpenMM_State_Energy;     /*for pot. energy (expensive)*/
    }
    /* Forces are also available (and cheap). */

    /* State object is created here and must be explicitly destroyed below. */
    state = OpenMM_Context_getState(omm->context, infoMask, 0);
    *timeInPs = OpenMM_State_getTime(state); /* OpenMM time is in ps already. */

    /* Positions are maintained as a Vec3Array inside the State. This will give
     * us access, but don't destroy it yourself -- it will go away with the State. */
    posArrayInNm = OpenMM_State_getPositions(state);
    for (n=0; *atoms[n].pdb; ++n)
        /* Sets atoms[n].pos = posArray[n] * Angstroms/nm. */
        *(OpenMM_Vec3*)atoms[n].posInAng = 
            OpenMM_Vec3_scale(*OpenMM_Vec3Array_get(posArrayInNm, n),
                              OpenMM_AngstromsPerNm); 

    /* If energy has been requested, obtain it and convert from kJ to kcal. */
    *energyInKcal = 0;
    if (wantEnergy) 
        *energyInKcal = (   OpenMM_State_getPotentialEnergy(state) 
                          + OpenMM_State_getKineticEnergy(state))
                        * OpenMM_KcalPerKJ;

    OpenMM_State_destroy(state);
}


// -----------------------------------------------------------------------------
//                     TAKE MULTIPLE STEPS USING OpenMM 
// -----------------------------------------------------------------------------
static void 
myStepWithOpenMM(MyOpenMMData* omm, int numSteps) {
    OpenMM_Integrator_step(omm->integrator, numSteps);
}

// -----------------------------------------------------------------------------
//                     DEALLOCATE OpenMM OBJECTS
// -----------------------------------------------------------------------------
static void 
myTerminateOpenMM(MyOpenMMData* omm) {
    /* Clean up top-level heap allocated objects that we're done with now. */
    OpenMM_Context_destroy(omm->context);
    OpenMM_Integrator_destroy(omm->integrator);
    OpenMM_System_destroy(omm->system);
    free(omm);
}


