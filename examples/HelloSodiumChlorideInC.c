/* -----------------------------------------------------------------------------
 *            OpenMM(tm) HelloSodiumChloride example in C (May 2009)
 * -----------------------------------------------------------------------------
 * This is a complete, self-contained "hello world" example demonstrating 
 * GPU-accelerated constant energy simulation of a very simple system with just 
 * nonbonded forces, consisting of several sodium (Na+) and chloride (Cl-) ions. 
 * A multi-frame PDB file is written to stdout which can be read by VMD or other 
 * visualization tool to produce an animation of the resulting trajectory.
 *
 * Pay particular attention to the handling of units in this example. Incorrect
 * handling of units is a very common error; this example shows how you can
 * continue to work with Amber-style units of Angstroms and kCals while correctly
 * communicating with OpenMM in nanometers and kJoules.
 *
 * This example is written entirely in ANSI C, using a set of wrappers which
 * are NOT an official part of OpenMM.
 * -------------------------------------------------------------------------- */
#include "wrappers/OpenMM_CWrapper.h"

#include <stdio.h>

/* -----------------------------------------------------------------------------
 *                   MODELING AND SIMULATION PARAMETERS
 * -------------------------------------------------------------------------- */
const double StepSizeInFs        = 2;       // integration step size (fs)
const double ReportIntervalInFs  = 10;      // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 100;     // total simulation time (ps)

static void simulateNaCl();
static void writePDB(const OpenMM_Context*); // PDB file writer; see below.

/* -----------------------------------------------------------------------------
 *                                MAIN PROGRAM
 * ----------------------------------------------------------------------------- */
int main() {
    /* TODO: what about thrown exceptions? */
    OpenMM_Platform_loadPluginsFromDirectory
       (OpenMM_Platform_getDefaultPluginsDirectory());

    simulateNaCl();

    return 0;
}

/* -----------------------------------------------------------------------------
 *                          ATOM AND FORCE FIELD DATA
 * -----------------------------------------------------------------------------
 * This is not part of OpenMM; just a struct we can use to collect
 * atom parameters for this example. Normally atom parameters would
 * come from the force field's parameterization file.
 * We're going to use data in Angstrom and Kilocalorie units and
 * show how to safely convert to OpenMM's internal unit system
 * which uses nanometers and kilojoules.
 */
struct AtomInfo { 
    const char* pdb; 
    double      mass, charge, vdwRadiusInAng, vdwEnergyInKcal;
    OpenMM_Vec3 initPosInAngstroms;
} atoms[] = {
    /* pdb   mass   charge   vdwRadius  vdwEnergy   initPos */
    {" NA ", 22.99,   1,     1.8680,    0.00277,    { 8, 0,  0}},
    {" CL ", 35.45,  -1,     2.4700,    0.1000,     {-8, 0,  0}},
    {" NA ", 22.99,   1,     1.8680,    0.00277,    { 0, 9,  0}},
    {" CL ", 35.45,  -1,     2.4700,    0.1000,     { 0,-9,  0}},
    {" NA ", 22.99,   1,     1.8680,    0.00277,    { 0, 0,-10}},
    {" CL ", 35.45,  -1,     2.4700,    0.1000,     { 0, 0, 10}},
    {""} /* end of list */
};

/* -----------------------------------------------------------------------------
 *                              NaCl SIMULATION
 * ----------------------------------------------------------------------------- */
void simulateNaCl() {
    OpenMM_Integrator*  integrator;
    OpenMM_Context*     context;
    OpenMM_Vec3Array*   initialPositionsInNm;
    const OpenMM_Vec3   a = {5,0,0}, b = {0,5,0}, c = {0,0,5};
    const int           NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);
    int n;
    /* -------------------------------------------------------------------------
     * Create System and Force object. Add the Force to the System. The System
     * takes over ownership, so you should not destroy the Force object yourself.
     * ------------------------------------------------------------------------- */
    OpenMM_System*           system  = OpenMM_System_create();
    OpenMM_NonbondedForce*   nonbond = OpenMM_NonbondedForce_create();
    OpenMM_System_addForce(system, nonbond);

    /* -------------------------------------------------------------------------
     * Specify a periodic box and cutoff distance.
     * ------------------------------------------------------------------------- */
    OpenMM_NonbondedForce_setNonbondedMethod(nonbond, OpenMM_NonbondedForce_CutoffPeriodic);
    OpenMM_NonbondedForce_setCutoffDistance(nonbond, 2);
    OpenMM_NonbondedForce_setPeriodicBoxVectors(nonbond, a, b, c);

    /* -------------------------------------------------------------------------
     * Specify the atoms and their properties:
     *  (1) System needs to know the masses.
     *  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
     *  (3) Collect default positions for initializing the simulation later.
     * ------------------------------------------------------------------------- */
    initialPositionsInNm = OpenMM_Vec3Array_create(0);
    for (n=0; *atoms[n].pdb; ++n) {
        const struct AtomInfo* atom = &atoms[n];
        OpenMM_Vec3 posInNm;

        OpenMM_System_addParticle(system, atom->mass);
        OpenMM_NonbondedForce_addParticle
           (nonbond, atom->charge,
                     atom->vdwRadiusInAng  * OpenMM_NmPerAngstrom 
                                           * OpenMM_SigmaPerVdwRadius,
                     atom->vdwEnergyInKcal * OpenMM_KJPerKcal);

        /* Convert the initial position to nm and append to the array. */
        OpenMM_Vec3_scale(atom->initPosInAngstroms, OpenMM_NmPerAngstrom, posInNm);
        OpenMM_Vec3Array_append(initialPositionsInNm, posInNm);
    }

    /* -------------------------------------------------------------------------
     * Choose an Integrator for advancing time, and a Context connecting the
     * System with the Integrator for simulation. Let the Context choose the
     * best available Platform. Initialize the configuration from the default
     * positions we collected above. Initial velocities will be zero.
     * ------------------------------------------------------------------------- */
    integrator = (OpenMM_Integrator*)OpenMM_VerletIntegrator_create(StepSizeInFs * OpenMM_PsPerFs);
    context    = OpenMM_Context_create(system, integrator);
    OpenMM_Context_setPositions(context, initialPositionsInNm);

    /* -------------------------------------------------------------------------
     * Run the simulation:
     *  (1) Write the first line of the PDB file and the initial configuration.
     *  (2) Run silently entirely within OpenMM between reporting intervals.
     *  (3) Write a PDB frame when the time comes.
     * ------------------------------------------------------------------------- */
    printf( "REMARK  Using OpenMM platform %s\n", OpenMM_Context_getPlatform(context));
    writePDB(context);

    do {
        OpenMM_Integrator_step(integrator, NumSilentSteps);
        writePDB(context);
    } while (OpenMM_Context_getTime(context) < SimulationTimeInPs);

    /* Clean up top-level heap allocated objects that we're done with now. */
    OpenMM_Vec3Array_destroy(initialPositionsInNm);
    OpenMM_Context_destroy(context);
    OpenMM_Integrator_destroy(integrator);
}

/* -----------------------------------------------------------------------------
 *                               PDB FILE WRITER
 * ----------------------------------------------------------------------------- */
static void
writePDB(const OpenMM_Context* context) {
    static int modelFrameNumber = 0; /*numbering for MODEL records in pdb output*/

    /* Caution: at the moment asking for energy requires use of slow Reference 
       platform calculation. */

    /* Don't forget to destroy this State when you're done with it. */
    OpenMM_State* state = OpenMM_Context_createState(context,
        OpenMM_State_Positions | OpenMM_State_Velocities | OpenMM_State_Energy);

    const double energy =   OpenMM_State_getPotentialEnergy(state) 
                          + OpenMM_State_getKineticEnergy(state);

    /* Positions are maintained as a Vec3Array inside the State. This will give
     * us access, but don't destroy it yourself -- it will go away with the State.
     */
    const OpenMM_Vec3Array* posArray  = OpenMM_State_getPositions(state);
    const OpenMM_Vec3*      positions = OpenMM_Vec3Array_getAsVec3Ptr(posArray);
    const int               npos      = OpenMM_Vec3Array_size(posArray); 
    int i;

    /* write out in PDB format */
    modelFrameNumber++;
    printf("MODEL     %d\n", modelFrameNumber);
    printf("REMARK 250 time=%.3f picoseconds; Energy = %.3f kilojoules/mole\n", 
           OpenMM_State_getTime(state), energy);

    for (i=0; i < npos; ++i) {
        OpenMM_Vec3 pos;
        OpenMM_Vec3_scale(positions[i], OpenMM_AngstromsPerNm, pos);
        printf("ATOM  %5d %4s SLT     1    %8.3f%8.3f%8.3f  1.00  0.00            \n", 
            i+1, atoms[i].pdb, pos[0], pos[1], pos[2]);
    }
    printf("ENDMDL\n");

    OpenMM_State_destroy(state);
}

