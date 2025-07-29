/* -----------------------------------------------------------------------------
 *             OpenMM(tm) HelloArgon example in C (June 2009)
 * -----------------------------------------------------------------------------
 * This program demonstrates a simple molecular simulation using the OpenMM
 * API for GPU-accelerated molecular dynamics simulation. The primary goal is
 * to make sure you can compile, link, and run with OpenMM and view the output.
 * The example is available in C++, C, and Fortran 95.
 *
 * The system modeled here is a small number of argon atoms in a vacuum.
 * A multi-frame PDB file is written to stdout which  can be read by VMD or 
 * other visualization tool to produce an animation of the resulting trajectory.
 * -------------------------------------------------------------------------- */

#include "OpenMMCWrapper.h"
#include <stdio.h>

/* Forward declaration of routine for printing one frame of the
   trajectory, defined later in this source file. */
void writePdbFrame(int frameNum, const OpenMM_State*);

void simulateArgon()
{
    OpenMM_System*         system;
    OpenMM_Integrator*     integrator;
    OpenMM_Context*        context;
    OpenMM_Platform*       platform;
    OpenMM_NonbondedForce* nonbond; 
    OpenMM_Vec3Array*      initPosInNm;
    OpenMM_StringArray*    pluginList;
    int                    a, frameNum;

    /* Load any shared libraries containing GPU implementations. */
    pluginList = OpenMM_Platform_loadPluginsFromDirectory(
        OpenMM_Platform_getDefaultPluginsDirectory());
    OpenMM_StringArray_destroy(pluginList);

    /* Create a system with nonbonded forces. System takes ownership
       of Force; don't destroy it yourself. */
    system  = OpenMM_System_create();
    nonbond = OpenMM_NonbondedForce_create(); 
    OpenMM_System_addForce(system, (OpenMM_Force*)nonbond);

    /* Create three atoms. */
    initPosInNm = OpenMM_Vec3Array_create(3);
    for (a = 0; a < 3; ++a) 
    {
        const OpenMM_Vec3 posNm = {0.5*a, 0, 0};     /*location, nm*/
        OpenMM_Vec3Array_set(initPosInNm, a, posNm);

        OpenMM_System_addParticle(system, 39.95); /*mass of Ar, grams/mole*/

        /* charge, L-J sigma (nm), well depth (kJ) */
        OpenMM_NonbondedForce_addParticle(nonbond, 0.0, 0.3350, 0.996); /*vdWRad(Ar)=.188 nm*/
    }

    /* Create particular integrator, and recast to generic one. */
    integrator = (OpenMM_Integrator*)OpenMM_VerletIntegrator_create(0.004); /*step size in ps*/

    /* Let OpenMM Context choose best platform. */
    context = OpenMM_Context_create(system, integrator);
    platform = OpenMM_Context_getPlatform(context);
    printf( "REMARK  Using OpenMM platform %s\n", 
        OpenMM_Platform_getName(platform));

    /* Set starting positions of the atoms. Leave time and velocity zero. */
    OpenMM_Context_setPositions(context, initPosInNm);

    /* Simulate. */
   for (frameNum=1; ;++frameNum) {
        /* Output current state information. */
        OpenMM_State* state    = OpenMM_Context_getState(context, OpenMM_State_Positions, 0);
        const double  timeInPs = OpenMM_State_getTime(state);
        writePdbFrame(frameNum, state); /*output coordinates*/
        OpenMM_State_destroy(state);

        if (timeInPs >= 10.)
            break;

        /* Advance state many steps at a time, for efficient use of OpenMM. */
        OpenMM_Integrator_step(integrator, 10); /*(use a lot more than this normally)*/
    }

    /* Free heap space for all the objects created above. */
    OpenMM_Vec3Array_destroy(initPosInNm);
    OpenMM_Context_destroy(context);
    OpenMM_Integrator_destroy(integrator);
    OpenMM_System_destroy(system);
}

int main() 
{
    simulateArgon();
    return 0;
}

/* Handy homebrew PDB writer for quick-and-dirty trajectory output. */
void writePdbFrame(int frameNum, const OpenMM_State* state) 
{
    int a;

    /* Reference atomic positions in the OpenMM State. */
    const OpenMM_Vec3Array* posInNm = OpenMM_State_getPositions(state);

    /* Use PDB MODEL cards to number trajectory frames. */
    printf("MODEL     %d\n", frameNum);  /*start of frame*/
    for (a = 0; a < OpenMM_Vec3Array_getSize(posInNm); ++a)
    {
        OpenMM_Vec3 posInAng;
        /* "10" here converts nanometers to Angstroms */
        posInAng = OpenMM_Vec3_scale(*OpenMM_Vec3Array_get(posInNm, a), 10.);
        printf("ATOM  %5d  AR   AR     1    ", a+1); /*atom number*/
        printf("%8.3f%8.3f%8.3f  1.00  0.00\n",      /*coordinates*/
            posInAng.x, posInAng.y, posInAng.z);
    }
    printf("ENDMDL\n"); /*end of frame*/
}
