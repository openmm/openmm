/**
 * Scale the particle positions with each axis independent.
 */

KERNEL void scalePositions(float scaleX, float scaleY, float scaleZ, int numMolecules, real4 periodicBoxSize,
        real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, GLOBAL real4* RESTRICT posq,
        GLOBAL const int* RESTRICT moleculeAtoms, GLOBAL const int* RESTRICT moleculeStartIndex) {
    for (int index = GLOBAL_ID; index < numMolecules; index += GLOBAL_SIZE) {
        int first = moleculeStartIndex[index];
        int last = moleculeStartIndex[index+1];
        int numAtoms = last-first;

        // Find the center of each molecule.

        real3 center = make_real3(0, 0, 0);
        for (int atom = first; atom < last; atom++) {
            real4 pos = posq[moleculeAtoms[atom]];
            center.x += pos.x;
            center.y += pos.y;
            center.z += pos.z;
        }
        real invNumAtoms = RECIP((real) numAtoms);
        center.x *= invNumAtoms;
        center.y *= invNumAtoms;
        center.z *= invNumAtoms;

        // Now scale the position of the molecule center.

        real3 delta;
        delta.x = center.x*(scaleX-1);
        delta.y = center.y*(scaleY-1);
        delta.z = center.z*(scaleZ-1);
        for (int atom = first; atom < last; atom++) {
            real4 pos = posq[moleculeAtoms[atom]];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            posq[moleculeAtoms[atom]] = pos;
        }
    }
}

/**
 * Compute the kinetic energy of molecular centers of mass.  Depending on the value
 * of the COMPONENTS macro, this may compute either 1) the total kinetic energy,
 * 2) three components of kinetic energy corresponding to the three axes, or
 * 3) six components corresponding to all elements of the box tensor.
 */
KERNEL void computeMolecularKineticEnergy(int numMolecules, GLOBAL mixed4* RESTRICT velm, GLOBAL const int* RESTRICT moleculeAtoms,
        GLOBAL const int* RESTRICT moleculeStartIndex,
#if COMPONENTS == 1
        GLOBAL mixed* RESTRICT energyBuffer
#else
        GLOBAL mixed* RESTRICT energyBuffer1, GLOBAL mixed* RESTRICT energyBuffer2, GLOBAL mixed* RESTRICT energyBuffer3
#if COMPONENTS == 6
        , GLOBAL mixed* RESTRICT energyBuffer4, GLOBAL mixed* RESTRICT energyBuffer5, GLOBAL mixed* RESTRICT energyBuffer6
#endif
#endif
    ) {
#if COMPONENTS == 1
    GLOBAL mixed* buffers[] = {energyBuffer};
#elif COMPONENTS == 3
    GLOBAL mixed* buffers[] = {energyBuffer1, energyBuffer2, energyBuffer3};
#else
    GLOBAL mixed* buffers[] = {energyBuffer1, energyBuffer2, energyBuffer3, energyBuffer4, energyBuffer5, energyBuffer6};
#endif
    mixed ke[COMPONENTS];
    for (int i = 0; i < COMPONENTS; i++)
        ke[i] = 0;
    for (int index = GLOBAL_ID; index < numMolecules; index += GLOBAL_SIZE) {
        int first = moleculeStartIndex[index];
        int last = moleculeStartIndex[index+1];
        int numAtoms = last-first;
        mixed3 molVel = make_mixed3(0);
        float molMass = 0;
        for (int atom = first; atom < last; atom++) {
            mixed4 v = velm[moleculeAtoms[atom]];
            float mass = (v.w == 0 ? 0 : 1/v.w);
            molVel += mass*trimTo3(v);
            molMass += mass;
        }
        molVel *= RECIP((mixed) numAtoms);
#if COMPONENTS == 1
        ke[0] += 0.5f*molMass*dot(molVel, molVel);
#else
        ke[0] += 0.5f*molMass*molVel.x*molVel.x;
        ke[1] += 0.5f*molMass*molVel.y*molVel.y;
        ke[2] += 0.5f*molMass*molVel.z*molVel.z;
#if COMPONENTS == 6
        ke[3] += 0.5f*molMass*molVel.x*molVel.y;
        ke[4] += 0.5f*molMass*molVel.x*molVel.z;
        ke[5] += 0.5f*molMass*molVel.y*molVel.z;
#endif
#endif
    }

    // Sum the contributions from all the threads in this block and write the result.

    LOCAL mixed tempBuffer[WORK_GROUP_SIZE];
    for (int j = 0; j < COMPONENTS; j++) {
        tempBuffer[LOCAL_ID] = ke[j];
        for (int i = 1; i < WORK_GROUP_SIZE; i *= 2) {
            SYNC_THREADS;
            if (LOCAL_ID%(i*2) == 0 && LOCAL_ID+i < WORK_GROUP_SIZE)
                tempBuffer[LOCAL_ID] += tempBuffer[LOCAL_ID+i];
        }
        if (LOCAL_ID == 0)
            buffers[j][GROUP_ID] = tempBuffer[0];
        SYNC_THREADS;
    }
}
