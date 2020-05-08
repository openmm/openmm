// Propagates a Nose Hoover chain a full timestep
KERNEL void propagateNoseHooverChain(GLOBAL mixed2* RESTRICT chainData, GLOBAL const mixed2 * RESTRICT energySum, GLOBAL mixed2* RESTRICT scaleFactor,
                                     GLOBAL mixed* RESTRICT chainMasses, GLOBAL mixed* RESTRICT chainForces, int chainType, int chainLength, int numMTS,
                                     int numDOFs, float timeStep, mixed kT, float frequency){
    const mixed kineticEnergy = chainType == 0 ? energySum[0].x : energySum[0].y;
    mixed scale = 1;
    if(kineticEnergy < 1e-8) return;
    for (int bead = 0; bead < chainLength; ++bead) chainMasses[bead] = kT / (frequency * frequency);
    chainMasses[0] *= numDOFs;
    mixed KE2 = 2.0f * kineticEnergy;
    mixed timeOverMTS = timeStep / numMTS;
    chainForces[0] = (KE2 - numDOFs * kT) / chainMasses[0];
    for (int bead = 0; bead < chainLength - 1; ++bead) {
        chainForces[bead + 1] = (chainMasses[bead] * chainData[bead].y * chainData[bead].y - kT) / chainMasses[bead + 1];
    }
    for (int mts = 0; mts < numMTS; ++mts) {
        BEGIN_YS_LOOP
            mixed wdt = ys * timeOverMTS;
            chainData[chainLength-1].y += 0.5f * wdt * chainForces[chainLength-1];
            for (int bead = chainLength - 2; bead >= 0; --bead) {
                mixed aa = exp(-0.25f * wdt * chainData[bead + 1].y);
                chainData[bead].y = aa * (chainData[bead].y * aa + 0.5f * wdt * chainForces[bead]);
            }
            // update particle velocities
            scale *= (mixed) exp(-wdt * chainData[0].y);;
            // update the thermostat positions
            for (int bead = 0; bead < chainLength; ++bead) {
                chainData[bead].x += chainData[bead].y * wdt;
            }
            // update the forces
            chainForces[0] = (scale * scale * KE2 - numDOFs * kT) / chainMasses[0];
            // update thermostat velocities
            for (int bead = 0; bead < chainLength - 1; ++bead) {
                mixed aa = exp(-0.25f * wdt * chainData[bead + 1].y);
                chainData[bead].y = aa * (aa * chainData[bead].y + 0.5f * wdt * chainForces[bead]);
                chainForces[bead + 1] = (chainMasses[bead] * chainData[bead].y * chainData[bead].y - kT) / chainMasses[bead + 1];
            }
            chainData[chainLength-1].y += 0.5f * wdt * chainForces[chainLength-1];
        END_YS_LOOP
    } // MTS loop
    if (chainType == 0) {
        scaleFactor[0].x = scale;
    } else {
        scaleFactor[0].y = scale;
    }
}


/**
 * Compute total (potential + kinetic) energy of the Nose-Hoover beads
 */
KERNEL void computeHeatBathEnergy(GLOBAL mixed* RESTRICT heatBathEnergy, int chainLength, int numDOFs,
                                  mixed kT, float frequency, GLOBAL const mixed2* RESTRICT chainData){
    // Note that this is always incremented; make sure it's zeroed properly before the first call
    for(int i = 0; i < chainLength; ++i) {
        mixed prefac = i ? 1 : numDOFs;
        mixed mass = prefac * kT / (frequency * frequency);
        mixed velocity = chainData[i].y; 
        // The kinetic energy of this bead
        heatBathEnergy[0] += 0.5f * mass * velocity * velocity;
        // The potential energy of this bead
        mixed position = chainData[i].x;
        heatBathEnergy[0] += prefac * kT * position;
    }
}

KERNEL void computeAtomsKineticEnergy(GLOBAL mixed2 * RESTRICT energyBuffer, int numAtoms,
                                      GLOBAL const mixed4* RESTRICT velm, GLOBAL const int *RESTRICT atoms){
    mixed2 energy = make_mixed2(0,0);
    int index = GLOBAL_ID;
    while (index < numAtoms){
        int atom = atoms[index];
        mixed4 v = velm[atom];
        mixed mass = v.w == 0 ? 0 : 1 / v.w;
        energy.x += 0.5f * mass * (v.x*v.x + v.y*v.y + v.z*v.z);
        index += GLOBAL_SIZE;
    }
    energyBuffer[GLOBAL_ID] = energy;
}

KERNEL void computePairsKineticEnergy(GLOBAL mixed2 * RESTRICT energyBuffer, int numPairs,
                                      GLOBAL const mixed4* RESTRICT velm, GLOBAL const int2 *RESTRICT pairs){
    mixed2 energy = make_mixed2(0,0);
    int index = GLOBAL_ID;
    while (index < numPairs){
        int2 pair = pairs[index];
        int atom1 = pair.x;
        int atom2 = pair.y;
        mixed4 v1 = velm[atom1];
        mixed4 v2 = velm[atom2];
        mixed m1 = v1.w == 0 ? 0 : 1 / v1.w;
        mixed m2 = v2.w == 0 ? 0 : 1 / v2.w;
        mixed4 cv;
        cv.x = (m1*v1.x + m2*v2.x) / (m1 + m2);
        cv.y = (m1*v1.y + m2*v2.y) / (m1 + m2);
        cv.z = (m1*v1.z + m2*v2.z) / (m1 + m2);
        mixed4 rv;
        rv.x = v2.x - v1.x;
        rv.y = v2.y - v1.y;
        rv.z = v2.z - v1.z;
        energy.x += 0.5f * (m1 + m2) * (cv.x*cv.x + cv.y*cv.y + cv.z*cv.z);
        energy.y += 0.5f * (m1 * m2 / (m1 + m2)) * (rv.x*rv.x + rv.y*rv.y + rv.z*rv.z);
        index += GLOBAL_SIZE;
    }
    // The atoms version of this has been called already, so accumulate instead of assigning here
    energyBuffer[GLOBAL_ID].x += energy.x;
    energyBuffer[GLOBAL_ID].y += energy.y;
}

KERNEL void scaleAtomsVelocities(GLOBAL mixed2* RESTRICT scaleFactor, int numAtoms,
                                   GLOBAL mixed4* RESTRICT velm, GLOBAL const int *RESTRICT atoms){
    const mixed scale = scaleFactor[0].x;
    int index = GLOBAL_ID;
    while (index < numAtoms){
        int atom = atoms[index];
        velm[atom].x *= scale;
        velm[atom].y *= scale;
        velm[atom].z *= scale;
        index += GLOBAL_SIZE;
    }
}

KERNEL void scalePairsVelocities(GLOBAL mixed2 * RESTRICT scaleFactor, int numPairs,
                                 GLOBAL mixed4* RESTRICT velm, GLOBAL const int2 *RESTRICT pairs){
    int index = GLOBAL_ID;
    mixed comScale = scaleFactor[0].x;
    mixed relScale = scaleFactor[0].y;
    while (index < numPairs){
        int atom1 = pairs[index].x;
        int atom2 = pairs[index].y;
        mixed m1 = velm[atom1].w == 0 ? 0 : 1 / velm[atom1].w;
        mixed m2 = velm[atom2].w == 0 ? 0 : 1 / velm[atom2].w;
        mixed4 cv;
        cv.x = (m1*velm[atom1].x + m2*velm[atom2].x) / (m1 + m2);
        cv.y = (m1*velm[atom1].y + m2*velm[atom2].y) / (m1 + m2);
        cv.z = (m1*velm[atom1].z + m2*velm[atom2].z) / (m1 + m2);
        mixed4 rv;
        rv.x = velm[atom2].x - velm[atom1].x;
        rv.y = velm[atom2].y - velm[atom1].y;
        rv.z = velm[atom2].z - velm[atom1].z;
        velm[atom1].x = comScale * cv.x - relScale * rv.x * m2 / (m1 + m2);
        velm[atom1].y = comScale * cv.y - relScale * rv.y * m2 / (m1 + m2);
        velm[atom1].z = comScale * cv.z - relScale * rv.z * m2 / (m1 + m2);
        velm[atom2].x = comScale * cv.x + relScale * rv.x * m1 / (m1 + m2);
        velm[atom2].y = comScale * cv.y + relScale * rv.y * m1 / (m1 + m2);
        velm[atom2].z = comScale * cv.z + relScale * rv.z * m1 / (m1 + m2);
        index += GLOBAL_SIZE;
    }
}

/**
 * Sum the energy buffer containing a pair of energies stored as mixed2.  This is taken from the analogous customIntegrator code
 */
KERNEL void reduceEnergyPair(GLOBAL const mixed2* RESTRICT sumBuffer, GLOBAL mixed2* result, int bufferSize) {
    LOCAL mixed2 tempBuffer[WORK_GROUP_SIZE];
    const unsigned int thread = LOCAL_ID;
    mixed2 sum = make_mixed2(0,0);
    for (unsigned int index = thread; index < bufferSize; index += LOCAL_SIZE) {
        sum.x += sumBuffer[index].x;
        sum.y += sumBuffer[index].y;
    }
    tempBuffer[thread].x = sum.x;
    tempBuffer[thread].y = sum.y;
    for (int i = 1; i < WORK_GROUP_SIZE; i *= 2) {
        SYNC_THREADS;
        if (thread%(i*2) == 0 && thread+i < WORK_GROUP_SIZE) {
            tempBuffer[thread].x += tempBuffer[thread+i].x;
            tempBuffer[thread].y += tempBuffer[thread+i].y;
        }
    }
    if (thread == 0)
        *result = tempBuffer[0];
}
