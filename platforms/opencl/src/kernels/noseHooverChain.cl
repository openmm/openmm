
//#include <initializer_list>

__kernel void propagateNoseHooverChain(__global mixed2* restrict chainData, __global const mixed2 * restrict energySum, __global mixed2* restrict scaleFactor,
                                                    __global mixed* restrict chainMasses, __global mixed* restrict chainForces, 
                                                    int chainType, int chainLength, int numMTS, int numDOFs, float timeStep,
                                                    mixed kT, float frequency){
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
            chainData[chainLength-1].y += 0.25f * wdt * chainForces[chainLength-1];
            for (int bead = chainLength - 2; bead >= 0; --bead) {
                mixed aa = EXP(-0.125f * wdt * chainData[bead + 1].y);
                chainData[bead].y = aa * (chainData[bead].y * aa + 0.25f * wdt * chainForces[bead]);
            }
            // update particle velocities
            mixed aa = EXP(-0.5f * wdt * chainData[0].y);
            scale *= aa;
            // update the thermostat positions
            for (int bead = 0; bead < chainLength; ++bead) {
                chainData[bead].x += 0.5f * chainData[bead].y * wdt;
            }
            // update the forces
            chainForces[0] = (scale * scale * KE2 - numDOFs * kT) / chainMasses[0];
            // update thermostat velocities
            for (int bead = 0; bead < chainLength - 1; ++bead) {
                mixed aa = EXP(-0.125f * wdt * chainData[bead + 1].y);
                chainData[bead].y = aa * (aa * chainData[bead].y + 0.25f * wdt * chainForces[bead]);
                chainForces[bead + 1] = (chainMasses[bead] * chainData[bead].y * chainData[bead].y - kT) / chainMasses[bead + 1];
            }
            chainData[chainLength-1].y += 0.25f * wdt * chainForces[chainLength-1];
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
__kernel void computeHeatBathEnergy(__global mixed* restrict heatBathEnergy, int chainLength, int numDOFs,
                                                 mixed kT, float frequency, __global const mixed2* restrict chainData){
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

__kernel void computeAtomsKineticEnergy(__global mixed2 * restrict energyBuffer, int numAtoms,
                                        __global const mixed4* restrict velm, __global const int *restrict atoms){
    mixed2 energy = (mixed2) (0,0);
    //energy = 1; return;
    int index = get_global_id(0);
    while (index < numAtoms){
        int atom = atoms[index];
        mixed4 v = velm[atom];
        mixed mass = v.w == 0 ? 0 : 1 / v.w;
        energy.x += 0.5f * mass * (v.x*v.x + v.y*v.y + v.z*v.z);
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] = energy;
}

__kernel void computePairsKineticEnergy(__global mixed2 * restrict energyBuffer, int numPairs,
                                        __global const mixed4* restrict velm, __global const int2 *restrict pairs){
    mixed2 energy = (mixed2) (0,0);
    int index = get_global_id(0);
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
        index += get_global_size(0);
    }
    // The atoms version of this has been called already, so accumulate instead of assigning here
    energyBuffer[get_global_id(0)].xy += energy.xy;
}

__kernel void scaleAtomsVelocities(__global mixed2* restrict scaleFactor, int numAtoms,
                                   __global mixed4* restrict velm, __global const int *restrict atoms){
    const mixed scale = scaleFactor[0].x;
    int index = get_global_id(0);
    while (index < numAtoms){
        int atom = atoms[index];
        velm[atom].x *= scale;
        velm[atom].y *= scale;
        velm[atom].z *= scale;
        index += get_global_size(0);
    }
}

__kernel void scalePairsVelocities(__global mixed2 * restrict scaleFactor, int numPairs,
                                   __global mixed4* restrict velm, __global const int2 *restrict pairs){
    int index = get_global_id(0);
    while (index < numPairs){
        int atom1 = pairs[index].x;
        int atom2 = pairs[index].y;
        mixed m1 = velm[atom1].w == 0 ? 0 : 1 / velm[atom1].w;
        mixed m2 = velm[atom2].w == 0 ? 0 : 1 / velm[atom2].w;
        mixed4 cv;
        cv.xyz = (m1*velm[atom1].xyz + m2*velm[atom2].xyz) / (m1 + m2);
        mixed4 rv;
        rv.xyz = velm[atom2].xyz - velm[atom1].xyz;
        velm[atom1].x = scaleFactor[0].x * cv.x - scaleFactor[0].y * rv.x * m2 / (m1 + m2);
        velm[atom1].y = scaleFactor[0].x * cv.y - scaleFactor[0].y * rv.y * m2 / (m1 + m2);
        velm[atom1].z = scaleFactor[0].x * cv.z - scaleFactor[0].y * rv.z * m2 / (m1 + m2);
        velm[atom2].x = scaleFactor[0].x * cv.x + scaleFactor[0].y * rv.x * m1 / (m1 + m2);
        velm[atom2].y = scaleFactor[0].x * cv.y + scaleFactor[0].y * rv.y * m1 / (m1 + m2);
        velm[atom2].z = scaleFactor[0].x * cv.z + scaleFactor[0].y * rv.z * m1 / (m1 + m2);
        index += get_global_size(0);
    }
}

/**
 * Sum the energy buffer containing a pair of energies stored as mixed2.  This is copied from utilities.cu with small modifications
 */
__kernel void reduceEnergyPair(__global const mixed2* restrict energyBuffer, __global mixed2* restrict result, int bufferSize, int workGroupSize, __local mixed2* restrict tempBuffer) {
    const unsigned int thread = get_local_id(0);
    mixed2 sum = (mixed2) (0,0);
    for (unsigned int index = thread; index < bufferSize; index += get_local_size(0)) {
        sum.xy += energyBuffer[index].xy;
    }
    tempBuffer[thread].xy = sum.xy;
    for (int i = 1; i < workGroupSize; i *= 2) {
        barrier(CLK_LOCAL_MEM_FENCE);
        if (thread%(i*2) == 0 && thread+i < workGroupSize) {
            tempBuffer[thread].xy += tempBuffer[thread+i].xy;
        }
    }
    if (thread == 0) {
        *result = tempBuffer[0];
    }
}
