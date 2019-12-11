
#include <initializer_list>

extern "C" __global__ void propagateNoseHooverChain(mixed2* __restrict__ chainData, const mixed2 * __restrict__ energySum, mixed2* __restrict__ scaleFactor,
                                                    mixed* __restrict__ chainMasses, mixed* __restrict__ chainForces, 
                                                    int chainType, int chainLength, int numMTS, int numDOFs, float timeStep,
                                                    mixed kT, float frequency){
    const mixed & kineticEnergy = chainType ? energySum[0].y : energySum[0].x;
    mixed &scale = chainType ? scaleFactor[0].y : scaleFactor[0].x;
    scale = (mixed) 1;
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
                mixed aa = MIXEDEXP(-0.125f * wdt * chainData[bead + 1].y);
                chainData[bead].y = aa * (chainData[bead].y * aa + 0.25f * wdt * chainForces[bead]);
            }
            // update particle velocities
            mixed aa = MIXEDEXP(-0.5f * wdt * chainData[0].y);
            scale *= aa;
            // update the thermostat positions
            for (int bead = 0; bead < chainLength; ++bead) {
                chainData[bead].x += 0.5f * chainData[bead].y * wdt;
            }
            // update the forces
            chainForces[0] = (scale * scale * KE2 - numDOFs * kT) / chainMasses[0];
            // update thermostat velocities
            for (int bead = 0; bead < chainLength - 1; ++bead) {
                mixed aa = MIXEDEXP(-0.125f * wdt * chainData[bead + 1].y);
                chainData[bead].y = aa * (aa * chainData[bead].y + 0.25f * wdt * chainForces[bead]);
                chainForces[bead + 1] = (chainMasses[bead] * chainData[bead].y * chainData[bead].y - kT) / chainMasses[bead + 1];
            }
            chainData[chainLength-1].y += 0.25f * wdt * chainForces[chainLength-1];
        END_YS_LOOP
    } // MTS loop
}


/**
 * Compute total (potential + kinetic) energy of the Nose-Hoover beads
 */
extern "C" __global__ void computeHeatBathEnergy(mixed* __restrict__ heatBathEnergy, int chainLength, int numDOFs,
                                                 mixed kT, float frequency, const mixed2* __restrict__ chainData){
    // Note that this is always incremented; make sure it's zeroed properly before the first call
    mixed &energy = heatBathEnergy[0];

    for(int i = 0; i < chainLength; ++i) {
        mixed prefac = i ? 1 : numDOFs;
        mixed mass = prefac * kT / (frequency * frequency);
        mixed velocity = chainData[i].y; 
        // The kinetic energy of this bead
        energy += 0.5f * mass * velocity * velocity;
        // The potential energy of this bead
        mixed position = chainData[i].x;
        energy += prefac * kT * position;
    }
}

extern "C" __global__ void computeAtomsKineticEnergy(mixed2 * __restrict__ energyBuffer, int numAtoms,
                                                     const mixed4* __restrict__ velm, const int *__restrict__ atoms){
    mixed2 energy = make_mixed2(0,0);
    //energy = 1; return;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        int atom = atoms[index];
        mixed4 v = velm[atom];
        mixed mass = v.w == 0 ? 0 : 1 / v.w;
        energy.x += 0.5f * mass * (v.x*v.x + v.y*v.y + v.z*v.z);
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] = energy;
}

extern "C" __global__ void computePairsKineticEnergy(mixed2 * __restrict__ energyBuffer, int numPairs,
                                                     const mixed4* __restrict__ velm, const int2 *__restrict__ pairs){
    mixed2 energy = make_mixed2(0,0);
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numPairs; index += blockDim.x*gridDim.x) {
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
    }
    // The atoms version of this has been called already, so accumulate instead of assigning here
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x].x += energy.x;
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x].y += energy.y;
}

extern "C" __global__ void scaleAtomsVelocities(mixed2* __restrict__ scaleFactor, int numAtoms,
                                                mixed4* __restrict__ velm, const int *__restrict__ atoms){
    const mixed &scale = scaleFactor[0].x;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        int atom = atoms[index];
        mixed4 &v = velm[atom];
        v.x *= scale;
        v.y *= scale;
        v.z *= scale;
    }
}

extern "C" __global__ void scalePairsVelocities(mixed2 * __restrict__ scaleFactor, int numPairs,
                                                mixed4* __restrict__ velm, const int2 *__restrict__ pairs){
    const mixed &absScale = scaleFactor[0].x;
    const mixed &relScale = scaleFactor[0].y;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numPairs; index += blockDim.x*gridDim.x) {
        int atom1 = pairs[index].x;
        int atom2 = pairs[index].y;
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
        v1.x = absScale * cv.x - relScale * rv.x * m2 / (m1 + m2);
        v1.y = absScale * cv.y - relScale * rv.y * m2 / (m1 + m2);
        v1.z = absScale * cv.z - relScale * rv.z * m2 / (m1 + m2);
        v2.x = absScale * cv.x + relScale * rv.x * m1 / (m1 + m2);
        v2.y = absScale * cv.y + relScale * rv.y * m1 / (m1 + m2);
        v2.z = absScale * cv.z + relScale * rv.z * m1 / (m1 + m2);
        velm[atom1] = v1;
        velm[atom2] = v2;
    }
}

/**
 * Sum the energy buffer containing a pair of energies stored as mixed2.  This is copied from utilities.cu with small modifications
 */
extern "C" __global__ void reduceEnergyPair(const mixed2* __restrict__ energyBuffer, mixed2* __restrict__ result, int bufferSize, int workGroupSize) {
    __shared__ mixed2 tempBuffer[WORK_GROUP_SIZE];
    const unsigned int thread = threadIdx.x;
    mixed2 sum = make_mixed2(0,0);
    for (unsigned int idx = thread; idx < bufferSize; idx += blockDim.x) {
        sum.x += energyBuffer[idx].x;
        sum.y += energyBuffer[idx].y;
    }
    tempBuffer[thread] = sum;
    for (int i = 1; i < workGroupSize; i *= 2) {
        __syncthreads();
        if (thread%(i*2) == 0 && thread+i < workGroupSize) {
            tempBuffer[thread].x += tempBuffer[thread+i].x;
            tempBuffer[thread].y += tempBuffer[thread+i].y;
        }
    }
    if (thread == 0)
        *result = tempBuffer[0];
}
