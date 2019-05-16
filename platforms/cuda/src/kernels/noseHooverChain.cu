
#include <initializer_list>

/**
 * Propagate the Nose-Hoover chain with one yoshida-suzuki term
 */
extern "C" __global__ void propagateNoseHooverChain(mixed2* __restrict__ chainData, const mixed * __restrict__ energySum, mixed* __restrict__ scaleFactor,
                                                    mixed* __restrict__ chainMasses, mixed* __restrict__ chainForces, 
                                                    int chainLength, int numMTS, int numDOFs, float timeStep,
                                                    mixed kT, float frequency){
    mixed &scale = scaleFactor[0];
    scale = (mixed) 1;
    const mixed & kineticEnergy = energySum[0];
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
}


/**
 * Compute total (potential + kinetic) energy of the Nose-Hoover beads
 */
extern "C" __global__ void computeHeatBathEnergy(mixed* __restrict__ heatBathEnergy, int chainLength, int numDOFs, 
                                                 mixed kT, float frequency, const mixed2* __restrict__ chainData){

    mixed &energy = heatBathEnergy[0];
    energy = (mixed) 0;

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

extern "C" __global__ void computeMaskedKineticEnergy(mixed * __restrict__ energyBuffer, int paddedNumAtoms,
                                                      const mixed4* __restrict__ velm, const int *__restrict__ mask){
    mixed energy = 0;
    //energy = 1; return;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < paddedNumAtoms; index += blockDim.x*gridDim.x) {
        mixed4 v = velm[index];
        mixed mass = v.w == 0 ? 0 : 1 / v.w;
        if (mask[index] >= 0){
            const mixed4& vparent = velm[mask[index]];
            mixed massp = vparent.w == 0 ? 0 : 1/vparent.w;
            mass = (massp + mass) == 0 ? 0 : (massp * mass) / ( massp + mass );
            v.x -= vparent.x;
            v.y -= vparent.y;
            v.z -= vparent.z;
        }
        energy += 0.5f * mass * (v.x*v.x + v.y*v.y + v.z*v.z);
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] = energy;
}

extern "C" __global__ void scaleVelocities(mixed * __restrict__ scaleFactor, int paddedNumAtoms,
                                           mixed4* __restrict__ velm, const int *__restrict__ mask){
    const mixed &scale = scaleFactor[0];
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < paddedNumAtoms; index += blockDim.x*gridDim.x) {
        mixed4 &v = velm[index];
        mixed4 vparent = mask[index] >= 0 ? velm[mask[index]] : make_mixed4(0.0f, 0.0f, 0.0f, 0.0f);
        mixed maskedScale = mask[index] == index ? 1 : scale;
        v.x = vparent.x + maskedScale * (v.x - vparent.x);
        v.y = vparent.y + maskedScale * (v.y - vparent.y);
        v.z = vparent.z + maskedScale * (v.z - vparent.z);
    }
}
