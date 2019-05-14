
#include <initializer_list>

/**
 * Propagate the Nose-Hoover chain with one yoshida-suzuki term
 */
extern "C" __global__ void propagateNoseHooverChain(real* __restrict__ scaleFactor, real* __restrict__ chainMasses, real* __restrict__ chainForces, 
                                                    int chainLength, int numMTS, int numDOFs, float timeStep,
                                                    real kT, float frequency, real kineticEnergy, real2* __restrict__ chainData){

     real &scale = scaleFactor[0];
     scale = (real) 1;
     for (int bead = 0; bead < chainLength; ++bead) chainMasses[bead] = kT / (frequency * frequency);
     chainMasses[0] *= numDOFs;
     real KE2 = 2.0f * kineticEnergy;
     real timeOverMTS = timeStep / numMTS;
     chainForces[0] = (KE2 - numDOFs * kT) / chainMasses[0];
     for (int bead = 0; bead < chainLength - 1; ++bead) {
         chainForces[bead + 1] = (chainMasses[bead] * chainData[bead].y * chainData[bead].y - kT) / chainMasses[bead + 1];
     }
     for (int mts = 0; mts < numMTS; ++mts) {
         BEGIN_YS_LOOP
             real wdt = ys * timeOverMTS;
             chainData[chainLength-1].y += 0.25f * wdt * chainForces[chainLength-1];
             for (int bead = chainLength - 2; bead >= 0; --bead) {
                 real aa = exp(-0.125f * wdt * chainData[bead + 1].y);
                 chainData[bead].y = aa * (chainData[bead].y * aa + 0.25f * wdt * chainForces[bead]);
             }
             // update particle velocities
             real aa = exp(-0.5f * wdt * chainData[0].y);
             scale *= aa;
             // update the thermostat positions
             for (int bead = 0; bead < chainLength; ++bead) {
                 chainData[bead].x += 0.5f * chainData[bead].y * wdt;
             }
             // update the forces
             chainForces[0] = (scale * scale * KE2 - numDOFs * kT) / chainMasses[0];
             // update thermostat velocities
             for (int bead = 0; bead < chainLength - 1; ++bead) {
                 real aa = exp(-0.125f * wdt * chainData[bead + 1].y);
                 chainData[bead].y = aa * (aa * chainData[bead].y + 0.25f * wdt * chainForces[bead]);
                 chainForces[bead + 1] = (chainMasses[bead] * chainData[bead].y * chainData[bead].y - kT) / chainMasses[bead + 1];
             }
             chainData[chainLength-1].y += 0.25f * wdt * chainForces[chainLength-1];
         END_YS_LOOP
     } // MTS loop
     //printf("SCALE %f\n", scale);

}

/**
 * Compute total (potential + kinetic) energy of the Nose-Hoover beads
 */
extern "C" __global__ void computeHeatBathEnergy(real* __restrict__ heatBathEnergy, int chainLength, int numDOFs, 
                                                 real kT, float frequency, real2* __restrict__ chainData){

    real &energy = heatBathEnergy[0];
    energy = (real) 0;

    for(int i = 0; i < chainLength; ++i) {
        real prefac = i ? 1 : numDOFs;
        real mass = prefac * kT / (frequency * frequency);
        real velocity = chainData[i].y; 
        // The kinetic energy of this bead
        energy += 0.5f * mass * velocity * velocity;
        // The potential energy of this bead
        real position = chainData[i].x;
        energy += prefac * kT * position;
    }
}
