/**
 * Perform the first step of Brownian integration.
 */

__kernel void integrateBrownianPart1(mixed tauDeltaT, mixed noiseAmplitude, __global const real4* restrict force,
        __global mixed4* restrict posDelta, __global const mixed4* restrict velm, __global const float4* restrict random, unsigned int randomIndex) {
    randomIndex += get_global_id(0);
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        mixed invMass = velm[index].w;
        if (invMass != 0) {
            posDelta[index] = (mixed4) (tauDeltaT*invMass*force[index].x + noiseAmplitude*sqrt(invMass)*random[randomIndex].x,
                                        tauDeltaT*invMass*force[index].y + noiseAmplitude*sqrt(invMass)*random[randomIndex].y,
                                        tauDeltaT*invMass*force[index].z + noiseAmplitude*sqrt(invMass)*random[randomIndex].z, 0);
        }
        randomIndex += get_global_size(0);
    }
}

/**
 * Perform the second step of Brownian integration.
 */

__kernel void integrateBrownianPart2(mixed oneOverDeltaT, __global real4* posq, __global real4* posqCorrection, __global mixed4* velm, __global const mixed4* restrict posDelta) {
    for (int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        if (velm[index].w != 0) {
            mixed4 delta = posDelta[index];
            velm[index].x = oneOverDeltaT*delta.x;
            velm[index].y = oneOverDeltaT*delta.y;
            velm[index].z = oneOverDeltaT*delta.z;
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = (mixed4) (pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
#ifdef USE_MIXED_PRECISION
            posq[index] = (real4) ((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = (real4) (pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
        }
    }
}
