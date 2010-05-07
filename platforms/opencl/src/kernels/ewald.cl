
float2 multofFloat2(float2 a, float2 b) {
    return (float2) (a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

/**
 * Precompute the cosine and sine sums which appear in each force term.
 */

__kernel void calculateEwaldCosSinSums(__global float* energyBuffer, __global float4* posq, __global float2* cosSinSum, float4 reciprocalPeriodicBoxSize, float reciprocalCoefficient) {
    const unsigned int ksizex = 2*KMAX_X-1;
    const unsigned int ksizey = 2*KMAX_Y-1;
    const unsigned int ksizez = 2*KMAX_Z-1;
    const unsigned int totalK = ksizex*ksizey*ksizez;
    unsigned int index = get_global_id(0);
    float energy = 0.0f;
    while (index < (KMAX_Y-1)*ksizez+KMAX_Z)
        index += get_global_size(0);
    while (index < totalK) {
        // Find the wave vector (kx, ky, kz) this index corresponds to.

        int rx = index/(ksizey*ksizez);
        int remainder = index - rx*ksizey*ksizez;
        int ry = remainder/ksizez;
        int rz = remainder - ry*ksizez - KMAX_Z + 1;
        ry += -KMAX_Y + 1;
        float kx = rx*reciprocalPeriodicBoxSize.x;
        float ky = ry*reciprocalPeriodicBoxSize.y;
        float kz = rz*reciprocalPeriodicBoxSize.z;

        // Compute the sum for this wave vector.

        float2 sum = 0.0f;
        for (int atom = 0; atom < NUM_ATOMS; atom++) {
            float4 apos = posq[atom];
            float phase = apos.x*kx;
            float2 structureFactor = (float2) (cos(phase), sin(phase));
            phase = apos.y*ky;
            structureFactor = multofFloat2(structureFactor, (float2) (cos(phase), sin(phase)));
            phase = apos.z*kz;
            structureFactor = multofFloat2(structureFactor, (float2) (cos(phase), sin(phase)));
            sum += apos.w*structureFactor;
        }
        cosSinSum[index] = sum;

        // Compute the contribution to the energy.

        float k2 = kx*kx + ky*ky + kz*kz;
        float ak = exp(k2*EXP_COEFFICIENT) / k2;
        energy += reciprocalCoefficient*ak*(sum.x*sum.x + sum.y*sum.y);
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy;
}

/**
 * Compute the reciprocal space part of the Ewald force, using the precomputed sums from the
 * previous routine.
 */

__kernel void calculateEwaldForces(__global float4* forceBuffers, __global float4* posq, __global float2* cosSinSum, float4 reciprocalPeriodicBoxSize, float reciprocalCoefficient) {
    unsigned int atom = get_global_id(0);
    while (atom < NUM_ATOMS) {
        float4 force = forceBuffers[atom];
        float4 apos = posq[atom];

        // Loop over all wave vectors.

        int lowry = 0;
        int lowrz = 1;
        for (int rx = 0; rx < KMAX_X; rx++) {
            float kx = rx*reciprocalPeriodicBoxSize.x;
            for (int ry = lowry; ry < KMAX_Y; ry++) {
                float ky = ry*reciprocalPeriodicBoxSize.y;
                float phase = apos.x*kx;
                float2 tab_xy = (float2) (cos(phase), sin(phase));
                phase = apos.y*ky;
                tab_xy = multofFloat2(tab_xy, (float2) (cos(phase), sin(phase)));
                for (int rz = lowrz; rz < KMAX_Z; rz++) {
                    float kz = rz*reciprocalPeriodicBoxSize.z;

                    // Compute the force contribution of this wave vector.

                    int index = rx*(KMAX_Y*2-1)*(KMAX_Z*2-1) + (ry+KMAX_Y-1)*(KMAX_Z*2-1) + (rz+KMAX_Z-1);
                    float k2 = kx*kx + ky*ky + kz*kz;
                    float ak = exp(k2*EXP_COEFFICIENT)/k2;
                    phase = apos.z*kz;
                    float2 structureFactor = multofFloat2(tab_xy, (float2) (cos(phase), sin(phase)));
                    float2 sum = cosSinSum[index];
                    float dEdR = 2*reciprocalCoefficient*ak*apos.w*(sum.x*structureFactor.y - sum.y*structureFactor.x);
                    force.x += dEdR*kx;
                    force.y += dEdR*ky;
                    force.z += dEdR*kz;
                    lowrz = 1 - KMAX_Z;
                }
                lowry = 1 - KMAX_Y;
            }
        }

        // Record the force on the atom.

        forceBuffers[atom] = force;
        atom += get_global_size(0);
    }
}
