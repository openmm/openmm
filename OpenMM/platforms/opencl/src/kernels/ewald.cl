real2 multofReal2(real2 a, real2 b) {
    return (real2) (a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

/**
 * Precompute the cosine and sine sums which appear in each force term.
 */

__kernel void calculateEwaldCosSinSums(__global mixed* restrict energyBuffer, __global const real4* restrict posq, __global real2* restrict cosSinSum, real4 reciprocalPeriodicBoxSize, real reciprocalCoefficient) {
    const unsigned int ksizex = 2*KMAX_X-1;
    const unsigned int ksizey = 2*KMAX_Y-1;
    const unsigned int ksizez = 2*KMAX_Z-1;
    const unsigned int totalK = ksizex*ksizey*ksizez;
    unsigned int index = get_global_id(0);
    mixed energy = 0;
    while (index < (KMAX_Y-1)*ksizez+KMAX_Z)
        index += get_global_size(0);
    while (index < totalK) {
        // Find the wave vector (kx, ky, kz) this index corresponds to.

        int rx = index/(ksizey*ksizez);
        int remainder = index - rx*ksizey*ksizez;
        int ry = remainder/ksizez;
        int rz = remainder - ry*ksizez - KMAX_Z + 1;
        ry += -KMAX_Y + 1;
        real kx = rx*reciprocalPeriodicBoxSize.x;
        real ky = ry*reciprocalPeriodicBoxSize.y;
        real kz = rz*reciprocalPeriodicBoxSize.z;

        // Compute the sum for this wave vector.

        real2 sum = 0.0f;
        for (int atom = 0; atom < NUM_ATOMS; atom++) {
            real4 apos = posq[atom];
            real phase = apos.x*kx;
            real2 structureFactor = (real2) (cos(phase), sin(phase));
            phase = apos.y*ky;
            structureFactor = multofReal2(structureFactor, (real2) (cos(phase), sin(phase)));
            phase = apos.z*kz;
            structureFactor = multofReal2(structureFactor, (real2) (cos(phase), sin(phase)));
            sum += apos.w*structureFactor;
        }
        cosSinSum[index] = sum;

        // Compute the contribution to the energy.

        real k2 = kx*kx + ky*ky + kz*kz;
        real ak = EXP(k2*EXP_COEFFICIENT) / k2;
        energy += reciprocalCoefficient*ak*(sum.x*sum.x + sum.y*sum.y);
        index += get_global_size(0);
    }
    energyBuffer[get_global_id(0)] += energy;
}

/**
 * Compute the reciprocal space part of the Ewald force, using the precomputed sums from the
 * previous routine.
 */

__kernel void calculateEwaldForces(__global real4* restrict forceBuffers, __global const real4* restrict posq, __global const real2* restrict cosSinSum, real4 reciprocalPeriodicBoxSize, real reciprocalCoefficient) {
    unsigned int atom = get_global_id(0);
    while (atom < NUM_ATOMS) {
        real4 force = forceBuffers[atom];
        real4 apos = posq[atom];

        // Loop over all wave vectors.

        int lowry = 0;
        int lowrz = 1;
        for (int rx = 0; rx < KMAX_X; rx++) {
            real kx = rx*reciprocalPeriodicBoxSize.x;
            for (int ry = lowry; ry < KMAX_Y; ry++) {
                real ky = ry*reciprocalPeriodicBoxSize.y;
                real phase = apos.x*kx;
                real2 tab_xy = (real2) (cos(phase), sin(phase));
                phase = apos.y*ky;
                tab_xy = multofReal2(tab_xy, (real2) (cos(phase), sin(phase)));
                for (int rz = lowrz; rz < KMAX_Z; rz++) {
                    real kz = rz*reciprocalPeriodicBoxSize.z;

                    // Compute the force contribution of this wave vector.

                    int index = rx*(KMAX_Y*2-1)*(KMAX_Z*2-1) + (ry+KMAX_Y-1)*(KMAX_Z*2-1) + (rz+KMAX_Z-1);
                    real k2 = kx*kx + ky*ky + kz*kz;
                    real ak = EXP(k2*EXP_COEFFICIENT)/k2;
                    phase = apos.z*kz;
                    real2 structureFactor = multofReal2(tab_xy, (real2) (cos(phase), sin(phase)));
                    real2 sum = cosSinSum[index];
                    real dEdR = 2*reciprocalCoefficient*ak*apos.w*(sum.x*structureFactor.y - sum.y*structureFactor.x);
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
