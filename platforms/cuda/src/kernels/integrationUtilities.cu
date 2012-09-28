/**
 * Generate random numbers
 */
extern "C" __global__ void generateRandomNumbers(int numValues, float4* __restrict__ random, uint4* __restrict__ seed) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    uint4 state = seed[index];
    unsigned int carry = 0;
    while (index < numValues) {
        float4 value;

        // Generate first value.

        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        unsigned int k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        unsigned int m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x1 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        x1 = sqrt(-2.0f * log(x1));
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x2 = (float)(state.x + state.y + state.w) / (float)0xffffffff;
        value.x = x1 * cos(2.0f * 3.14159265f * x2);

        // Generate second value.

        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x3 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        x3 = sqrt(-2.0f * log(x3));
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x4 = (float)(state.x + state.y + state.w) / (float)0xffffffff;
        value.y = x3 * cos(2.0f * 3.14159265f * x4);

        // Generate third value.

        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x5 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        x5 = sqrt(-2.0f * log(x5));
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x6 = (float)(state.x + state.y + state.w) / (float)0xffffffff;
        value.z = x5 * cos(2.0f * 3.14159265f * x6);

        // Generate fourth value.

        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x7 = (float)max(state.x + state.y + state.w, 0x00000001u) / (float)0xffffffff;
        state.x = state.x * 69069 + 1;
        state.y ^= state.y << 13;
        state.y ^= state.y >> 17;
        state.y ^= state.y << 5;
        x7 = sqrt(-2.0f * log(x7));
        k = (state.z >> 2) + (state.w >> 3) + (carry >> 2);
        m = state.w + state.w + state.z + carry;
        state.z = state.w;
        state.w = m;
        carry = k >> 30;
        float x8 = (float)(state.x + state.y + state.w) / (float)0xffffffff;
        value.w = x7 * cos(2.0f * 3.14159265f * x8);

        // Record the values.

        random[index] = value;
        index += blockDim.x*gridDim.x;
    }
    seed[blockIdx.x*blockDim.x+threadIdx.x] = state;
}

/**
 * Enforce constraints on SHAKE clusters
 */
extern "C" __global__ void applyShakeToPositions(int numClusters, real tol, const real4* __restrict__ oldPos, real4* __restrict__ posDelta, const int4* __restrict__ clusterAtoms, const float4* __restrict__ clusterParams) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 params = clusterParams[index];
        real4 pos = oldPos[atoms.x];
        real4 xpi = posDelta[atoms.x];
        real4 pos1 = oldPos[atoms.y];
        real4 xpj1 = posDelta[atoms.y];
        real4 pos2 = make_real4(0);
        real4 xpj2 = make_real4(0);
        real invMassCentral = params.x;
        real avgMass = params.y;
        float d2 = params.z;
        float invMassPeripheral = params.w;
        if (atoms.z != -1) {
            pos2 = oldPos[atoms.z];
            xpj2 = posDelta[atoms.z];
        }
        real4 pos3 = make_real4(0);
        real4 xpj3 = make_real4(0);
        if (atoms.w != -1) {
            pos3 = oldPos[atoms.w];
            xpj3 = posDelta[atoms.w];
        }

        // Precompute quantities.

        real3 rij1 = make_real3(pos.x-pos1.x, pos.y-pos1.y, pos.z-pos1.z);
        real3 rij2 = make_real3(pos.x-pos2.x, pos.y-pos2.y, pos.z-pos2.z);
        real3 rij3 = make_real3(pos.x-pos3.x, pos.y-pos3.y, pos.z-pos3.z);
        real rij1sq = rij1.x*rij1.x + rij1.y*rij1.y + rij1.z*rij1.z;
        real rij2sq = rij2.x*rij2.x + rij2.y*rij2.y + rij2.z*rij2.z;
        real rij3sq = rij3.x*rij3.x + rij3.y*rij3.y + rij3.z*rij3.z;
        real ld1 = d2-rij1sq;
        real ld2 = d2-rij2sq;
        real ld3 = d2-rij3sq;

        // Iterate until convergence.

        bool converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged) {
            converged = true;
            real3 rpij = make_real3(xpi.x-xpj1.x, xpi.y-xpj1.y, xpi.z-xpj1.z);
            real rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
            real rrpr = rij1.x*rpij.x + rij1.y*rpij.y + rij1.z*rpij.z;
            real diff = fabs(ld1-2.0f*rrpr-rpsqij) / (d2*tol);
            if (diff >= 1.0f) {
                real acor  = (ld1-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij1sq);
                real3 dr = rij1*acor;
                xpi.x += dr.x*invMassCentral;
                xpi.y += dr.y*invMassCentral;
                xpi.z += dr.z*invMassCentral;
                xpj1.x -= dr.x*invMassPeripheral;
                xpj1.y -= dr.y*invMassPeripheral;
                xpj1.z -= dr.z*invMassPeripheral;
                converged = false;
            }
            if (atoms.z != -1) {
                rpij = make_real3(xpi.x-xpj2.x, xpi.y-xpj2.y, xpi.z-xpj2.z);
                rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
                rrpr = rij2.x*rpij.x + rij2.y*rpij.y + rij2.z*rpij.z;
                diff = fabs(ld2-2.0f*rrpr-rpsqij) / (d2*tol);
                if (diff >= 1.0f) {
                    real acor  = (ld2 - 2.0f*rrpr - rpsqij)*avgMass / (rrpr + rij2sq);
                    real3 dr = rij2*acor;
                    xpi.x += dr.x*invMassCentral;
                    xpi.y += dr.y*invMassCentral;
                    xpi.z += dr.z*invMassCentral;
                    xpj2.x -= dr.x*invMassPeripheral;
                    xpj2.y -= dr.y*invMassPeripheral;
                    xpj2.z -= dr.z*invMassPeripheral;
                    converged = false;
                }
            }
            if (atoms.w != -1) {
                rpij = make_real3(xpi.x-xpj3.x, xpi.y-xpj3.y, xpi.z-xpj3.z);
                rpsqij = rpij.x*rpij.x + rpij.y*rpij.y + rpij.z*rpij.z;
                rrpr = rij3.x*rpij.x + rij3.y*rpij.y + rij3.z*rpij.z;
                diff = fabs(ld3 - 2.0f*rrpr - rpsqij) / (d2*tol);
                if (diff >= 1.0f) {
                    real acor  = (ld3-2.0f*rrpr-rpsqij)*avgMass / (rrpr+rij3sq);
                    real3 dr = rij3*acor;
                    xpi.x += dr.x*invMassCentral;
                    xpi.y += dr.y*invMassCentral;
                    xpi.z += dr.z*invMassCentral;
                    xpj3.x -= dr.x*invMassPeripheral;
                    xpj3.y -= dr.y*invMassPeripheral;
                    xpj3.z -= dr.z*invMassPeripheral;
                    converged = false;
                }
            }
            iteration++;
        }

        // Record the new positions.

        posDelta[atoms.x] = xpi;
        posDelta[atoms.y] = xpj1;
        if (atoms.z != -1)
            posDelta[atoms.z] = xpj2;
        if (atoms.w != -1)
            posDelta[atoms.w] = xpj3;
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Enforce velocity constraints on SHAKE clusters
 */
extern "C" __global__ void applyShakeToVelocities(int numClusters, real tol, const real4* __restrict__ oldPos, real4* __restrict__ posDelta, const int4* __restrict__ clusterAtoms, const float4* __restrict__ clusterParams) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float4 params = clusterParams[index];
        real4 pos = oldPos[atoms.x];
        real4 xpi = posDelta[atoms.x];
        real4 pos1 = oldPos[atoms.y];
        real4 xpj1 = posDelta[atoms.y];
        real4 pos2 = make_real4(0);
        real4 xpj2 = make_real4(0);
        real invMassCentral = params.x;
        real avgMass = params.y;
        float d2 = params.z;
        float invMassPeripheral = params.w;
        if (atoms.z != -1) {
            pos2 = oldPos[atoms.z];
            xpj2 = posDelta[atoms.z];
        }
        real4 pos3 = make_real4(0);
        real4 xpj3 = make_real4(0);
        if (atoms.w != -1) {
            pos3 = oldPos[atoms.w];
            xpj3 = posDelta[atoms.w];
        }

        // Precompute quantities.

        real3 rij1 = make_real3(pos.x-pos1.x, pos.y-pos1.y, pos.z-pos1.z);
        real3 rij2 = make_real3(pos.x-pos2.x, pos.y-pos2.y, pos.z-pos2.z);
        real3 rij3 = make_real3(pos.x-pos3.x, pos.y-pos3.y, pos.z-pos3.z);
        real rij1sq = rij1.x*rij1.x + rij1.y*rij1.y + rij1.z*rij1.z;
        real rij2sq = rij2.x*rij2.x + rij2.y*rij2.y + rij2.z*rij2.z;
        real rij3sq = rij3.x*rij3.x + rij3.y*rij3.y + rij3.z*rij3.z;
        real ld1 = d2-rij1sq;
        real ld2 = d2-rij2sq;
        real ld3 = d2-rij3sq;

        // Iterate until convergence.

        bool converged = false;
        int iteration = 0;
        while (iteration < 15 && !converged) {
            converged = true;
            real3 rpij = make_real3(xpi.x-xpj1.x, xpi.y-xpj1.y, xpi.z-xpj1.z);
            real rrpr = rpij.x*rij1.x + rpij.y*rij1.y + rpij.z*rij1.z;
            real delta = -2.0f*avgMass*rrpr/rij1sq;
            real3 dr = rij1*delta;
            xpi.x += dr.x*invMassCentral;
            xpi.y += dr.y*invMassCentral;
            xpi.z += dr.z*invMassCentral;
            xpj1.x -= dr.x*invMassPeripheral;
            xpj1.y -= dr.y*invMassPeripheral;
            xpj1.z -= dr.z*invMassPeripheral;
            if (fabs(delta) > tol)
                converged = false;
            if (atoms.z != -1) {
                rpij = make_real3(xpi.x-xpj2.x, xpi.y-xpj2.y, xpi.z-xpj2.z);
                rrpr = rpij.x*rij2.x + rpij.y*rij2.y + rpij.z*rij2.z;
                delta = -2.0f*avgMass*rrpr/rij2sq;
                dr = rij2*delta;
                xpi.x += dr.x*invMassCentral;
                xpi.y += dr.y*invMassCentral;
                xpi.z += dr.z*invMassCentral;
                xpj2.x -= dr.x*invMassPeripheral;
                xpj2.y -= dr.y*invMassPeripheral;
                xpj2.z -= dr.z*invMassPeripheral;
                if (fabs(delta) > tol)
                    converged = false;
            }
            if (atoms.w != -1) {
                rpij = make_real3(xpi.x-xpj3.x, xpi.y-xpj3.y, xpi.z-xpj3.z);
                rrpr = rpij.x*rij3.x + rpij.y*rij3.y + rpij.z*rij3.z;
                delta = -2.0f*avgMass*rrpr/rij3sq;
                dr = rij3*delta;
                xpi.x += dr.x*invMassCentral;
                xpi.y += dr.y*invMassCentral;
                xpi.z += dr.z*invMassCentral;
                xpj3.x -= dr.x*invMassPeripheral;
                xpj3.y -= dr.y*invMassPeripheral;
                xpj3.z -= dr.z*invMassPeripheral;
                if (fabs(delta) > tol)
                    converged = false;
            }
            iteration++;
        }

        // Record the new positions.

        posDelta[atoms.x] = xpi;
        posDelta[atoms.y] = xpj1;
        if (atoms.z != -1)
            posDelta[atoms.z] = xpj2;
        if (atoms.w != -1)
            posDelta[atoms.w] = xpj3;
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Enforce constraints on SETTLE clusters
 */
extern "C" __global__ void applySettleToPositions(int numClusters, float tol, const real4* __restrict__ oldPos, real4* __restrict__ posDelta, const real4* __restrict__ velm, const int4* __restrict__ clusterAtoms, const float2* __restrict__ clusterParams) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numClusters) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        float2 params = clusterParams[index];
        real4 apos0 = oldPos[atoms.x];
        real4 xp0 = posDelta[atoms.x];
        real4 apos1 = oldPos[atoms.y];
        real4 xp1 = posDelta[atoms.y];
        real4 apos2 = oldPos[atoms.z];
        real4 xp2 = posDelta[atoms.z];
        real m0 = RECIP(velm[atoms.x].w);
        real m1 = RECIP(velm[atoms.y].w);
        real m2 = RECIP(velm[atoms.z].w);

        // Apply the SETTLE algorithm.

        real xb0 = apos1.x-apos0.x;
        real yb0 = apos1.y-apos0.y;
        real zb0 = apos1.z-apos0.z;
        real xc0 = apos2.x-apos0.x;
        real yc0 = apos2.y-apos0.y;
        real zc0 = apos2.z-apos0.z;

        real invTotalMass = RECIP(m0+m1+m2);
        real xcom = (xp0.x*m0 + (xb0+xp1.x)*m1 + (xc0+xp2.x)*m2) * invTotalMass;
        real ycom = (xp0.y*m0 + (yb0+xp1.y)*m1 + (yc0+xp2.y)*m2) * invTotalMass;
        real zcom = (xp0.z*m0 + (zb0+xp1.z)*m1 + (zc0+xp2.z)*m2) * invTotalMass;

        real xa1 = xp0.x - xcom;
        real ya1 = xp0.y - ycom;
        real za1 = xp0.z - zcom;
        real xb1 = xb0 + xp1.x - xcom;
        real yb1 = yb0 + xp1.y - ycom;
        real zb1 = zb0 + xp1.z - zcom;
        real xc1 = xc0 + xp2.x - xcom;
        real yc1 = yc0 + xp2.y - ycom;
        real zc1 = zc0 + xp2.z - zcom;

        real xaksZd = yb0*zc0 - zb0*yc0;
        real yaksZd = zb0*xc0 - xb0*zc0;
        real zaksZd = xb0*yc0 - yb0*xc0;
        real xaksXd = ya1*zaksZd - za1*yaksZd;
        real yaksXd = za1*xaksZd - xa1*zaksZd;
        real zaksXd = xa1*yaksZd - ya1*xaksZd;
        real xaksYd = yaksZd*zaksXd - zaksZd*yaksXd;
        real yaksYd = zaksZd*xaksXd - xaksZd*zaksXd;
        real zaksYd = xaksZd*yaksXd - yaksZd*xaksXd;

        real axlng = SQRT(xaksXd*xaksXd + yaksXd*yaksXd + zaksXd*zaksXd);
        real aylng = SQRT(xaksYd*xaksYd + yaksYd*yaksYd + zaksYd*zaksYd);
        real azlng = SQRT(xaksZd*xaksZd + yaksZd*yaksZd + zaksZd*zaksZd);
        real trns11 = xaksXd / axlng;
        real trns21 = yaksXd / axlng;
        real trns31 = zaksXd / axlng;
        real trns12 = xaksYd / aylng;
        real trns22 = yaksYd / aylng;
        real trns32 = zaksYd / aylng;
        real trns13 = xaksZd / azlng;
        real trns23 = yaksZd / azlng;
        real trns33 = zaksZd / azlng;

        real xb0d = trns11*xb0 + trns21*yb0 + trns31*zb0;
        real yb0d = trns12*xb0 + trns22*yb0 + trns32*zb0;
        real xc0d = trns11*xc0 + trns21*yc0 + trns31*zc0;
        real yc0d = trns12*xc0 + trns22*yc0 + trns32*zc0;
        real za1d = trns13*xa1 + trns23*ya1 + trns33*za1;
        real xb1d = trns11*xb1 + trns21*yb1 + trns31*zb1;
        real yb1d = trns12*xb1 + trns22*yb1 + trns32*zb1;
        real zb1d = trns13*xb1 + trns23*yb1 + trns33*zb1;
        real xc1d = trns11*xc1 + trns21*yc1 + trns31*zc1;
        real yc1d = trns12*xc1 + trns22*yc1 + trns32*zc1;
        real zc1d = trns13*xc1 + trns23*yc1 + trns33*zc1;

        //                                        --- Step2  A2' ---

        float rc = 0.5f*params.y;
        float rb = SQRT(params.x*params.x-rc*rc);
        real ra = rb*(m1+m2)*invTotalMass;
        rb -= ra;
        real sinphi = za1d/ra;
        real cosphi = SQRT(1-sinphi*sinphi);
        real sinpsi = (zb1d-zc1d) / (2*rc*cosphi);
        real cospsi = SQRT(1-sinpsi*sinpsi);

        real ya2d =   ra*cosphi;
        real xb2d = - rc*cospsi;
        real yb2d = - rb*cosphi - rc*sinpsi*sinphi;
        real yc2d = - rb*cosphi + rc*sinpsi*sinphi;
        real xb2d2 = xb2d*xb2d;
        real hh2 = 4.0f*xb2d2 + (yb2d-yc2d)*(yb2d-yc2d) + (zb1d-zc1d)*(zb1d-zc1d);
        real deltx = 2.0f*xb2d + SQRT(4.0f*xb2d2 - hh2 + params.y*params.y);
        xb2d -= deltx*0.5f;

        //                                        --- Step3  al,be,ga ---

        real alpha = (xb2d*(xb0d-xc0d) + yb0d*yb2d + yc0d*yc2d);
        real beta = (xb2d*(yc0d-yb0d) + xb0d*yb2d + xc0d*yc2d);
        real gamma = xb0d*yb1d - xb1d*yb0d + xc0d*yc1d - xc1d*yc0d;

        real al2be2 = alpha*alpha + beta*beta;
        real sintheta = (alpha*gamma - beta*SQRT(al2be2 - gamma*gamma)) / al2be2;

        //                                        --- Step4  A3' ---

        real costheta = SQRT(1-sintheta*sintheta);
        real xa3d = - ya2d*sintheta;
        real ya3d =   ya2d*costheta;
        real za3d = za1d;
        real xb3d =   xb2d*costheta - yb2d*sintheta;
        real yb3d =   xb2d*sintheta + yb2d*costheta;
        real zb3d = zb1d;
        real xc3d = - xb2d*costheta - yc2d*sintheta;
        real yc3d = - xb2d*sintheta + yc2d*costheta;
        real zc3d = zc1d;

        //                                        --- Step5  A3 ---

        real xa3 = trns11*xa3d + trns12*ya3d + trns13*za3d;
        real ya3 = trns21*xa3d + trns22*ya3d + trns23*za3d;
        real za3 = trns31*xa3d + trns32*ya3d + trns33*za3d;
        real xb3 = trns11*xb3d + trns12*yb3d + trns13*zb3d;
        real yb3 = trns21*xb3d + trns22*yb3d + trns23*zb3d;
        real zb3 = trns31*xb3d + trns32*yb3d + trns33*zb3d;
        real xc3 = trns11*xc3d + trns12*yc3d + trns13*zc3d;
        real yc3 = trns21*xc3d + trns22*yc3d + trns23*zc3d;
        real zc3 = trns31*xc3d + trns32*yc3d + trns33*zc3d;

        xp0.x = xcom + xa3;
        xp0.y = ycom + ya3;
        xp0.z = zcom + za3;
        xp1.x = xcom + xb3 - xb0;
        xp1.y = ycom + yb3 - yb0;
        xp1.z = zcom + zb3 - zb0;
        xp2.x = xcom + xc3 - xc0;
        xp2.y = ycom + yc3 - yc0;
        xp2.z = zcom + zc3 - zc0;

        // Record the new positions.

        posDelta[atoms.x] = xp0;
        posDelta[atoms.y] = xp1;
        posDelta[atoms.z] = xp2;
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Enforce velocity constraints on SETTLE clusters
 */
extern "C" __global__ void applySettleToVelocities(int numClusters, float tol, const real4* __restrict__ oldPos, real4* __restrict__ posDelta, real4* __restrict__ velm, const int4* __restrict__ clusterAtoms, const float2* __restrict__ clusterParams) {
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numClusters; index += blockDim.x*gridDim.x) {
        // Load the data for this cluster.

        int4 atoms = clusterAtoms[index];
        real4 apos0 = oldPos[atoms.x];
        real4 apos1 = oldPos[atoms.y];
        real4 apos2 = oldPos[atoms.z];
        real4 v0 = velm[atoms.x];
        real4 v1 = velm[atoms.y];
        real4 v2 = velm[atoms.z];
        
        // Compute intermediate quantities: the atom masses, the bond directions, the relative velocities,
        // and the angle cosines and sines.
        
        real mA = RECIP(v0.w);
        real mB = RECIP(v1.w);
        real mC = RECIP(v2.w);
        real3 eAB = make_real3(apos1.x-apos0.x, apos1.y-apos0.y, apos1.z-apos0.z);
        real3 eBC = make_real3(apos2.x-apos1.x, apos2.y-apos1.y, apos2.z-apos1.z);
        real3 eCA = make_real3(apos0.x-apos2.x, apos0.y-apos2.y, apos0.z-apos2.z);
        eAB *= RECIP(SQRT(eAB.x*eAB.x + eAB.y*eAB.y + eAB.z*eAB.z));
        eBC *= RECIP(SQRT(eBC.x*eBC.x + eBC.y*eBC.y + eBC.z*eBC.z));
        eCA *= RECIP(SQRT(eCA.x*eCA.x + eCA.y*eCA.y + eCA.z*eCA.z));
        real vAB = (v1.x-v0.x)*eAB.x + (v1.y-v0.y)*eAB.y + (v1.z-v0.z)*eAB.z;
        real vBC = (v2.x-v1.x)*eBC.x + (v2.y-v1.y)*eBC.y + (v2.z-v1.z)*eBC.z;
        real vCA = (v0.x-v2.x)*eCA.x + (v0.y-v2.y)*eCA.y + (v0.z-v2.z)*eCA.z;
        real cA = -(eAB.x*eCA.x + eAB.y*eCA.y + eAB.z*eCA.z);
        real cB = -(eAB.x*eBC.x + eAB.y*eBC.y + eAB.z*eBC.z);
        real cC = -(eBC.x*eCA.x + eBC.y*eCA.y + eBC.z*eCA.z);
        real s2A = 1-cA*cA;
        real s2B = 1-cB*cB;
        real s2C = 1-cC*cC;
        
        // Solve the equations.  These are different from those in the SETTLE paper (JCC 13(8), pp. 952-962, 1992), because
        // in going from equations B1 to B2, they make the assumption that mB=mC (but don't bother to mention they're
        // making that assumption).  We allow all three atoms to have different masses.
        
        real mABCinv = RECIP(mA*mB*mC);
        real denom = (((s2A*mB+s2B*mA)*mC+(s2A*mB*mB+2*(cA*cB*cC+1)*mA*mB+s2B*mA*mA))*mC+s2C*mA*mB*(mA+mB))*mABCinv;
        real tab = ((cB*cC*mA-cA*mB-cA*mC)*vCA + (cA*cC*mB-cB*mC-cB*mA)*vBC + (s2C*mA*mA*mB*mB*mABCinv+(mA+mB+mC))*vAB)/denom;
        real tbc = ((cA*cB*mC-cC*mB-cC*mA)*vCA + (s2A*mB*mB*mC*mC*mABCinv+(mA+mB+mC))*vBC + (cA*cC*mB-cB*mA-cB*mC)*vAB)/denom;
        real tca = ((s2B*mA*mA*mC*mC*mABCinv+(mA+mB+mC))*vCA + (cA*cB*mC-cC*mB-cC*mA)*vBC + (cB*cC*mA-cA*mB-cA*mC)*vAB)/denom;
        v0.x += (tab*eAB.x - tca*eCA.x)*v0.w;
        v0.y += (tab*eAB.y - tca*eCA.y)*v0.w;
        v0.z += (tab*eAB.z - tca*eCA.z)*v0.w;
        v1.x += (tbc*eBC.x - tab*eAB.x)*v1.w;
        v1.y += (tbc*eBC.y - tab*eAB.y)*v1.w;
        v1.z += (tbc*eBC.z - tab*eAB.z)*v1.w;
        v2.x += (tca*eCA.x - tbc*eBC.x)*v2.w;
        v2.y += (tca*eCA.y - tbc*eBC.y)*v2.w;
        v2.z += (tca*eCA.z - tbc*eBC.z)*v2.w;
        velm[atoms.x] = v0;
        velm[atoms.y] = v1;
        velm[atoms.z] = v2;
    }
}

/**
 * Compute the direction each CCMA constraint is pointing in.  This is called once at the beginning of constraint evaluation.
 */
extern "C" __global__ void computeCCMAConstraintDirections(const int2* __restrict__ constraintAtoms, real4* __restrict__ constraintDistance, const real4* __restrict__ atomPositions) {
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CCMA_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        // Compute the direction for this constraint.

        int2 atoms = constraintAtoms[index];
        real4 dir = constraintDistance[index];
        real4 oldPos1 = atomPositions[atoms.x];
        real4 oldPos2 = atomPositions[atoms.y];
        dir.x = oldPos1.x-oldPos2.x;
        dir.y = oldPos1.y-oldPos2.y;
        dir.z = oldPos1.z-oldPos2.z;
        constraintDistance[index] = dir;
    }
}

/**
 * Compute the force applied by each CCMA position constraint.
 */
extern "C" __global__ void computeCCMAPositionConstraintForce(const int2* __restrict__ constraintAtoms, const real4* __restrict__ constraintDistance, const real4* __restrict__ atomPositions,
        const real* __restrict__ reducedMass, real* __restrict__ delta1, int* __restrict__ converged, float tol, int iteration) {
    __shared__ int groupConverged;
    if (converged[1-iteration%2]) {
        if (blockIdx.x == 0 && threadIdx.x == 0)
            converged[iteration%2] = 1;
        return; // The constraint iteration has already converged.
    }
    if (threadIdx.x == 0)
        groupConverged = 1;
    __syncthreads();
    real lowerTol = 1.0f-2.0f*tol+tol*tol;
    real upperTol = 1.0f+2.0f*tol+tol*tol;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CCMA_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        real4 dir = constraintDistance[index];
        real4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
        rp_ij.x += dir.x;
        rp_ij.y += dir.y;
        rp_ij.z += dir.z;
        real rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        real d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
        real rp2 = rp_ij.x*rp_ij.x + rp_ij.y*rp_ij.y + rp_ij.z*rp_ij.z;
        real dist2 = dir.w*dir.w;
        real diff = dist2 - rp2;
        delta1[index] = (rrpr > d_ij2*1e-6f ? reducedMass[index]*diff/rrpr : 0.0f);

        // See whether it has converged.

        if (groupConverged && (rp2 < lowerTol*dist2 || rp2 > upperTol*dist2)) {
            groupConverged = 0;
            converged[iteration%2] = 0;
        }
    }
}

/**
 * Compute the force applied by each CCMA velocity constraint.
 */
extern "C" __global__ void computeCCMAVelocityConstraintForce(const int2* __restrict__ constraintAtoms, const real4* __restrict__ constraintDistance, const real4* __restrict__ atomPositions,
        const real* __restrict__ reducedMass, real* __restrict__ delta1, int* __restrict__ converged, float tol, int iteration) {
    __shared__ int groupConverged;
    if (converged[1-iteration%2]) {
        if (blockIdx.x == 0 && threadIdx.x == 0)
            converged[iteration%2] = 1;
        return; // The constraint iteration has already converged.
    }
    if (threadIdx.x == 0)
        groupConverged = 1;
    __syncthreads();
    real lowerTol = 1.0f-2.0f*tol+tol*tol;
    real upperTol = 1.0f+2.0f*tol+tol*tol;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CCMA_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        // Compute the force due to this constraint.

        int2 atoms = constraintAtoms[index];
        real4 dir = constraintDistance[index];
        real4 rp_ij = atomPositions[atoms.x]-atomPositions[atoms.y];
        real rrpr = rp_ij.x*dir.x + rp_ij.y*dir.y + rp_ij.z*dir.z;
        real d_ij2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
        delta1[index] = -2.0f*reducedMass[index]*rrpr/d_ij2;

        // See whether it has converged.

        if (groupConverged && fabs(delta1[index]) > tol) {
            groupConverged = 0;
            converged[iteration%2] = 0;
        }
    }
}

/**
 * Multiply the vector of CCMA constraint forces by the constraint matrix.
 */
extern "C" __global__ void multiplyByCCMAConstraintMatrix(const real* __restrict__ delta1, real* __restrict__ delta2, const int* __restrict__ constraintMatrixColumn,
        const real* __restrict__ constraintMatrixValue, const int* __restrict__ converged, int iteration) {
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.

    // Multiply by the inverse constraint matrix.

    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_CCMA_CONSTRAINTS; index += blockDim.x*gridDim.x) {
        real sum = 0.0f;
        for (int i = 0; ; i++) {
            int element = index+i*NUM_CCMA_CONSTRAINTS;
            int column = constraintMatrixColumn[element];
            if (column >= NUM_CCMA_CONSTRAINTS)
                break;
            sum += delta1[column]*constraintMatrixValue[element];
        }
        delta2[index] = sum;
    }
}

/**
 * Update the atom positions based on CCMA constraint forces.
 */
extern "C" __global__ void updateCCMAAtomPositions(const int* __restrict__ numAtomConstraints, const int* __restrict__ atomConstraints, const real4* __restrict__ constraintDistance,
        real4* __restrict__ atomPositions, const real4* __restrict__ velm, const real* __restrict__ delta1, const real* __restrict__ delta2, int* __restrict__ converged, int iteration) {
    if (blockIdx.x == 0 && threadIdx.x == 0)
        converged[1-iteration%2] = 1;
    if (converged[iteration%2])
        return; // The constraint iteration has already converged.
    real damping = (iteration < 2 ? 0.5f : 1.0f);
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_ATOMS; index += blockDim.x*gridDim.x) {
        // Compute the new position of this atom.

        real4 atomPos = atomPositions[index];
        real invMass = velm[index].w;
        int num = numAtomConstraints[index];
        for (int i = 0; i < num; i++) {
            int constraint = atomConstraints[index+i*NUM_ATOMS];
            bool forward = (constraint > 0);
            constraint = (forward ? constraint-1 : -constraint-1);
            real constraintForce = damping*invMass*delta2[constraint];
            constraintForce = (forward ? constraintForce : -constraintForce);
            real4 dir = constraintDistance[constraint];
            atomPos.x += constraintForce*dir.x;
            atomPos.y += constraintForce*dir.y;
            atomPos.z += constraintForce*dir.z;
        }
        atomPositions[index] = atomPos;
    }
}

/**
 * Compute the positions of virtual sites
 */
extern "C" __global__ void computeVirtualSites(real4* __restrict__ posq, const int4* __restrict__ avg2Atoms, const real2* __restrict__ avg2Weights,
        const int4* __restrict__ avg3Atoms, const real4* __restrict__ avg3Weights,
        const int4* __restrict__ outOfPlaneAtoms, const real4* __restrict__ outOfPlaneWeights) {
    
    // Two particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_2_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        real4 pos = posq[atoms.x];
        real4 pos1 = posq[atoms.y];
        real4 pos2 = posq[atoms.z];
        pos.x = pos1.x*weights.x + pos2.x*weights.y;
        pos.y = pos1.y*weights.x + pos2.y*weights.y;
        pos.z = pos1.z*weights.x + pos2.z*weights.y;
        posq[atoms.x] = pos;
    }
    
    // Three particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_3_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        real4 pos = posq[atoms.x];
        real4 pos1 = posq[atoms.y];
        real4 pos2 = posq[atoms.z];
        real4 pos3 = posq[atoms.w];
        pos.x = pos1.x*weights.x + pos2.x*weights.y + pos3.x*weights.z;
        pos.y = pos1.y*weights.x + pos2.y*weights.y + pos3.y*weights.z;
        pos.z = pos1.z*weights.x + pos2.z*weights.y + pos3.z*weights.z;
        posq[atoms.x] = pos;
    }
    
    // Out of plane sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_OUT_OF_PLANE; index += blockDim.x*gridDim.x) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        real4 pos = posq[atoms.x];
        real4 pos1 = posq[atoms.y];
        real4 pos2 = posq[atoms.z];
        real4 pos3 = posq[atoms.w];
        real4 v12 = pos2-pos1;
        real4 v13 = pos3-pos1;
        real3 cr = cross(v12, v13);
        pos.x = pos1.x + v12.x*weights.x + v13.x*weights.y + cr.x*weights.z;
        pos.y = pos1.y + v12.y*weights.x + v13.y*weights.y + cr.y*weights.z;
        pos.z = pos1.z + v12.z*weights.x + v13.z*weights.y + cr.z*weights.z;
        posq[atoms.x] = pos;
    }
}

inline __device__ real3 loadForce(int index, long long* __restrict__ force) {
    real scale = RECIP((real) 0xFFFFFFFF);
    return make_real3(scale*force[index], scale*force[index+PADDED_NUM_ATOMS], scale*force[index+PADDED_NUM_ATOMS*2]);
}

inline __device__ void storeForce(int index, long long* __restrict__ force, real3 value) {
    force[index] = (long long) (value.x*0xFFFFFFFF);
    force[index+PADDED_NUM_ATOMS] = (long long) (value.y*0xFFFFFFFF);
    force[index+2*PADDED_NUM_ATOMS] = (long long) (value.z*0xFFFFFFFF);
}

/**
 * Distribute forces from virtual sites to the atoms they are based on.
 */
extern "C" __global__ void distributeVirtualSiteForces(const real4* __restrict__ posq, long long* __restrict__ force,
        const int4* __restrict__ avg2Atoms, const real2* __restrict__ avg2Weights,
        const int4* __restrict__ avg3Atoms, const real4* __restrict__ avg3Weights,
        const int4* __restrict__ outOfPlaneAtoms, const real4* __restrict__ outOfPlaneWeights) {
    
    // Two particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_2_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg2Atoms[index];
        real2 weights = avg2Weights[index];
        real3 f = loadForce(atoms.x, force);
        real3 f1 = loadForce(atoms.y, force);
        real3 f2 = loadForce(atoms.z, force);
        f1 += f*weights.x;
        f2 += f*weights.y;
        storeForce(atoms.y, force, f1);
        storeForce(atoms.z, force, f2);
    }
    
    // Three particle average sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_3_AVERAGE; index += blockDim.x*gridDim.x) {
        int4 atoms = avg3Atoms[index];
        real4 weights = avg3Weights[index];
        real3 f = loadForce(atoms.x, force);
        real3 f1 = loadForce(atoms.y, force);
        real3 f2 = loadForce(atoms.z, force);
        real3 f3 = loadForce(atoms.w, force);
        f1 += f*weights.x;
        f2 += f*weights.y;
        f3 += f*weights.z;
        storeForce(atoms.y, force, f1);
        storeForce(atoms.z, force, f2);
        storeForce(atoms.w, force, f3);
    }
    
    // Out of plane sites.
    
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < NUM_OUT_OF_PLANE; index += blockDim.x*gridDim.x) {
        int4 atoms = outOfPlaneAtoms[index];
        real4 weights = outOfPlaneWeights[index];
        real4 pos1 = posq[atoms.y];
        real4 pos2 = posq[atoms.z];
        real4 pos3 = posq[atoms.w];
        real4 v12 = pos2-pos1;
        real4 v13 = pos3-pos1;
        real3 f = loadForce(atoms.x, force);
        real3 f1 = loadForce(atoms.y, force);
        real3 f2 = loadForce(atoms.z, force);
        real3 f3 = loadForce(atoms.w, force);
        real3 fp2 = make_real3(weights.x*f.x - weights.z*v13.z*f.y + weights.z*v13.y*f.z,
                   weights.z*v13.z*f.x + weights.x*f.y - weights.z*v13.x*f.z,
                  -weights.z*v13.y*f.x + weights.z*v13.x*f.y + weights.x*f.z);
        real3 fp3 = make_real3(weights.y*f.x + weights.z*v12.z*f.y - weights.z*v12.y*f.z,
                  -weights.z*v12.z*f.x + weights.y*f.y + weights.z*v12.x*f.z,
                   weights.z*v12.y*f.x - weights.z*v12.x*f.y + weights.y*f.z);
        f1 += f-fp2-fp3;
        f2 += fp2;
        f3 += fp3;
        storeForce(atoms.y, force, f1);
        storeForce(atoms.z, force, f2);
        storeForce(atoms.w, force, f3);
    }
}
